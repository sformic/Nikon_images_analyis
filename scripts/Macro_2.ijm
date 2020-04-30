my_path = "Y:/Sara/Julia/Data/confocal/191216 Sara-Julia Nikon A1/";
setBatchMode(true);
list_of_files = getFileList(my_path);
print(list_of_files.length);

n = 0;
while (n<list_of_files.length) {
	print(n); /// for debugging
	name_radix = replace(list_of_files[n], ".nd2", "");
	name_radix_new = replace(name_radix, " ", "_"); /// necessary for later Bioformat export command
	//// importing image
	run("Bio-Formats Importer", "open=[" + my_path + list_of_files[n] + "] autoscale color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
	selectImage(my_path + name_radix + ".nd2 - C=0");
	run("Duplicate...", "title=intensity_ch0 duplicate");
	selectImage(my_path + name_radix + ".nd2 - C=1");
	run("Duplicate...", "title=intensity_ch1 duplicate");
	close(my_path + name_radix + ".nd2 - C=0");
	close(my_path + name_radix + ".nd2 - C=1");
	close(my_path + name_radix + ".nd2 - C=2");
	
	//// creating DAPI BINARY MASK STACK 1, which is needed for multiplication with ch1 intensity images in order to remove the cytosol for successive quantification
	selectImage("intensity_ch0");
	run("Duplicate...", "title=blur duplicate");
	run("Gaussian Blur...", "sigma=10 stack");
	run("Duplicate...", "title=binary duplicate");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	run("Convert to Mask", "method=Default background=Dark calculate black");
	close("blur");
	/// applying erosion to each slice of DAPI binary mask stack 1 and reconverting it to a stack after erosion. This step is needed to exclude nuclear membrane and cytosol immediately adjacent which get highly stained in channel 1 and can counfound ch1 intensity measurement in the nuclei
	selectWindow("binary");
	run("Stack to Images");
	close("binary");
	binaryMask_titles = newArray(0);
	erodedBinaryMask_titles = newArray(0);
	for (i=1; i<=nImages(); i++) {
		selectImage(i);
		if(startsWith(getTitle(), "c:")) {
			rename("binary_" + (i-2)); /// (i - 2) needed because there are also the two intensity images open
			binaryMask_titles = Array.concat(binaryMask_titles, getTitle());
			run("Morphological Filters", "operation=Erosion element=Square radius=5");
		}
	}
	for (i=1; i<=nImages(); i++) { /// for debugging
		selectImage(i);
		if (endsWith(getTitle(), "Erosion")) {
			erodedBinaryMask_titles = Array.concat(erodedBinaryMask_titles, getTitle());
		}
	}
	Array.show(binaryMask_titles); /// for debugging
	Array.show(erodedBinaryMask_titles); /// for debugging
	run("Images to Stack", "name=erodedBMstack title=Erosion use");
	close("binary_*");

	//// creating a DAPI LABEL MASK 2D from the max z-projection of a second DAPI BINARY MASK STACK (DAPI BINARY MASK STACK 2), because this represents the max area of each nucleus across planes, which is needed it to calculate DAPI mean intensity inside it and thus be able to say which plane is the one where each nucleus shows max DAPI intensity ("best nuclear plane" for that cell)'
	/// creating DAPI BINARY MASK STACK 2
	selectWindow("intensity_ch0");
	run("Duplicate...", "title=blur2 duplicate");
	run("Gaussian Blur...", "sigma=10 stack");
	run("Duplicate...", "title=binary2 duplicate");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	run("Convert to Mask", "method=RenyiEntropy background=Dark calculate black"); /// this time using RenyiEntropy algorithm because no erosion will follow and this algithm seems better at focusing on DAPI area only (although not perfect, still it includes part of cytosol immediately adyacent to nuclear membrane, but it is not so important in this case since we this mask will be used on the DAPI channel)  
	/// MAX z-projection
	run("Z Project...", "projection=[Max Intensity]");
	/// creating the LABEL MASK from the z-projection ('MAX_labelMask')
	run("Distance Transform Watershed 3D", "distances=[Quasi-Euclidean (1,1.41,1.73)] output=[16 bits] normalize dynamic=1 connectivity=26");
	run("Duplicate...", "title=MAX_labelMask duplicate");
	/// saving the DAPI MAX_labelMask for an eventual visual inspetion
	run("Bio-Formats Exporter", "save=Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/images/" + name_radix_new + "_MAX_labelMask.tif " + "export compression=Uncompressed");
	close("MAX_binary2dist-watershed");
	close("binary2");
	close("MAX_binary2");

	//// measuring DAPI intensity in original DAPI intensity image in each plane using as label mask DAPI MAX_labelMask
	selectImage("intensity_ch0");
	run("Stack to Images"); /// output images are called e.g. "c:1/2 z:[1-12]/12 - 5 HGAC85.nd2 (series 1)"
	close("intensity_ch0");
	all_titles = newArray(0);
	intensity_ch0_titles = newArray(0);
	for (i=1; i<=nImages(); i++) {
		selectImage(i);
		all_titles = Array.concat(all_titles, getTitle());
		if(startsWith(getTitle(), "c:")) {
			plane_number = replace(replace(getTitle(), "^.*z:", ""), "/.*", "");
			new_name = "ich0_" + plane_number; /// needed because otherwise the imageTitles are too long and the stupid ImageJ gets confused when it does the intensity measurements 
			rename(new_name);
			intensity_ch0_titles = Array.concat(intensity_ch0_titles, new_name);
		}
	}
	Array.show(intensity_ch0_titles); /// for debugging
	Array.show(all_titles); /// for debugging
	for (i = 0; i<intensity_ch0_titles.length; i++) {
		plane_number = replace(intensity_ch0_titles[i], "^.*_", "");
		run("Intensity Measurements 2D/3D", "input=" + intensity_ch0_titles[i] + " labels=" + "MAX_labelMask" + " mean stddev max min median mode numberofvoxels volume");
		saveAs("Results", "Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/nuclei_measurements/ch0/" + name_radix_new + "_2Dmeasurements_plane" + plane_number + ".csv");			close(name_radix_new + "_2Dmeasurements_plane" + plane_number + ".csv");
	}
	close("c:*");
	close("ich0_*");
	
	//// multiplying DAPI BINARY STACK 1 for the original intensity image of channel 1
	selectWindow("erodedBMstack");
	run("Subtract...", "value=254 stack"); /// necessary because the binary stack has values [0-255] and not already [0-1]
	imageCalculator("Multiply create stack", "intensity_ch1","erodedBMstack");
	close("erodedBMstack");
	selectWindow("Result of intensity_ch1");
	rename("ch1_nucleiOnly");
	/// saving the resunting image for any eventual visual inspection
	run("Bio-Formats Exporter", "save=Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/images/" + name_radix_new + "_ch1_nucleiOnly.tif " + "export compression=Uncompressed");
	
	//// using this ch1_nucleiOnly stack for ch1 intensity measurements inside all DAPI MAX_labelMask nuclei areas and then a posteriori select only the values corresponding to the best nuclear planes in R script 'ch1_at_best_nuclear_plane.R'
	run("Stack to Images"); /// output images are called e.g. "c:2/2 z:[1-12]/12 - 5 HGAC85.nd2 (series 1)"
	all_titles_2 = newArray(0);
	ch1_nucleiOnly_titles = newArray(0);
	for (i=1; i<=nImages(); i++) {
		selectImage(i);
		all_titles_2 = Array.concat(all_titles_2, getTitle());
		if(startsWith(getTitle(), "c:")) {
			plane_number = replace(replace(getTitle(), "^.*z:", ""), "/.*", "");
			new_name = "ich1_" + plane_number; /// needed because otherwise the imageTitles are too long and the stupid ImageJ gets confused when it does the intensity measurements 
			rename(new_name);
			ch1_nucleiOnly_titles = Array.concat(ch1_nucleiOnly_titles, new_name);
		}
	}
	Array.show(ch1_nucleiOnly_titles); /// for debugging
	Array.show(all_titles_2); /// for debugging
	for (i = 0; i<ch1_nucleiOnly_titles.length; i++) {
			plane_number = replace(ch1_nucleiOnly_titles[i],"^.*_", "");
			run("Intensity Measurements 2D/3D", "input=" + ch1_nucleiOnly_titles[i] + " labels=" + "MAX_labelMask" + " mean stddev max min median mode numberofvoxels volume");
			saveAs("Results", "Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/nuclei_measurements/ch1/" + name_radix_new + "_2Dmeasurements_plane" + plane_number + ".csv");
			close(name_radix_new + "_2Dmeasurements_plane" + plane_number + ".csv");
	}
	close("MAX_labelMask");
	close("c:*");
	close("ich1_*");

	//// creating DAPI BINARY STACK 3 and multiply it with ch1 intensity images in order to measure ch1 intensity in the cytosol only - this is needed to control for differential overall staining in ch1
	selectWindow("blur2");
	setAutoThreshold("Huang");
	//run("Threshold...");
	run("Convert to Mask", "method=Huang background=Light calculate"); /// this time Huang algorithm is better in order to obtain cytosol only with no nucleus or nuclear membrane in any of the plane. This time background of binary has to be 'light' because I need the complementary of the mask as outcome of the multiplication
	run("Invert LUT");
	run("Subtract...", "value=254 stack");
	imageCalculator("Multiply create stack", "intensity_ch1","blur2");
	close("blur2");
	close("intensity_ch1");
	selectWindow("Result of intensity_ch1");
	/// saving the resunting image for any eventual visual inspection
	run("Bio-Formats Exporter", "save=Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/images/" + name_radix_new + "_ch1_cytoOnly.tif " + "export compression=Uncompressed");
	/// measuring the max (for each pixel among all planes) mean ch1 intensity over the whole resulting image. This is enough for the purpose of controlling for overall ch1 intensity
	run("Z Project...", "projection=[Max Intensity]");
	run("Measure");
	saveAs("Results", "Y:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/Macro_2_outputs/cyto_measurements/" + name_radix_new + ".csv");
	close("Results");
	close("*");
	n = n + 1;
	print(n);
}