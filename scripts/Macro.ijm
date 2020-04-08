my_path = "Z:/Sara/Julia/Data/confocal/191216 Sara-Julia Nikon A1/";
setBatchMode(true);
list_of_files = getFileList(my_path);
Array.show(list_of_files);

for (i = 0; i < list_of_files.length; i++) {
	print(i);
	name_radix = replace(list_of_files[i], ".nd2", "");
	run("Bio-Formats Importer", "open=[" + my_path + list_of_files[i] + "] autoscale color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
	selectImage(my_path + name_radix + ".nd2 - C=0");
	run("Duplicate...", "title=intensity duplicate");
	selectImage(my_path + name_radix + ".nd2 - C=0");
	run("Duplicate...", "title=blur duplicate");
	run("Gaussian Blur...", "sigma=2.70 scaled stack");
	run("Duplicate...", "title=binary duplicate");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Dark calculate black");
	run("Distance Transform Watershed 3D", "distances=[Quasi-Euclidean (1,1.41,1.73)] output=[16 bits] normalize dynamic=1 connectivity=26");
	run("Duplicate...", "title=labelMask duplicate");
	close("binarydist-watershed");
	close("binary");
	close("blur");
	close(my_path + name_radix + ".nd2 - C=0");
	selectImage("intensity");
	run("Stack to Images");
	//images are called e.g. "c:1/2 z:[1-12]/12 - 5 HGAC85.nd2 (series 1)"
	selectWindow("labelMask");
	run("Stack to Images");
	close("labelMask");
	close("intensity");
	//images are called "labelMask-00[01-12]"
	intensity_titles = newArray(0);
	labelMask_titles = newArray(0);
	/// aim of the next two for loops: for each plane of the stack, measuring DAPI fluorescence intensity in the corresponding plane's DAPI labels	
	for (i=1; i<=nImages(); i++) {
		selectImage(i);
		if(startsWith(getTitle(), "c:")) {
			intensity_titles = Array.concat(intensity_titles, getTitle());
		} else if (startsWith(getTitle(), "labelMask")) {
			labelMask_titles = Array.concat(labelMask_titles, getTitle());
		}
	}
	Array.show(intensity_titles);
	Array.show(labelMask_titles);
	for (i = 0; i < intensity_titles.length; i++) {
		plane_number = replace(replace(intensity_titles[i], "^.*z:", ""), "/.*", "");
		run("Intensity Measurements 2D/3D", "input=" + intensity_titles[i] + " labels=" + labelMask_titles[i] + " mean stddev max median mode numberofvoxels volume");
		saveAs("Results", "Z:/Sara/Julia/Data/confocal/Nikon_A1_image_analysis/nuclei_measurements/" + name_radix + "_2Dmeasurements_plane" + plane_number + ".csv");
		close(name_radix + "_2Dmeasurements_plane" + plane_number + ".csv");
	}
	
	close("*");
}
