#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

my_path = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/nuclei_measurements/"

# get character vector with all fine names in folder specified in my_path
all_files <- list.files(path = my_path)
# function to get character vecotr of all sample names
get_sample_names <- function (file_name) {
  return(strsplit(file_name, "_")[[1]][1])
}
sample_names <- unique(unlist(lapply(all_files, FUN = get_sample_names)))

# function to get - for a specific sample - the vector of plane numbers
get_planes <- function (file_name_2) {
  return(gsub(".csv", "", strsplit(file_name_2, "plane")[[1]][2]))
}

# function to concatenate all measurements' tables belonging to a sample in one table; the plane number is stored in a column
concat_sample_tables <- function (one_sample_name) {
  sample_files <- grep(one_sample_name, all_files, value = TRUE)
  planes <- sort(unlist(lapply(sample_files, FUN = get_planes)))
  sample_measurements_list <- lapply(sample_files, FUN = function(f){fread(paste0(my_path, f))})
  names(sample_measurements_list) <- planes
  plane_col <- rep(planes, lapply(sample_measurements_list, nrow))
  sample_measurements_table <- do.call("rbind", sample_measurements_list)
  sample_measurements_table$plane <- plane_col
  return(sample_measurements_table)
}

# applying the function to all samples, result is a list of table, one element per sample
all_samples_meas_list <- lapply(sample_names, FUN = concat_sample_tables)

# trying out some plotting on one element of the table in order to understand relationship between measurements
ciao <- all_samples_meas_list[[1]]
nuclei_dim <- ciao[, .(maxVol = max(Volume)), by = .(Label)]
summary(nuclei_dim$maxVol)
bigger_nuclei <- nuclei_dim[maxVol > quantile(maxVol, 0.75)]$Label

ciao <- ciao[order(Label, Mean)]
#ciao[order(Label, Mean), "real_plane" := seq(1, .N), by = Label]

to_plot <- ciao[Label %in% sample(bigger_nuclei, 10)]

library(PerformanceAnalytics)
chart.Correlation(to_plot[, .(Mean, Max, Median, Mode, StdDev, Volume)])

### This two plots show the problem: since the mean is computed inside the "plane piece" of the 3D nucleus and this peace is smaller where the nucleus is less visible, the mean changes too.
ggplot(data = ciao[Label == 41]) +
  geom_line(aes(x = plane, y = Volume, group = as.character(Label), colour = as.character(Label))) +
  labs(colour = "Label")
ggplot(data = ciao[Label == 41]) +
  geom_line(aes(x = plane, y = Mean, group = as.character(Label), colour = as.character(Label))) +
  labs(colour = "Label")




