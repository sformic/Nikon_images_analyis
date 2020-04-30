#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)

### Nuclei Differential Intensity analysis
nuclei_path_ch0 = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/Macro_2_outputs/nuclei_measurements/ch0/"
nuclei_path_ch1 = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/Macro_2_outputs/nuclei_measurements/ch1/"

# getting character vector with all fine names in folder specified in path
all_files <- list.files(path = nuclei_path_ch0)
# function to get character vector of all sample names
get_sample_names <- function (file_name) {
  return(strsplit(file_name, "_2D")[[1]][1])
}
sample_names <- sort(unique(unlist(lapply(all_files, FUN = get_sample_names))))
# getting character vector of condition names
condition_names <- unique(sort(vapply(strsplit(sample_names, "_"), function(v) {v[1]}, FUN.VALUE = as.character(1))))

# function to get - for a specific sample - the vector of plane numbers
get_planes <- function (file_name_2) {
  return(gsub(".csv", "", strsplit(file_name_2, "plane")[[1]][2]))
}

# function to concatenate all measurements' tables belonging to a sample in one table; the plane number is stored in a column
concat_sample_tables <- function (one_sample_name, path) {
  sample_files <- grep(one_sample_name, all_files, value = TRUE)
  planes <- sort(unlist(lapply(sample_files, FUN = get_planes)))
  sample_measurements_list <- lapply(sample_files, FUN = function(f){fread(paste0(path, f))})
  plane_col <- rep(planes, lapply(sample_measurements_list, nrow))
  sample_measurements_table <- do.call("rbind", sample_measurements_list)
  sample_measurements_table$plane <- plane_col
  return(sample_measurements_table)
}

# applying the function to ch0 measurements for all samples, result is a list of table, one element per sample
ch0_meas_list <- lapply(sample_names, FUN = concat_sample_tables, path = nuclei_path_ch0)
names(ch0_meas_list) <- sample_names
# defining measurements' columns
meas_cols <- grep("Label|plane", names(ch0_meas_list[[1]]), value = TRUE, invert = TRUE)
# addying channel info to measurements' columns
ch0_meas_list <- lapply(ch0_meas_list, function(e) {names(e)[names(e) %in% meas_cols] <- paste0(meas_cols, "_DAPI"); return(e)})

# applying the function to ch1 measurements for all samples, result is a list of table, one element per sample
ch1_meas_list <- lapply(sample_names, FUN = concat_sample_tables, path = nuclei_path_ch1)
names(ch1_meas_list) <- sample_names
# addying channel info to measurements' columns and taking only those columns
ch1_meas_list <- lapply(ch1_meas_list, function(e) {colnames(e)[colnames(e) %in% meas_cols] <- paste0(meas_cols, "_ch1"); return(e[, .SD, .SDcols = patterns("_ch1")])})

# binding ch0 and ch1 measurements in same table for each sample
meas_list <- Map("cbind", ch0_meas_list, ch1_meas_list)

# adding sample name in a column for plotting purposes
meas_list <- lapply(1:length(meas_list), function (i) {cbind(meas_list[[i]], sample = rep(names(meas_list)[i], nrow(meas_list[[i]])))})
names(meas_list) <- sample_names

# sanity check correlation plot - area of DAPI and ch1 mask should be the same and Mean intensities of DAPI and ch1 should be quite correlated
library(PerformanceAnalytics)
pdf(file = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/scripts/images/corr_nuclei_ch0_ch1.pdf")
lapply(meas_list, function (t) {chart.Correlation(t[, .(Mean_DAPI, Max_DAPI, Volume_DAPI, Mean_ch1, Max_ch1, Volume_ch1)], main = unique(t$sample))})  
dev.off() ##Note: the 'main' parameter does not work i.e. there are no titles in the pdf

# some plotting on each sample table in order to understand and check relationship between measurements
plot_some_meas <- function (t) {
  pdf(file = paste0("/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/scripts/images/distrib_ch0_ch1_", unique(t$sample),".pdf"), onefile = TRUE)
  # for each nucleus, computing the highest DAPI mean intensity across planes
  t[, max_DAPI_Mean := max(Mean_DAPI), by = .(Label)]
  hist(t$max_DAPI_Mean, xlab = "max DAPI mean intensity", main = "")
  # keeping only nuclei whose max DAPI mean intensity is above the median
  t_f <- t[max_DAPI_Mean > median(max_DAPI_Mean)]
  g1 <- ggplot(data = t_f) +
    geom_line(aes(x = plane, y = Mean_DAPI, group = as.character(Label), colour = as.character(Label))) +
    labs(colour = "Label") +
    ylab("DAPI mean intensity") +
    theme(plot.title = element_text(size = 20), 
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  plot(g1)
  # for each nucleus, selecting only the plane where DAPI shows max mean intensity
  t_f_bestPlaneOnly <- t_f[Mean_DAPI == max_DAPI_Mean]
  hist(t_f_bestPlaneOnly$Mean_ch1, xlab = "ch1 mean intensity at best nuclear plane", main = "", breaks = 11)
  # labeling those nuclei with max mean DAPI intensity above the median
  t_f_bestPlaneOnly[, is_maxDAPI_aboveMedian := ifelse(max_DAPI_Mean > median(max_DAPI_Mean), "well in focus", "less in focus")]
  hist(t_f_bestPlaneOnly[is_maxDAPI_aboveMedian == "well in focus"]$Mean_ch1, xlab = "ch1 mean intensity at best nuclear plane", main = "nuclei well in focus", breaks = 11)
  g2 <- ggplot(data = t_f_bestPlaneOnly) +
    geom_violin(aes(factor(is_maxDAPI_aboveMedian), Mean_ch1)) +
    xlab("")
  plot(g2)
  dev.off()
}
lapply(meas_list, plot_some_meas)

# function to select only the plane where DAPI shows max mean intensity and only those nuclei with max mean DAPI intensity above the median
bestPlaneOnly_f <- function (sample_table) {
  sample_table[, max_DAPI_Mean := max(Mean_DAPI), by = .(Label)]
  sample_table_f <- sample_table[max_DAPI_Mean > median(max_DAPI_Mean)]
  return(sample_table_f[Mean_DAPI == max_DAPI_Mean])
}

# applying the function to all samples
bestPlaneOnly_meas_list <- lapply(meas_list, bestPlaneOnly_f)

# splitting the list by condition
c1_list <- bestPlaneOnly_meas_list[grepl(paste0("^", condition_names[1], "_"), names(bestPlaneOnly_meas_list))]
names(c1_list) <- gsub("^._", "", names(c1_list))
c2_list <- bestPlaneOnly_meas_list[grepl(paste0("^", condition_names[2], "_"), names(bestPlaneOnly_meas_list))]
names(c2_list) <- gsub("^._", "", names(c2_list))

# function to compute t-test between ch1 mean intensity
#diffInt_test <- function (c1_t, c2_t) {
#  t.test(x = c1_t$Mean_ch1, y = c2_t$Mean_ch1,
#         alternative = "two.sided",
#         mu = 0, paired = FALSE, var.equal = FALSE,
#         conf.level = 0.99)
#}
#DI_test_list <- Map(diffInt_test, c1_list, c2_list)

### Plotting ch1 intensity measurements' distribution differences and significance

# joining conditions for same antibody for plotting purpose
join_cond <- function (c1_t, c2_t) {
  # sampling the same number of nuclei from both conditions of same antibody 
  min_n <- min(nrow(c1_t), nrow(c2_t))
  c1_t <- c1_t[sample(.N,min_n)]
  c2_t <- c2_t[sample(.N,min_n)]
  return(data.frame(Mean_ch1 = c(c1_t$Mean_ch1, c2_t$Mean_ch1), antibody = gsub("^._|_{1,}.*$", "", c(c1_t$sample, c2_t$sample)), condition = rep(condition_names, c(nrow(c1_t), nrow(c2_t)))))
}
DI_list <- Map(join_cond, c1_list, c2_list)
DI_table <- do.call("rbind", DI_list)

# table of unpaired t-test results for each antibody
ttest_result <- compare_means(Mean_ch1 ~ condition, data = DI_table, method = "t.test", paired = FALSE,
                              group.by = "antibody",
                              alternative = "two.sided", mu = 0, var.equal = FALSE)
# changing names of conditions for plotting purposes
DI_table$condition <- factor(DI_table$condition, levels=c("8", "5"))
DI_table$antibody <- gsub("MBS", "MBS8245799", gsub("OGA", "anti-OGA", DI_table$antibody))
DI_table$antibody <- factor(DI_table$antibody, levels = c("ab177941", "ab96718", "MBS8245799", "HGAC85", "anti-OGA"))
v <- ggviolin(DI_table, x = "condition", y = "Mean_ch1",
              color = "condition",
              add = "jitter",
              ylab = "ch1 mean intensity",
              xlab = "",
              facet.by = "antibody", short.panel.labs = TRUE, panel.labs.font = list(size = 16)) +
      ggtitle("Nuclei\n") +
      scale_x_discrete(labels=c("8" = "wt", "5" = expression(paste("NLS", ""^0, ""^"/", ""^-{})))) +
      theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        #legend.text = element_text(size = 19),
        #legend.title = element_text(size = 19),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA),
        legend.key.size = unit(1, "cm"))

# Use only p.format as label. Remove method name.
pdf(file = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/scripts/images/ch1_nuclei_comp_plot.pdf")
v + stat_compare_means(method = "t.test", paired = FALSE, method.args = list(alternative = "two.sided", mu = 0, var.equal = FALSE),
                       label =  "p.format", label.x = 1, label.y = ggplot_build(v)$layout$panel_params[[5]]$y.range[2], size = 5)
dev.off()

### Cytosol Differential Intensity Analysis

cyto_path <- "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/Macro_2_outputs/cyto_measurements/"

# function to read - for a specific sample - the table of ch1 measurements in the cytosol
get_cyto_sample_tables <- function (one_sample_name, path) {
  fread(paste0(path, one_sample_name, ".csv"))
}
# applying the function to all samples
cyto_meas_df <- as.data.frame(t(sapply(sample_names, get_cyto_sample_tables, path = cyto_path)))

# plotting ch1 mean for the different antibodies, comparing conditions
cyto_meas_df$condition <- vapply(strsplit(rownames(cyto_meas_df), "_"), function(v) {v[1]}, FUN.VALUE = as.character(1))
cyto_meas_df$condition <- factor(cyto_meas_df$condition, levels = c("8", "5"))
cyto_meas_df$antibody <- vapply(strsplit(rownames(cyto_meas_df), "_"), function(v) {v[2]}, FUN.VALUE = as.character(1))
cyto_meas_df$antibody <- gsub("MBS", "MBS8245799", gsub("OGA", "anti-OGA", cyto_meas_df$antibody))
cyto_meas_df$antibody <- factor(cyto_meas_df$antibody, levels = c("ab177941", "ab96718", "MBS8245799", "HGAC85", "anti-OGA"))

pdf(file = "/home/sara/PhD/OGT_NLS_project/Nikon_images_analysis/scripts/images/ch1_cyto_comp_plot.pdf")
ggplot(data = cyto_meas_df, aes(x = antibody, y = Mean, fill = condition)) +
  geom_col() +
  scale_fill_discrete(labels=c("8" = "wt", "5" = expression(paste("NLS", ""^0, ""^"/", ""^-{})))) +
  ylab("ch1 mean intensity") +
  xlab("") +
  ggtitle("Cytosol\n") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        legend.text = element_text(size = 19),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill=NA),
        legend.key.size = unit(1, "cm"))
dev.off()