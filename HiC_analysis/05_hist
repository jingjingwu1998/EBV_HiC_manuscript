# hist_all_in_one

# EV histograph only for all samples
setwd('/Users/80030577/Desktop/HiC_analysis/EV_normalization')
getwd()
load("normalized_ev.100k_multi_samples_06_03_2025.rda")

# List of all NORMALIZED EV objects for easy iteration
normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

# Define sample names based on the normalized data list
sample_names <- names(normalized_ev_lists)

# Define standard human chromosome names in correct order.
chroms <- paste0('chr',c(1:22,'X'))

# --- Data Preparation: Concatenate EV values for each sample ---
# This section combines all chromosome data for each sample into a single vector.
concatenated_normalized_evs <- list()
for (sample_name in sample_names) {
  current_sample_data <- normalized_ev_lists[[sample_name]]
  
  temp_vec <- c()
  for (chr in chroms) {
    # Assuming all chromosomes exist for all samples.
    temp_vec <- c(temp_vec, current_sample_data[[chr]])
  }
  concatenated_normalized_evs[[sample_name]] <- temp_vec
}


my_eigen_data <- concatenated_normalized_evs

# --- 1. Combine all concatenated_normalized_evs across samples ---
all_concatenated_normalized_evs <- unlist(my_eigen_data)

cat("Combined data created. Total number of eigenvector values:", length(all_concatenated_normalized_evs), "\n\n")


# --- 2. Determine Histogram Breaks and Y-axis Limit for the combined data ---

overall_min_value <- min(all_concatenated_normalized_evs)
overall_max_value <- max(all_concatenated_normalized_evs)

combined_hist_breaks <- seq(from = -0.2, to = 0.2, by = 0.008)

# Calculate the histogram (without plotting) to get frequencies
temp_combined_hist <- hist(all_concatenated_normalized_evs,
                           breaks = combined_hist_breaks, plot = FALSE)

# Get the maximum frequency for the Y-axis limit
max_freq_combined <- max(temp_combined_hist$counts)
common_ylim_combined <- c(0, max_freq_combined * 1.1) # Add a 10% buffer

cat("--- Histogram Settings for Combined Data ---\n")
cat(sprintf("X-axis range: [%.3f, %.3f], Bin Width: %.4f\n",
            min(combined_hist_breaks), max(combined_hist_breaks), diff(combined_hist_breaks)[1]))
cat(sprintf("Maximum frequency: %d\n", max_freq_combined))
cat(sprintf("Y-axis range: [%.0f, %.0f]\n\n",
            common_ylim_combined[1], common_ylim_combined[2]))


# --- 3. Plot the Histogram for all_concatenated_normalized_evs ---

# Define the output PDF file
pdf('all_concatenated_eigenvector_histogram.pdf', width = 8, height = 6) # Standard size for a single plot

# Set up plotting parameters (reset to default for a single plot)
par(mfrow = c(1, 1),
    oma = c(0, 0, 0, 0),
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3, 1, 0))

hist(all_concatenated_normalized_evs,
     breaks = combined_hist_breaks,
     xlim = range(combined_hist_breaks),
     ylim = common_ylim_combined,
     main = "Histogram of All Combined Normalized Eigenvector Values",
     xlab = "Normalized Eigenvector Value",
     ylab = "Frequency",
     col = "lightcoral", # Choose a different color for distinction
     border = "black",
     cex.main = 1.2,
     cex.lab = 1.1,
     cex.axis = 1.0
)

# Close the PDF device
dev.off()
