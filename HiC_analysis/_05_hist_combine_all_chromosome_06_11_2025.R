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


# Assign your data to a variable for easier use
my_eigen_data <- concatenated_normalized_evs

# --- Step 1: Determine Histogram Breaks and X-axis Range ---

# Find the smallest and largest values across ALL your samples.
# This helps decide a good starting range for your histogram.
overall_min_value <- min(unlist(my_eigen_data))
overall_max_value <- max(unlist(my_eigen_data))
cat("Overall minimum eigenvector value found:", overall_min_value, "\n")
cat("Overall maximum eigenvector value found:", overall_max_value, "\n\n")

my_hist_breaks <- seq(from = -0.2, to = 0.2, by = 0.008)
hist_xlim <- range(my_hist_breaks)
hist_xlim
cat("--- Your Chosen Histogram Settings ---\n")
cat("X-axis range (hist_xlim): [", hist_xlim[1], ", ", hist_xlim[2], "]\n", sep = "")
cat("Bin width ('by' value):", diff(my_hist_breaks)[1], "\n")
cat("Number of bins:", length(my_hist_breaks) - 1, "\n\n")

# --- Step 2: Calculate Maximum Frequency for Y-axis Range ---

max_frequencies <- c()
for (sample_name in names(my_eigen_data)) {
  # Calculate histogram for each sample without plotting it
  temp_hist <- hist(my_eigen_data[[sample_name]], breaks = my_hist_breaks, plot = FALSE)
  max_frequencies <- c(max_frequencies, max(temp_hist$counts))
}

# Find the single highest frequency among all samples
global_max_frequency <- max(max_frequencies)
common_ylim <- c(0, global_max_frequency * 1.1)
global_max_frequency

cat("Maximum frequency across all samples (for Y-axis):", global_max_frequency, "\n")
cat("Common Y-axis range (ylim): [", common_ylim[1], ", ", common_ylim[2], "]\n\n", sep = "")

# --- Step 3: Plot Histograms for Each Sample ---

par(mfrow = c(2, 3), mar = c(4, 4, 3, 2) + 0.1) # Adjust margins: bottom, left, top, right

# Loop through each sample and create a histogram
for (sample_name in names(my_eigen_data)) {
  hist(my_eigen_data[[sample_name]],
       breaks = my_hist_breaks,      # Use the same breaks for all plots
       xlim = hist_xlim,             # Use the same X-axis range for all plots
       ylim = common_ylim,           # Use the same Y-axis range for all plots
       main = paste("Eigenvector -", sample_name), # Title for each plot
       xlab = "Eigenvector Value",   # X-axis label
       ylab = "Frequency",           # Y-axis label
       col = "skyblue",              # Bar color
       border = "black"              # Bar border color
  )
}
