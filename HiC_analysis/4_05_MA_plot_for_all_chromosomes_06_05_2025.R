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

concatenated_normalized_evs <- list()
for (sample_name in sample_names) {
  current_sample_data <- normalized_ev_lists[[sample_name]]
  temp_vec <- c()
  for (chr in chroms) {
    temp_vec <- c(temp_vec, current_sample_data[[chr]])
  }
  concatenated_normalized_evs[[sample_name]] <- temp_vec
}
  concatenated_normalized_evs[[sample_name]] <- do.call(c, sample_ev_vectors)
}


# --- Matrix Plot Generation ---

# Define the PDF file name
pdf_file_name <- 'ev.ma_pairwise_plot.pdf'

# Define fixed global axis limits as requested
x_lim_matrix <- c(-0.1, 0.1)
y_lim_matrix <- c(-0.15, 0.15)

# Define plot layout parameters based on the number of samples
num_samples_N <- length(sample_names)
plot_rows <- num_samples_N
plot_cols <- num_samples_N

# Open PDF device with dimensions adjusted for N x N grid
pdf(file = pdf_file_name,
    width = plot_cols * 5, # e.g., 6 columns * 5 inches/column = 30 inches wide
    height = plot_rows * 5) # e.g., 6 rows * 5 inches/row = 30 inches tall

# Set general plot parameters for the entire matrix
par(mfrow = c(plot_rows, plot_cols), # Arrange plots in N rows by N columns (e.g., 6x6 grid)
    font.lab = 2, cex.lab = 1.0, # Bold labels, size 1.0
    mar = c(3, 3, 2, 1), # Margins for each individual plot: bottom, left, top, right (in lines)
    oma = c(2, 2, 2, 2), # Outer margins for the entire figure (for common titles/labels)
    mgp = c(1.5, 0.5, 0), # Axis title and label positions relative to axis line
    xaxs = 'i', yaxs = 'i' # Ensures plotting region extends exactly to axis limits
)

# Loop through all possible combinations for an N x N matrix plot
# Outer loop for rows of the matrix (representing Sample 1)
for (i in 1:num_samples_N) {
  s_row_name <- sample_names[i]
  data_row <- concatenated_normalized_evs[[s_row_name]]
  
  # Inner loop for columns of the matrix (representing Sample 2)
  for (j in 1:num_samples_N) {
    s_col_name <- sample_names[j]
    data_col <- concatenated_normalized_evs[[s_col_name]]
    
    # Check if samples are the same (diagonal)
    if (s_row_name == s_col_name) {
      # For diagonal plots (sample vs. itself), plot a placeholder or just the sample name
      plot(NA, xlim = x_lim_matrix, ylim = y_lim_matrix, type = "n", xlab = "", ylab = "",
           main = paste0(s_row_name), cex.main = 1.5, font.main = 2) # Just the sample name
      text(mean(x_lim_matrix), mean(y_lim_matrix), labels = s_row_name, cex = 3, col = "gray") # Larger sample name
    } else {
      # Calculate Difference (M) and Mean (A) for off-diagonal plots
      # M = data_col - data_row
      # A = (data_row + data_col) / 2
      M_val <- data_col - data_row
      A_val <- (data_row + data_col) / 2
      
      # Plot the MA plot for the current pair
      plot(A_val, M_val,
           pch = 16, cex = 0.25, # Point type and size for dense scatter plots
           col = rgb(0, 0, 0, alpha = 0.1), # Use semi-transparent black for point color to show density
           main = paste0(s_col_name, " vs ", s_row_name), # Title for individual plot
           xlab = "Mean EV", ylab = "EV Difference", # Axis labels for individual plot
           xlim = x_lim_matrix, # Use fixed global x-limits
           ylim = y_lim_matrix, # Use fixed global y-limits
           cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8) # Adjust font sizes within plot
      
      abline(h = 0, lty = "dashed", col = "red") # Horizontal line at EV Difference = 0
    }
  }
}

# Optional: Add an overall title to the entire matrix plot
mtext("Pairwise MA Plots for Normalized EV Data", side = 3, line = 0, outer = TRUE, cex = 2, font = 2)

dev.off() # Close the PDF device
cat(paste0("\nFull pairwise MA plot matrix generated and saved to '", pdf_file_name, "'.\n"))
