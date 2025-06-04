#-----------------------------------
# This script is used to generate a comprehensive MA plot matrix
# for normalized HiC Eigenvector (EV) data at 100k resolution across 6 samples.
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges) # For runmean, if used in original EV calculation, though not directly in this script.

# output directory
out_dir <- 'normalized_ev/MA_plot_for_all_chromosomes' 
# Load normalized EV data
load(file.path(out_dir, '../normalized_ev.100k_multi_samples_06_03_2025.rda'))

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


# --- 1. Concatenate All Chromosomes per Sample into a Single Long Vector ---
# This creates a flattened representation of each sample's EV data.
# The loop explicitly uses the 'chroms' vector to ensure correct order of concatenation across bins.
concatenated_normalized_evs <- list()

for (sample_name in names(normalized_ev_lists)) {
  sample_ev_vectors <- list()
  
  # Loop through 'chroms' ensures correct order for concatenating chromosome data
  for (chr_name in chroms) { 
    ev_vec <- normalized_ev_lists[[sample_name]][[chr_name]] # Get normalized EV vector
    
    if (!is.null(ev_vec)) {
      # Add to list, maintaining order as defined by the 'chroms' loop
      sample_ev_vectors[[chr_name]] <- ev_vec 
    }
  }
  # Concatenate vectors into a single long vector for the current sample
  concatenated_normalized_evs[[sample_name]] <- do.call(c, sample_ev_vectors)
}

# --- 2. Save the Concatenated Data to a Single File (.rda) ---
# This addresses your request to have "one file" containing all the concatenated data.
save(concatenated_normalized_evs, file = file.path(out_dir, "concatenated_normalized_evs_100k_multi_samples.rda"))
cat("Concatenated normalized EV data saved to:", file.path(out_dir, "concatenated_normalized_evs_100k_multi_samples.rda"), "\n\n")


# --- Removed: Calculation of optimal global axis limits ---
# Define fixed limits for X and Y axes.
x_lim_matrix <- c(-0.1, 0.1)  # Default fixed range for Mean EV (A)
y_lim_matrix <- c(-0.2, 0.2)  # Default fixed range for EV Difference (M)
cat("Fixed MA Plot Matrix X-axis limits (Mean EV): [", round(x_lim_matrix[1], 4), ", ", round(x_lim_matrix[2], 4), "]\n", sep="")
cat("Fixed MA Plot Matrix Y-axis limits (EV Difference): [", round(y_lim_matrix[1], 4), ", ", round(y_lim_matrix[2], 4), "]\n", sep="")
cat("--------------------------------------------------\n\n")


# --- 3. Generate N x N MA Plot Matrix (Single PDF for all samples, all chromosomes combined) ---
# This creates a grid of plots where each row represents a 'reference' sample and
# each column represents a 'comparison' sample.
library(grDevices) # Ensure pdf() is available

# Define plot layout: N rows by N columns
num_samples_N <- length(sample_names)
plot_rows <- num_samples_N
plot_cols <- num_samples_N

pdf_file_name <- file.path(out_dir, "normalized_ev.100k_MA_plot_matrix_all_chroms_06_04_2025.pdf")
pdf(file = pdf_file_name,
    width = plot_cols * 5, # e.g., 6 columns * 5 inches/column = 30 inches wide
    height = plot_rows * 5) # e.g., 6 rows * 5 inches/row = 30 inches tall
par(mfrow = c(plot_rows, plot_cols), # Arrange plots in N rows by N columns (e.g., 6x6 grid)
    font.lab = 2, cex.lab = 1.0, # Bold labels, size 1.0
    mar = c(3, 3, 2, 1), # Margins for each individual plot: bottom, left, top, right (in lines)
    oma = c(2, 2, 2, 2), # Outer margins for the entire figure (for common titles/labels)
    mgp = c(1.5, 0.5, 0), # Axis title and label positions relative to axis line
    xaxs = 'i', yaxs = 'i' # Ensures plotting region extends exactly to axis limits
)

# Outer loop for rows (Sample on the Y-axis/Reference for difference calculation)
for (s_row_name in sample_names) {
  # Inner loop for columns (Sample on the X-axis/Comparison for difference calculation)
  for (s_col_name in sample_names) {
    
    # Get concatenated EV vectors for the current pair of samples
    ev_row <- concatenated_normalized_evs[[s_row_name]] # EV data for the sample corresponding to the row
    ev_col <- concatenated_normalized_evs[[s_col_name]] # EV data for the sample corresponding to the column
    
    # Handle NAs for common genomic bins
    common_idx_matrix <- !is.na(ev_row) & !is.na(ev_col)
    ev_row_valid <- ev_row[common_idx_matrix]
    ev_col_valid <- ev_col[common_idx_matrix]
    
    if (s_row_name == s_col_name) {
      # Plot diagonal elements: Simply display the sample name, as self-comparison is trivial
      plot(NA, xlim = x_lim_matrix, ylim = y_lim_matrix, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "", cex.main = 1.5, font.main = 2) # Set main title to blank, will add text manually
      text(x = mean(x_lim_matrix), y = mean(y_lim_matrix), labels = s_row_name, cex = 3, col = "gray40", font = 2)
      
    } else if (length(ev_row_valid) > 0 && length(ev_col_valid) > 0) {
      # Calculate M (Difference) and A (Average)
      # M = (Column Sample EV) - (Row Sample EV)
      # This convention makes M > 0 when the column sample has a higher EV than the row sample.
      M_val_matrix <- ev_col_valid - ev_row_valid
      A_val_matrix <- (ev_row_valid + ev_col_valid) / 2
      
      plot(A_val_matrix, M_val_matrix,
           pch = 16, cex = 0.25, # Point type and size for dense scatter plots
           col = rgb(0, 0, 0, alpha = 0.1), # Use semi-transparent black for point color to show density
           main = paste0(s_col_name, " vs ", s_row_name), # Title for individual plot
           xlab = "Mean EV", ylab = "EV Difference", # Axis labels for individual plot
           xlim = x_lim_matrix, # Use fixed global x-limits
           ylim = y_lim_matrix, # Use fixed global y-limits
           cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8) # Adjust font sizes within plot
      abline(h = 0, lty = "dashed", col = "red") # Horizontal line at EV Difference = 0
    } else {
      # Plot a placeholder if no common valid data points exist for this pair
      plot(NA, xlim = x_lim_matrix, ylim = y_lim_matrix, type = "n", xlab = "", ylab = "",
           main = paste0(s_col_name, " vs ", s_row_name, " (No Data)"), cex.main = 1.0)
    }
  }
}

# Add overall X-axis label for the entire figure
mtext("Mean EV (A)", side = 1, line = 0.5, outer = TRUE, cex = 1.5, font = 2) 
# Add overall Y-axis label for the entire figure
mtext("EV Difference (M)", side = 2, line = 0.5, outer = TRUE, cex = 1.5, font = 2) 
# Add an overall title for the entire figure
mtext("MA Plot Matrix of Normalized EVs Across Samples (All Chromosomes Combined)", side = 3, line = 0.5, outer = TRUE, cex = 1.8, font = 2) 

dev.off() # Close the PDF device
cat("\nMA Plot matrix (all samples, all chromosomes combined) saved to:", pdf_file_name, "\n")
