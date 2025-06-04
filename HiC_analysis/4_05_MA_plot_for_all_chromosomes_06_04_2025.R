#-----------------------------------
# Robust concatenation that forces chromosome order (chr1, chr2, ...
# This script is used to generate a comprehensive MA plot matrix
# for normalized HiC Eigenvector (EV) data at 100k resolution across 6 samples.
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges) # For runmean functionality (if used in prior steps)
library(grDevices) # Ensure pdf() is available for plotting

# --- 1. Setup and Data Loading ---

# Define output directory for plots and saved data
# IMPORTANT: This should be a relative path from your current working directory (getwd())
# e.g., if your script is in /path/to/project/ and you want output in /path/to/project/figures/MA_plot/, set 'figures/MA_plot'
out_dir <- 'normalized_ev/MA_plot_for_all_chromosomes' 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # Create directory if it doesn't exist

# Load normalized EV data
# This file contains individual objects like EV.rbl, EV.lcl, etc., with normalized values.
# Adjust the path relative to your script's location.
# Example: If your .rda file is directly in the parent directory of 'out_dir'
load(file.path(out_dir, '../normalized_ev.100k_multi_samples_06_03_2025.rda'))

# List all NORMALIZED EV objects for easy iteration
# These are the *normalized* EV values loaded from the .rda file.
normalized_ev_lists <- list(
    RBL = EV.rbl,
    LCL = EV.lcl,
    GCBC = EV.gcbc,
    MBC = EV.mbc,
    NBC = EV.nbc,
    PC = EV.pc
)

# Define sample names (derived directly from the loaded data structure)
sample_names <- names(normalized_ev_lists) 

# Define standard human chromosome names in correct order.
# This vector dictates the canonical order for concatenation.
chroms <- paste0('chr',c(1:22,'X'))


# --- 2. Concatenate All Chromosomes per Sample into a Single Long Vector (Robust Fix) ---
# This block creates a flattened representation of each sample's EV data,
# ensuring correct chromosome order and consistent length across all samples.

# Step 2A: Pre-compute expected lengths for each chromosome's EV vector.
# This is crucial for correctly padding with NAs if a chromosome is missing for a sample.
# We use 'RBL' as a template; ensure it's a sample expected to have all chromosomes with data.
chromosome_ev_lengths <- list()
for (chr_name in chroms) {
    if (!is.null(normalized_ev_lists$RBL[[chr_name]]) && length(normalized_ev_lists$RBL[[chr_name]]) > 0) {
        chromosome_ev_lengths[[chr_name]] <- length(normalized_ev_lists$RBL[[chr_name]])
    } else {
        # IMPORTANT: If a chromosome is genuinely missing/empty even in the template,
        # its length needs to be set appropriately (e.g., 0, or a typical length if it's just temporarily missing).
        # Adjust '0' below if such chromosomes are expected to have a specific length.
        warning(paste("Chromosome", chr_name, "is NULL or empty in template sample RBL. Assuming 0 length for padding if missing from other samples."))
        chromosome_ev_lengths[[chr_name]] <- 0 # Default to 0 length for padding if truly absent.
    }
}
cat("Calculated expected chromosome EV lengths for robust alignment.\n")


# Step 2B: Perform concatenation, forcing chromosome order and handling missing data.
concatenated_normalized_evs <- list()

for (sample_name in names(normalized_ev_lists)) {
  # Create a list to hold chromosome vectors IN THE DESIRED ORDER (chr1, chr2, ...)
  # This list will be exactly 'length(chroms)' long.
  ordered_chromosome_vectors_for_concat <- vector("list", length(chroms)) 
  names(ordered_chromosome_vectors_for_concat) <- chroms # Assign names (chr1, chr2...) to the list elements
  
  for (k in seq_along(chroms)) { # Loop using index 'k' to explicitly guarantee order
    chr_name <- chroms[k]
    ev_vec <- normalized_ev_lists[[sample_name]][[chr_name]] # Access chromosome data by name from source list
    
    if (!is.null(ev_vec) && length(ev_vec) > 0) {
      # Assign the actual EV vector if it's present and not empty
      ordered_chromosome_vectors_for_concat[[k]] <- ev_vec 
    } else {
      # If the chromosome is missing (NULL or empty) for this sample,
      # replace it with a vector of NAs of the expected length.
      # This is CRITICAL for maintaining consistent length and alignment across all samples.
      expected_len <- chromosome_ev_lengths[[chr_name]]
      
      if (is.null(expected_len) || expected_len == 0) {
          # If length couldn't be determined or is 0 (for the chromosome),
          # insert an empty numeric vector. This means such chromosomes are truly ignored in length calculation.
          ordered_chromosome_vectors_for_concat[[k]] <- numeric(0) 
      } else {
          # Fill with NAs of the expected length for this chromosome
          ordered_chromosome_vectors_for_concat[[k]] <- rep(NA_real_, expected_len)
      }
    }
  }
  
  # Concatenate the explicitly ordered and padded list of chromosome vectors.
  # The 'do.call(c, ...)' will now correctly combine elements in 'chroms' order.
  concatenated_vector_for_sample <- do.call(c, ordered_chromosome_vectors_for_concat)
  
  # IMPORTANT: Remove the misleading 'names' attribute.
  # This attribute is typically created by do.call(c, ...) from the names of the list elements
  # (e.g., "chr1", "chr2"), which would generate names like "chr1", "chr1", "chr1", etc. for bins.
  # Setting names to NULL ensures a clean numeric vector for MA plotting.
  names(concatenated_vector_for_sample) <- NULL 
  
  concatenated_normalized_evs[[sample_name]] <- concatenated_vector_for_sample
}
cat("Concatenated normalized EV data with enforced chromosome order and NA padding.\n\n")

# Final check to confirm all concatenated vectors have the exact same length
sample_lengths_check <- sapply(concatenated_normalized_evs, length)
if (length(unique(sample_lengths_check)) == 1) {
  cat("SUCCESS: All concatenated sample EV vectors have the same length:", unique(sample_lengths_check), "\n")
} else {
  cat("ERROR: Concatenated sample EV vectors have different lengths! Alignment will be incorrect.\n")
  print(sample_lengths_check)
}


# --- 3. Save the Concatenated Data to a Single File (.rda) ---
# This saves the 'concatenated_normalized_evs' object (your unified data for MA plots) to disk.
save(concatenated_normalized_evs, file = file.path(out_dir, "concatenated_normalized_evs_100k_multi_samples.rda"))
cat("Concatenated normalized EV data saved to:", file.path(out_dir, "concatenated_normalized_evs_100k_multi_samples.rda"), "\n\n")


# --- 4. Generate N x N MA Plot Matrix (Single PDF for all samples, all chromosomes combined) ---
# This creates a grid of plots where each row represents a 'reference' sample and
# each column represents a 'comparison' sample.

# Define fixed limits for X and Y axes.
# These values are arbitrary and might need adjustment based on your data's actual range.
# If you want dynamic limits, re-add the 'Calculate optimal global axis limits' block.
x_lim_matrix <- c(-0.1, 0.1)  
y_lim_matrix <- c(-0.2, 0.2)  
cat("Fixed MA Plot Matrix X-axis limits (Mean EV): [", round(x_lim_matrix[1], 4), ", ", round(x_lim_matrix[2], 4), "]\n", sep="")
cat("Fixed MA Plot Matrix Y-axis limits (EV Difference): [", round(y_lim_matrix[1], 4), ", ", round(y_lim_matrix[2], 4), "]\n", sep="")
cat("--------------------------------------------------\n\n")

# Define plot layout for N rows by N columns
num_samples_N <- length(sample_names)
plot_rows <- num_samples_N
plot_cols <- num_samples_N

pdf_file_name <- file.path(out_dir, "normalized_ev.100k_MA_plot_matrix_all_chroms_06_04_2025.pdf")
pdf(file = pdf_file_name,
    width = plot_cols * 5, 
    height = plot_rows * 5) 

# Set graphical parameters for plotting multiple plots on one page
par(mfrow = c(plot_rows, plot_cols), 
    font.lab = 2, cex.lab = 1.0, 
    mar = c(3, 3, 2, 1), 
    oma = c(2, 2, 2, 2), 
    mgp = c(1.5, 0.5, 0), 
    xaxs = 'i', yaxs = 'i' 
)

# Outer loop for rows (Sample on the Y-axis/Reference for difference calculation)
for (s_row_name in sample_names) {
  # Inner loop for columns (Sample on the X-axis/Comparison for difference calculation)
  for (s_col_name in sample_names) {
    
    ev_row <- concatenated_normalized_evs[[s_row_name]] 
    ev_col <- concatenated_normalized_evs[[s_col_name]] 
    
    # Handle NAs for common genomic bins
    common_idx_matrix <- !is.na(ev_row) & !is.na(ev_col)
    ev_row_valid <- ev_row[common_idx_matrix]
    ev_col_valid <- ev_col[common_idx_matrix]
    
    # --- Plotting Logic for each cell in the matrix ---
    if (s_row_name == s_col_name) {
      # Diagonal plots: Self-comparison (display sample name)
      plot(NA, xlim = x_lim_matrix, ylim = y_lim_matrix, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "", cex.main = 1.5, font.main = 2) 
      text(x = mean(x_lim_matrix), y = mean(y_lim_matrix), labels = s_row_name, cex = 3, col = "gray40", font = 2)
      
    } else if (length(ev_row_valid) > 0 && length(ev_col_valid) > 0) {
      # Off-diagonal plots with data: Pairwise MA comparison
      M_val_matrix <- ev_col_valid - ev_row_valid # M = (Column Sample EV) - (Row Sample EV)
      A_val_matrix <- (ev_row_valid + ev_col_valid) / 2
      
      plot(A_val_matrix, M_val_matrix,
           pch = 16, cex = 0.25, 
           col = rgb(0, 0, 0, alpha = 0.1), 
           main = paste0(s_col_name, " vs ", s_row_name), 
           xlab = "Mean EV", ylab = "EV Difference", 
           xlim = x_lim_matrix, 
           ylim = y_lim_matrix, 
           cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8) 
      abline(h = 0, lty = "dashed", col = "red") 
      
    } else {
      # Off-diagonal plots with no data: Placeholder
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
