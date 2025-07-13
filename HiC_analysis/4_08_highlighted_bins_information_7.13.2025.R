
library(dplyr)
# --- Initial Setup ---
setwd('/Users/80030577/Desktop/HiC_analysis/EV_normalization')
getwd()
load("normalized_ev.100k_multi_samples_06_03_2025.rda")

normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

sample_names <- names(normalized_ev_lists)

chroms <- paste0('chr',c(1:22,'X'))

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

# Assign the concatenated data to my_eigen_data for subsequent analysis
my_eigen_data <- concatenated_normalized_evs
sample_lengths <- sapply(my_eigen_data, length)
print("Lengths of concatenated eigenvector data for each sample:")
print(sample_lengths)

# Verify data points of concatenated_normalized_evs
all_M_values_across_all_pairs_and_chrs <- c()
sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)

for (pair in sample_pairs_to_compare) {
  sample1_name <- pair[1]
  sample2_name <- pair[2]
  
  for (chr in chroms) {
    ev_values_sample1_chr <- normalized_ev_lists[[sample1_name]][[chr]]
    ev_values_sample2_chr <- normalized_ev_lists[[sample2_name]][[chr]]
    M_values_chr <- ev_values_sample2_chr - ev_values_sample1_chr # Sample2 - Sample1
    all_M_values_across_all_pairs_and_chrs <- c(all_M_values_across_all_pairs_and_chrs, M_values_chr)
  }
}
str(all_M_values_across_all_pairs_and_chrs) # Should be 30,376 (total bins) * 15 (pairs) = 455,640
cat("Total M-values calculated across all pairs and chromosomes:", length(all_M_values_across_all_pairs_and_chrs), "\n")


####################################################################
### Core Data Preparation and Filtering
####################################################################

### Pre-calculate Chromosome Bin Mapping for Location Tracking
cat("--- Pre-calculating Chromosome Bin Mapping ---\n")
chr_bin_map <- data.frame(
  chromosome = character(),
  start_idx = numeric(),
  end_idx = numeric(),
  bin_in_chr_start = numeric(), # First bin number within that chromosome
  stringsAsFactors = FALSE
)
current_global_idx = 1
for (chr in chroms) {
  # Assuming all samples have the same binning structure, use the first sample to map
  chr_length = length(normalized_ev_lists[[sample_names[1]]][[chr]])
  if (chr_length > 0) {
    chr_bin_map <- rbind(chr_bin_map, data.frame(
      chromosome = chr,
      start_idx = current_global_idx,
      end_idx = current_global_idx + chr_length - 1,
      bin_in_chr_start = 1, # Bin numbering starts from 1 for each chr
      stringsAsFactors = FALSE
    ))
    current_global_idx = current_global_idx + chr_length
  }
}
cat("Chromosome bin mapping complete.\n")
print("Chromosome Bin Map:")
print(chr_bin_map)


# --- Compute the Global Cutoff Value for Highlighting ---
all_pairwise_sd_of_differences <- c()
# Note: sample_pairs_to_compare is already defined above, but re-defining for clarity
# sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)

for (pair in sample_pairs_to_compare) {
  sample1_name <- pair[1]
  sample2_name <- pair[2]
  
  ev_values_sample1 <- my_eigen_data[[sample1_name]]
  ev_values_sample2 <- my_eigen_data[[sample2_name]]
  
  M_values <- ev_values_sample1 - ev_values_sample2
  sd_of_current_M <- sd(M_values, na.rm = TRUE)
  all_pairwise_sd_of_differences <- c(all_pairwise_sd_of_differences, sd_of_current_M)
}

cutoff_value <- mean(all_pairwise_sd_of_differences, na.rm = TRUE)
cutoff_value <- 2 * cutoff_value # 2 times the mean SD as the cutoff
print(paste0("Calculated Mean SD of M-values: ", mean(all_pairwise_sd_of_differences, na.rm = TRUE)))
cat("-------------------------------------------------------------------\n")
cat("Calculated Cutoff Value for identifying heavily changed bins:", cutoff_value, "\n")
cat("-------------------------------------------------------------------\n\n")

num_samples_N <- length(sample_names) # 6 samples


# These matrices will store the R-squared values
r_squared_all_points_matrix <- matrix(NA, nrow = num_samples_N, ncol = num_samples_N,
                                      dimnames = list(sample_names, sample_names))
r_squared_highlighted_points_matrix <- matrix(NA, nrow = num_samples_N, ncol = num_samples_N,
                                              dimnames = list(sample_names, sample_names))

# NEW: Matrix to store counts of highlighted bins
highlighted_bin_counts_matrix <- matrix(NA, nrow = num_samples_N, ncol = num_samples_N,
                                        dimnames = list(sample_names, sample_names))


# Minimum number of points required for meaningful R-squared calculation
min_points_for_lm <- 2 # At least 2 points are needed for lm(), 3 or more for robust fitting

cat("\n--- Step 1: Calculating R-squared for ALL points ---\n")
# Loop through all sample pairs to calculate R-squared for all points
for (i in 1:num_samples_N) {
  s_row_name <- sample_names[i]
  data_row <- my_eigen_data[[s_row_name]]
  
  for (j in 1:num_samples_N) {
    s_col_name <- sample_names[j]
    data_col <- my_eigen_data[[s_col_name]]
    
    # Handle diagonal cases (R-squared is 1.0 for self-comparison)
    if (s_row_name == s_col_name) {
      r_squared_all_points_matrix[s_row_name, s_col_name] <- 1.0
    } else {
      # Identify valid (non-NA, non-Inf) indices for both samples
      valid_indices_full <- !is.na(data_row) & !is.na(data_col) & is.finite(data_row) & is.finite(data_col)
      
      if (sum(valid_indices_full) >= min_points_for_lm) {
        lm_all_points <- lm(data_col[valid_indices_full] ~ data_row[valid_indices_full])
        r_squared_all_points_matrix[s_row_name, s_col_name] <- summary(lm_all_points)$r.squared
      } else {
        r_squared_all_points_matrix[s_row_name, s_col_name] <- NA # Not enough valid points
      }
    }
  }
}
cat("R-squared for all points calculation complete.\n")


cat("\n--- Step 2: Identifying and Locating Highlighted Points, and Calculating their R-squared ---\n")

# Initialize a list to store highlighted points' locations across all pairs
all_highlighted_points_locations <- list()

# Loop through all sample pairs to identify highlighted points, calculate their R-squared, and locate them
for (i in 1:num_samples_N) {
  s_row_name <- sample_names[i]
  data_row <- my_eigen_data[[s_row_name]]
  
  for (j in 1:num_samples_N) {
    s_col_name <- sample_names[j]
    data_col <- my_eigen_data[[s_col_name]]
    
    # Handle diagonal cases (R-squared is 1.0 for self-comparison)
    if (s_row_name == s_col_name) {
      r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- 1.0
      highlighted_bin_counts_matrix[s_row_name, s_col_name] <- 0 # No highlighted bins for self-comparison
    } else {
      # Identify valid (non-NA, non-Inf) indices for both samples
      valid_indices_full <- !is.na(data_row) & !is.na(data_col) & is.finite(data_row) & is.finite(data_col)
      
      # Calculate M-values for the valid subset to identify highlighted points
      M_val_for_highlight_check <- data_col[valid_indices_full] - data_row[valid_indices_full]
      
      # Identify points beyond the cutoff from the *valid* subset
      points_to_highlight_indices_in_valid_subset <- which(M_val_for_highlight_check > cutoff_value | M_val_for_highlight_check < -cutoff_value)
      
      # NEW: Count the number of highlighted bins for this pairwise comparison
      num_highlighted_bins <- length(points_to_highlight_indices_in_valid_subset)
      highlighted_bin_counts_matrix[s_row_name, s_col_name] <- num_highlighted_bins
      
      # --- Calculate R-squared for HIGHLIGHTED points ---
      highlighted_data_row_values <- data_row[valid_indices_full][points_to_highlight_indices_in_valid_subset]
      highlighted_data_col_values <- data_col[valid_indices_full][points_to_highlight_indices_in_valid_subset]
      
      if (length(highlighted_data_row_values) >= min_points_for_lm) {
        # Ensure there are enough finite points after subsetting
        valid_lm_subset_indices <- is.finite(highlighted_data_row_values) & is.finite(highlighted_data_col_values)
        
        if (sum(valid_lm_subset_indices) >= min_points_for_lm) {
          # Linear regression between the eigenvector values themselves for highlighted points
          lm_highlighted <- lm(highlighted_data_col_values[valid_lm_subset_indices] ~ highlighted_data_row_values[valid_lm_subset_indices])
          r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- summary(lm_highlighted)$r.squared
        } else {
          r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- NA # Not enough valid points
        }
      } else {
        r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- NA # Assign NA if no highlighted points
      }
      
      # --- Locate highlighted points ---
      if (length(points_to_highlight_indices_in_valid_subset) > 0) {
        # Get the actual global indices of the valid points from the original concatenated vector
        global_valid_indices <- which(valid_indices_full)
        
        # Now, get the original global indices of the highlighted points
        original_global_highlighted_indices <- global_valid_indices[points_to_highlight_indices_in_valid_subset]
        
        # Get the actual normalized values for the highlighted points from the original data
        # These are the values from the full concatenated vectors, indexed by their original global positions
        sample1_highlighted_values <- data_row[original_global_highlighted_indices]
        sample2_highlighted_values <- data_col[original_global_highlighted_indices]
        
        # Calculate average (A-value) for highlighted points
        average_highlighted_values <- (sample1_highlighted_values + sample2_highlighted_values) / 2
        
        # Create a temporary data frame for this pair's highlighted points
        temp_highlighted_df <- data.frame(
          sample1 = s_row_name,
          sample2 = s_col_name,
          global_index = original_global_highlighted_indices,
          M_value = M_val_for_highlight_check[points_to_highlight_indices_in_valid_subset], # M-value is already calculated
          sample1_value = sample1_highlighted_values, # NEW: Normalized value for sample1
          sample2_value = sample2_highlighted_values, # NEW: Normalized value for sample2
          average_value = average_highlighted_values, # NEW: Average value (A-value)
          chromosome = character(length(original_global_highlighted_indices)),
          bin_in_chromosome = numeric(length(original_global_highlighted_indices)),
          stringsAsFactors = FALSE
        )
        
        # Determine chromosome and bin for each highlighted global index
        for (k in 1:length(original_global_highlighted_indices)) {
          current_global_idx <- original_global_highlighted_indices[k]
          
          # Find which chromosome this global index belongs to using chr_bin_map
          # We search for the row in chr_bin_map where current_global_idx falls within [start_idx, end_idx]
          chr_info_row <- chr_bin_map[current_global_idx >= chr_bin_map$start_idx & current_global_idx <= chr_bin_map$end_idx, ]
          
          if (nrow(chr_info_row) == 1) {
            temp_highlighted_df$chromosome[k] <- chr_info_row$chromosome
            # Calculate bin number within the chromosome: global_index - start_of_chr_global_index + 1
            temp_highlighted_df$bin_in_chromosome[k] <- current_global_idx - chr_info_row$start_idx + chr_info_row$bin_in_chr_start
          } else {
            # This case indicates an issue with mapping or data if chr_info_row is not found or multiple rows match
            temp_highlighted_df$chromosome[k] <- NA
            temp_highlighted_df$bin_in_chromosome[k] <- NA
            warning(paste("Could not map global index", current_global_idx, "to a unique chromosome bin."))
          }
        }
        # Add to the list of all highlighted points locations
        all_highlighted_points_locations[[paste0(s_row_name, "_vs_", s_col_name)]] <- temp_highlighted_df
      }
    }
  }
}
cat("Highlighted point identification and R-squared calculation complete.\n")

# Combine all highlighted points into a single data frame
final_highlighted_points_df <- do.call(rbind, all_highlighted_points_locations)

# Print the results
cat("\n--- Locations of Highlighted Points (M-value beyond 2 * Mean SD) ---\n")
if (!is.null(final_highlighted_points_df) && nrow(final_highlighted_points_df) > 0) {
  print(final_highlighted_points_df)
} else {
  cat("No highlighted points found based on the current cutoff across all sample pairs.\n")
}

# You can now print these matrices or use them for heatmaps as discussed previously
print("\nR-squared matrix for all points:")
print(r_squared_all_points_matrix)

print("\nR-squared matrix for highlighted points:")
print(r_squared_highlighted_points_matrix)

# NEW: Print the matrix of highlighted bin counts
print("\nNumber of Highlighted Bins per Pairwise Comparison:")
print(highlighted_bin_counts_matrix)


# --- Sanity Check: Verify M_value against cutoff ---
# This script assumes 'final_highlighted_points_df' and 'cutoff_value'
# have already been generated by the main R script.

cat("\n--- Sanity Check: M_value vs. Cutoff Value ---\n")

# Check if final_highlighted_points_df exists and has rows
if (!is.null(final_highlighted_points_df) && nrow(final_highlighted_points_df) > 0) {
  # Check if all M_values are indeed outside the cutoff range
  # abs() is used because points can be highlighted if M_value > cutoff or M_value < -cutoff
  m_values_check <- abs(final_highlighted_points_df$M_value) > cutoff_value
  
  if (all(m_values_check)) {
    cat("Sanity Check Passed: All M_values in 'final_highlighted_points_df' are correctly beyond the 2*SD cutoff.\n")
  } else {
    cat("WARNING: Sanity Check FAILED. Some M_values in 'final_highlighted_points_df' are NOT beyond the 2*SD cutoff.\n")
    # Optionally, print the rows that failed the check for detailed inspection
    print(final_highlighted_points_df[!m_values_check, c("sample1", "sample2", "global_index", "M_value", "chromosome", "bin_in_chromosome")])
  }
} else {
  cat("Sanity Check Skipped: No highlighted points to check (final_highlighted_points_df is empty or NULL).\n")
}
cat("-------------------------------------------------------------------\n")
