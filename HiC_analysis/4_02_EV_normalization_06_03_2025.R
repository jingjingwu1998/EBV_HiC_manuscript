library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges) # Added for runmean

# output directory
out_dir <- 'figures' # Changed output directory to differentiate
# Load saved EVs
load(file.path(out_dir, 'ev.100k_multi_samples_raw.rda'))

# List of all EV objects
all_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

sample_names <- names(all_ev_lists)
chroms <- paste0("chr", c(1:22, "X"))

# Initialize output structure
median_results <- list()

# Loop over samples and chromosomes
for (sample in sample_names) {
  pos_medians <- numeric(length(chroms))
  neg_medians <- numeric(length(chroms))
  
  for (i in seq_along(chroms)) {
    chr <- chroms[i]
    ev <- all_ev_lists[[sample]][[chr]]
    
    if (!is.null(ev)) {
      pos_vals <- ev[ev > 0]
      neg_vals <- ev[ev < 0]
      
      pos_medians[i] <- if (length(pos_vals) > 0) median(pos_vals, na.rm = TRUE) else NA
      neg_medians[i] <- if (length(neg_vals) > 0) median(neg_vals, na.rm = TRUE) else NA
    } else {
      pos_medians[i] <- NA
      neg_medians[i] <- NA
    }
  }
  
  median_results[[sample]] <- list(
    chroms = chroms,
    positive_medians = pos_medians,
    negative_medians = neg_medians
  )
}

# Example: view median results for RBL
median_results$RBL


# Step 1: Compute Global Mean of chromosome-level Median for Each Sample

global_means <- list()

for (sample in names(median_results)) {
  pos_vals <- median_results[[sample]]$positive_medians
  neg_vals <- median_results[[sample]]$negative_medians
  
  
  # Compute mean
  pos_mean <- if (length(pos_vals) > 0) mean(pos_vals) else NA
  neg_mean <- if (length(neg_vals) > 0) mean(neg_vals) else NA
  
  # Save to list
  global_means[[sample]] <- list(
    positive_global_mean = pos_mean,
    negative_global_mean = neg_mean
  )
}


####### step 2 need to double check and then need to get normalized EV and plot MA plot and histrogram
# Step 2: Normalize EV Values as per your specified strategy
# X(A) = mean of the median value = ( A + A’ + A’’ + A’’’ + …)/n
# Eg. scare down chr1 in sample 1 by A/X(A)
# a1 ÷ A/X(A), a2 ÷ A/X(A), etc.

normalized_ev_lists <- list()

# Loop through each sample
for (sample_name in names(all_ev_lists)) {
  # Get original EV data for the current sample
  original_sample_evs <- all_ev_lists[[sample_name]]
  
  # Get the pre-calculated per-chromosome medians for the current sample
  # These are A, A', A'', ... and B, B', B'', ... for the current sample
  current_sample_pos_medians <- median_results[[sample_name]]$positive_medians
  current_sample_neg_medians <- median_results[[sample_name]]$negative_medians
  
  # Get the sample-specific mean of medians (X(A) and X(B) for the current sample)
  sample_X_A <- global_means[[sample_name]]$positive_global_mean
  sample_X_B <- global_means[[sample_name]]$negative_global_mean
  
  # Initialize a list to store normalized chromosomes for the current sample
  norm_chr_list <- list()
  
  # Loop through chromosomes within the current sample
  for (chr_name in names(original_sample_evs)) {
    ev_vector <- original_sample_evs[[chr_name]] # Get the EV vector for this chromosome
    
    # Check if EV data exists for this chromosome
    if (!is.null(ev_vector)) {
      # Initialize the normalized vector with original values first
      # This ensures NA/other values not covered by positive/negative are preserved
      normalized_ev_vector <- ev_vector 
      
      # Find the index of the current chromosome in the 'chroms' list
      # This is crucial for correctly linking ev_vector to its corresponding median_results entry
      chr_idx <- which(chroms == chr_name) 
      
      # Retrieve the specific median for this chromosome and sample (A or B)
      median_A_for_this_chr_sample <- current_sample_pos_medians[chr_idx]
      median_B_for_this_chr_sample <- current_sample_neg_medians[chr_idx]
      
      # --- Handle Positive EV values ---
      # Check if scaling for positive values is possible (median is valid and non-zero)
      if (!is.na(median_A_for_this_chr_sample) && median_A_for_this_chr_sample != 0 && !is.na(sample_X_A)) {
        # Calculate the scaling factor for positive EVs
        scale_factor_pos <- sample_X_A / median_A_for_this_chr_sample
        
        # Identify positive values in the original EV vector
        pos_indices <- which(ev_vector > 0)
        
        # Apply scaling only to positive values using vectorized operation
        if (length(pos_indices) > 0) {
          normalized_ev_vector[pos_indices] <- ev_vector[pos_indices] * scale_factor_pos
        }
      } else {
        # If scaling is not possible for positive values, set them to NA
        pos_indices <- which(ev_vector > 0)
        if (length(pos_indices) > 0) {
          normalized_ev_vector[pos_indices] <- NA
        }
      }
      
      # --- Handle Negative EV values ---
      # Check if scaling for negative values is possible
      if (!is.na(median_B_for_this_chr_sample) && median_B_for_this_chr_sample != 0 && !is.na(sample_X_B)) {
        # Calculate the scaling factor for negative EVs
        scale_factor_neg <- sample_X_B / median_B_for_this_chr_sample
        
        # Identify negative values in the original EV vector
        neg_indices <- which(ev_vector < 0)
        
        # Apply scaling only to negative values using vectorized operation
        if (length(neg_indices) > 0) {
          normalized_ev_vector[neg_indices] <- ev_vector[neg_indices] * scale_factor_neg
        }
      } else {
        # If scaling is not possible for negative values, set them to NA
        neg_indices <- which(ev_vector < 0)
        if (length(neg_indices) > 0) {
          normalized_ev_vector[neg_indices] <- NA
        }
      }
      
      # --- Handle Zero EV values ---
      # Zeros should remain zero and are not affected by the above positive/negative filters.
      # They are already preserved from the initial `normalized_ev_vector <- ev_vector` line.
      
      # --- Any original NA values in ev_vector are also preserved as they are not affected by `which()`
      
      norm_chr_list[[chr_name]] <- normalized_ev_vector
    } else {
      # If original EV data for chromosome was NULL, keep it NULL
      norm_chr_list[[chr_name]] <- NULL 
    }
  } # End of inner loop (chromosomes)
  
  normalized_ev_lists[[sample_name]] <- norm_chr_list
} # End of outer loop (samples)


# --- ADDITIONAL LINES TO SAVE NORMALIZED EV DATA IN RAW INPUT FORMAT ---
# (This block is placed here, AFTER the normalization loop is complete,
# to ensure all samples are processed before saving the final output.)

# Create individual EV objects (e.g., EV.rbl, EV.lcl) in the global environment
# from the 'normalized_ev_lists' object.
for (sample_name_full in names(normalized_ev_lists)) {
  # Construct the desired object name (e.g., EV.rbl, EV.lcl by converting "RBL" to "rbl")
  object_name <- paste0("EV.", tolower(sample_name_full))
  
  # Assign the normalized EV list for the current sample to the newly created object name
  assign(object_name, normalized_ev_lists[[sample_name_full]], envir = .GlobalEnv)
}

# Define the list of object names to save to the .rda file.
# This ensures that only the desired 'EV.sample_name' objects are saved,
# similar to how your raw input data file is structured.
objects_to_save <- paste0("EV.", tolower(sample_names)) 

# Save these individual EV objects to the specified .rda file
save(list = objects_to_save, file = file.path(out_dir, "normalized_ev.100k_multi_samples_06_03_2025.rda"))
cat("Normalized EV lists saved to:", file.path(out_dir, "normalized_ev.100k_multi_samples_06_03_2025.rda"), "\n")
  
