setwd("/Users/80030577/Desktop/HiC_analysis/EV_normalization")
load("normalized_ev.100k_multi_samples_06_03_2025.rda")

library(GenomicRanges)
library(IRanges)
library(S4Vectors)

normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

sample_names <- names(normalized_ev_lists)
chroms <- paste0("chr", c(1:22, "X"))
bin_size <- 100000  # 100kb

########################################################################
######## Step 1: Create GRanges of Eigenvectors per Sample ############
########################################################################

make_ev_gr <- function(ev_by_chr, sample_name, bin_size = 100000) {
  gr_list <- lapply(names(ev_by_chr), function(chr) {
    ev_vec <- ev_by_chr[[chr]]
    starts <- seq(0, (length(ev_vec) - 1) * bin_size, by = bin_size)
    ends <- starts + bin_size
    
    GRanges(
      seqnames = chr,
      ranges = IRanges(start = starts + 1, end = ends),
      eigen = ev_vec,
      sample = sample_name
    )
  })
  names(gr_list) <- names(ev_by_chr)
  return(unlist(GRangesList(gr_list)))
}

# Create GRanges for all samples
gr_by_sample <- lapply(sample_names, function(sname) {
  make_ev_gr(normalized_ev_lists[[sname]], sname)
})
names(gr_by_sample) <- sample_names

str(gr_by_sample, max.level = 1)

########################################################################
####### Step 2: Compute ΔEV for Each Directed Sample Pair #############
########################################################################

pair_deltas <- list()

for (s1 in sample_names) {
  for (s2 in sample_names) {
    if (s1 != s2) {
      gr1 <- gr_by_sample[[s1]]
      gr2 <- gr_by_sample[[s2]]
      
      # Confirm bin alignment
      stopifnot(identical(seqnames(gr1), seqnames(gr2)))
      stopifnot(identical(start(gr1), start(gr2)))
      
      delta <- gr2$eigen - gr1$eigen
      gr_delta <- gr1  # structure/template from sample 1
      gr_delta$deltaEV <- delta
      gr_delta$pair <- paste0(s1, "_to_", s2)
      pair_deltas[[paste0(s1, "_to_", s2)]] <- gr_delta
    }
  }
}
str(pair_deltas, max.level = 1)

########################################################################
### Step 3: Compare Transitions (RBL→GCBC vs NBC→GCBC) ################
########################################################################

delta_rbl_gcbc <- pair_deltas[["RBL_to_GCBC"]]
delta_nbc_gcbc <- pair_deltas[["NBC_to_GCBC"]]

# Assign direction signs to each bin
delta_rbl_gcbc$sign <- sign(delta_rbl_gcbc$deltaEV)
delta_nbc_gcbc$sign <- sign(delta_nbc_gcbc$deltaEV)

# Identify bins with different sign direction
same_direction <- delta_rbl_gcbc$sign == delta_nbc_gcbc$sign
different <- delta_rbl_gcbc[!same_direction]

########################################################################
### Step 4: Group Consecutive Bins with Same Sign #####################
########################################################################

group_consecutive_sign <- function(gr, sign_col = "sign") {
  # e.g. for: +1, +1, +1, -1, -1, +1, +1, you get:
  # runLength = [3, 2, 2]
  # runValue = [+1, -1, +1]
  signs <- mcols(gr)[[sign_col]]
  rle_signs <- Rle(signs)

  
  # Safe group assignment using inverse.rle
  # This uses inverse.rle() to expand a sequence like:
  # values = [1, 2, 3]
  # lengths = [3, 2, 2]
  # Into a vector like: 1, 1, 1, 2, 2, 3, 3
  # This gives each consecutive run of identical signs a unique group ID.
  group_ids <- inverse.rle(list(
    values = seq_along(runLength(rle_signs)),
    lengths = runLength(rle_signs)
  ))
  
  # Assign group IDs to bins
  gr$group <- group_ids
  
  # Merge bins by group using reduce()
  # split(gr, gr$group) splits the GRanges object into a list, one element per group.
  # reduce() merges consecutive genomic ranges in each group into a single interval.
  # unlist() flattens the list back into one GRanges object.
  merged <- reduce(split(gr, gr$group))
  merged <- unlist(merged)
  
  # Assign sign back to merged bins
  # Restore the sign label to the merged bins
  merged$sign <- signs[match(names(merged), names(gr))]
  
  return(merged)
}

# Apply grouping function
combined_rbl_gcbc <- group_consecutive_sign(delta_rbl_gcbc)
combined_nbc_gcbc <- group_consecutive_sign(delta_nbc_gcbc)

########################################################################
# Step 5: Output Bins with Directional Changes Between RBL→GCBC and NBC→GCBC
########################################################################

# Identify bins with same or different ΔEV sign
same_sign_rbl_gcbc_vs_nbc_gcbc <- delta_rbl_gcbc$sign == delta_nbc_gcbc$sign

# ------------------- DIFFERENTLY ALTERED BINS -----------------------

# Bins where sign differs
diff_bins_rbl_gcbc_vs_nbc_gcbc <- delta_rbl_gcbc[!same_sign_rbl_gcbc_vs_nbc_gcbc]

# Merge consecutive bins with same sign
merged_diff_bins_rbl_gcbc_vs_nbc_gcbc <- group_consecutive_sign(diff_bins_rbl_gcbc_vs_nbc_gcbc)

# Convert to data frame
df_diff_rbl_gcbc_vs_nbc_gcbc <- as.data.frame(merged_diff_bins_rbl_gcbc_vs_nbc_gcbc, row.names = NULL)
head(df_diff_rbl_gcbc_vs_nbc_gcbc[, c("seqnames", "start", "end", "sign")])

# ------------------- SIMILARLY ALTERED BINS -------------------------

# Bins where sign is the same
same_bins_rbl_gcbc_vs_nbc_gcbc <- delta_rbl_gcbc[same_sign_rbl_gcbc_vs_nbc_gcbc]

# Merge consecutive bins with same sign
merged_same_bins_rbl_gcbc_vs_nbc_gcbc <- group_consecutive_sign(same_bins_rbl_gcbc_vs_nbc_gcbc)

# Convert to data frame
df_same_rbl_gcbc_vs_nbc_gcbc <- as.data.frame(merged_same_bins_rbl_gcbc_vs_nbc_gcbc, row.names = NULL)
head(df_same_rbl_gcbc_vs_nbc_gcbc[, c("seqnames", "start", "end", "sign")])



########################################################################
# Step 6: Output Bins with Directional Changes Between RBL→LCL and NBC→LCL
########################################################################

# Extract GRanges for ΔEV
delta_rbl_lcl <- pair_deltas[["RBL_to_LCL"]]
delta_nbc_lcl <- pair_deltas[["NBC_to_LCL"]]

# Compute the sign of ΔEV
delta_rbl_lcl$sign <- sign(delta_rbl_lcl$deltaEV)
delta_nbc_lcl$sign <- sign(delta_nbc_lcl$deltaEV)

# Identify bins with same or different ΔEV sign
same_sign_rbl_lcl_vs_nbc_lcl <- delta_rbl_lcl$sign == delta_nbc_lcl$sign

# ------------------- DIFFERENTLY ALTERED BINS -----------------------

# Bins where sign differs
diff_bins_rbl_lcl_vs_nbc_lcl <- delta_rbl_lcl[!same_sign_rbl_lcl_vs_nbc_lcl]

# Merge consecutive bins with same sign
merged_diff_bins_rbl_lcl_vs_nbc_lcl <- group_consecutive_sign(diff_bins_rbl_lcl_vs_nbc_lcl)

# Convert to data frame
df_diff_rbl_lcl_vs_nbc_lcl <- as.data.frame(merged_diff_bins_rbl_lcl_vs_nbc_lcl, row.names = NULL)
head(df_diff_rbl_lcl_vs_nbc_lcl[, c("seqnames", "start", "end", "sign")])

# ------------------- SIMILARLY ALTERED BINS -------------------------

# Bins where sign is the same
same_bins_rbl_lcl_vs_nbc_lcl <- delta_rbl_lcl[same_sign_rbl_lcl_vs_nbc_lcl]

# Merge consecutive bins with same sign
merged_same_bins_rbl_lcl_vs_nbc_lcl <- group_consecutive_sign(same_bins_rbl_lcl_vs_nbc_lcl)

# Convert to data frame
df_same_rbl_lcl_vs_nbc_lcl <- as.data.frame(merged_same_bins_rbl_lcl_vs_nbc_lcl, row.names = NULL)
head(df_same_rbl_lcl_vs_nbc_lcl[, c("seqnames", "start", "end", "sign")])
