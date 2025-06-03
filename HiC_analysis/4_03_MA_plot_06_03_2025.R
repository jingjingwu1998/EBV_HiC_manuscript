#-----------------------------------
# This script is used to calculate MA plot for each chromosome from HiC data at 100k resolution across 6 samples.
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges) # Added for runmean

# output directory
out_dir <- 'normalized_ev/MA_plot' # Changed output directory to differentiate
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(out_dir, '../normalized_ev.100k_multi_samples_06_03_2025.rda'))

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


# MA plot: pairwise comparison across all samples
# This will now generate a separate PDF for each chromosome.
# Each PDF will contain (N * (N-1) / 2) plots for N samples.
# For 6 samples, this is 15 plots per chromosome PDF.

num_samples <- length(sample_names)
num_pairwise_plots <- num_samples * (num_samples - 1) / 2 

plots_per_row <- 3 # Desired number of columns in the plot grid
rows_per_page <- ceiling(num_pairwise_plots / plots_per_row) # Calculate required rows based on plot count



for(chr_name in chroms){ # Loop by chromosome name for consistent access
    # Adjust PDF dimensions to make individual plots roughly square in a 5x3 grid
    # A width of 15 inches for 3 columns gives 5 inches per column.
    # A height of 25 inches for 5 rows gives 5 inches per row.
    # This results in roughly 5x5 inch plots, which are square.
    pdf(file = file.path(out_dir, paste0(chr_name, '.ev.ma_pairwise_plots.pdf')),
        width = 15, # 3 columns * ~5 inches/column
        height = 25) # 5 rows * ~5 inches/row

    # Set mfrow to arrange 5 rows and 3 columns of plots on the page
    par(mfrow=c(rows_per_page, plots_per_row),
        font.lab=2, cex.lab=1.2,
        mar=c(4, 4, 3, 1), # Added margins (bottom, left, top, right)
        oma=c(0, 0, 0, 0), # Outer margins
        mgp = c(2.5, 0.8, 0), # Adjusted axis title and label positions
        xaxs='i', yaxs='i' # Ensures plotting region extends to axis limits
    ) # <--- ADDED THIS CLOSING PARENTHESIS FOR par()

    # The rest of your code block follows below, unchanged
    for (s1_idx in 1:(length(sample_names) - 1)) {
        for (s2_idx in (s1_idx + 1):length(sample_names)) {
            sample1_name <- sample_names[s1_idx]
            sample2_name <- sample_names[s2_idx]

            ev1 <- all_ev_lists[[sample1_name]][[chr_name]]
            ev2 <- all_ev_lists[[sample2_name]][[chr_name]]

            # Ensure both EV vectors are of the same length and contain valid numbers
            common_idx <- which(!is.na(ev1) & !is.na(ev2))
            
            if (length(common_idx) > 0) {
                EV_diff <- ev2[common_idx] - ev1[common_idx]
                EV_mean <- (ev1[common_idx] + ev2[common_idx]) / 2
                
                # Determine appropriate x and y limits for the plot
                x_limit <- c(-0.1, 0.1) # Mean EV ranges from -0.1 to 0.1
                y_limit <- c(-0.2, 0.2) # Diff EV ranges from -0.2 to 0.2

                plot(EV_mean, EV_diff, 
                     pch=16, cex=0.25, # Smaller points for dense plots
                     main=paste0(chr_name, ": ", sample2_name, " vs ", sample1_name),
                     xlab = "Mean EV", ylab = "EV Difference",
                     xlim = x_limit, # Set consistent x-axis limits
                     ylim = y_limit, # Set consistent y-axis limits
                     cex.main=1.0, cex.lab=1.0, cex.axis=0.8) # Adjust font sizes within plot
                abline(h=0,lty="dashed", col="red") # Make horizontal line red for visibility
            } else {
                # Plot placeholder if no common data points exist
                plot(NA, xlim=c(min(x_limit), max(x_limit)), ylim=c(min(y_limit), max(y_limit)), # Use consistent limits for empty plots too
                     type="n", xlab="", ylab="",
                     main=paste0(chr_name, ": ", sample2_name, " vs ", sample1_name, " (No common data)"),
                     cex.main=1.0)
            }
        }
    }
    dev.off() # Close the PDF for the current chromosome
}
