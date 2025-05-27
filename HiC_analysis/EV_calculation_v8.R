#-----------------------------------
# This script is used to calculate EV value and identify compartment A and B of HiC at 100k resolution.
# And the generation of the compartment tracks in fig1a,fig1b and fig1c
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges) # Added for runmean

# output directory
out_dir <- 'figures/ab_fig_2_multi_samples_v8' # Changed output directory to differentiate
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# input
rbl <- importC('matrix/Day0/iced/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/Day0/raw/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/Day0/raw/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)
lcl <- importC('matrix/Day28/iced/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/Day28/raw/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/Day28/raw/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)
gcbc <- importC('matrix/GCBC/iced/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/GCBC/raw/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/GCBC/raw/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)
mbc <- importC('matrix/MBC/iced/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/MBC/raw/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/MBC/raw/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)
nbc <- importC('matrix/NBC/iced/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/NBC/raw/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/NBC/raw/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)
pc <- importC('matrix/PC/iced/100000/PC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',xgi.bed="matrix/PC/raw/100000/PC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",ygi.bed="matrix/PC/raw/100000/PC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",rm.trans=T)


binsize <- 100000

# manually set up EV sign
manual_sign_rbl <- list(
  chr1 = +1,
  chr2 = -1,
  chr3 = -1,
  chr4 = -1,
  chr5 = -1,
  chr6 = +1,
  chr7 = -1,
  chr8 = +1,
  chr9 = -1,
  chr10 = +1,
  chr11 = +1,
  chr12 = -1,
  chr13 = +1,
  chr14 = +1,
  chr15 = +1,
  chr16 = -1,
  chr17 = +1,
  chr18 = -1,
  chr19 = +1,
  chr20 = +1,
  chr21 = +1,
  chr22 = -1,
  chrX = +1
)

# calculate EV for each chromosome
EV.rbl <- EV.lcl <- EV.gcbc <- EV.mbc <- EV.nbc <- EV.pc <- list()
for(i in 23:1){ # Original loop for chromosomes (23:1)
    ## 1st approach (as per template)
    # Process RBL
    x <- normPerExpected(rbl[[i]], method="loess",stdev=TRUE)
    rbl.xdata <- as.matrix(intdata(forceSymmetric(x)))
    cat("RBL NA ratio chr", seqlevels(rbl[[i]]), ":", sum(is.na(rbl.xdata))/length(rbl.xdata), '\n')
    rbl.xdata[is.na(rbl.xdata)] = 0
    rbl.ev <- as.numeric(runmean(Rle(pca(rbl.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Manually flip RBL sign based on manual_sign_rbl
    chr <- seqlevels(rbl[[i]])
    sign_factor <- ifelse(!is.null(manual_sign_rbl[[chr]]), manual_sign_rbl[[chr]], 1)
    rbl.ev <- rbl.ev * sign_factor

    # Process LCL
    y <- normPerExpected(lcl[[i]], method="loess",stdev=TRUE)
    lcl.xdata <- as.matrix(intdata(forceSymmetric(y)))
    cat("LCL NA ratio chr", seqlevels(lcl[[i]]), ":", sum(is.na(lcl.xdata))/length(lcl.xdata), '\n')
    lcl.xdata[is.na(lcl.xdata)] = 0
    lcl.ev <- as.numeric(runmean(Rle(pca(lcl.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Process GCBC
    z <- normPerExpected(gcbc[[i]], method="loess",stdev=TRUE)
    gcbc.xdata <- as.matrix(intdata(forceSymmetric(z)))
    cat("GCBC NA ratio chr", seqlevels(gcbc[[i]]), ":", sum(is.na(gcbc.xdata))/length(gcbc.xdata), '\n')
    gcbc.xdata[is.na(gcbc.xdata)] = 0 # Corrected
    gcbc.ev <- as.numeric(runmean(Rle(pca(gcbc.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Process MBC
    m <- normPerExpected(mbc[[i]], method="loess",stdev=TRUE)
    mbc.xdata <- as.matrix(intdata(forceSymmetric(m)))
    cat("MBC NA ratio chr", seqlevels(mbc[[i]]), ":", sum(is.na(mbc.xdata))/length(mbc.xdata), '\n')
    mbc.xdata[is.na(mbc.xdata)] = 0
    mbc.ev <- as.numeric(runmean(Rle(pca(mbc.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Process NBC
    n <- normPerExpected(nbc[[i]], method="loess",stdev=TRUE)
    nbc.xdata <- as.matrix(intdata(forceSymmetric(n)))
    cat("NBC NA ratio chr", seqlevels(nbc[[i]]), ":", sum(is.na(nbc.xdata))/length(nbc.xdata), '\n')
    nbc.xdata[is.na(nbc.xdata)] = 0 # Corrected
    nbc.ev <- as.numeric(runmean(Rle(pca(nbc.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Process PC
    p <- normPerExpected(pc[[i]], method="loess",stdev=TRUE)
    pc.xdata <- as.matrix(intdata(forceSymmetric(p)))
    cat("PC NA ratio chr", seqlevels(pc[[i]]), ":", sum(is.na(pc.xdata))/length(pc.xdata), '\n')
    pc.xdata[is.na(pc.xdata)] = 0 # Corrected
    pc.ev <- as.numeric(runmean(Rle(pca(pc.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))

    # Flip signs of other samples based on RBL data
    # Use "complete.obs" for correlation to handle any remaining NAs gracefully.
    valid_idx_lcl <- which(!is.na(rbl.ev) & !is.na(lcl.ev))
    if(length(valid_idx_lcl) > 10 && cor(rbl.ev[valid_idx_lcl], lcl.ev[valid_idx_lcl], use="complete.obs") < 0) lcl.ev <- -lcl.ev

    valid_idx_gcbc <- which(!is.na(rbl.ev) & !is.na(gcbc.ev))
    if(length(valid_idx_gcbc) > 10 && cor(rbl.ev[valid_idx_gcbc], gcbc.ev[valid_idx_gcbc], use="complete.obs") < 0) gcbc.ev <- -gcbc.ev

    valid_idx_mbc <- which(!is.na(rbl.ev) & !is.na(mbc.ev))
    if(length(valid_idx_mbc) > 10 && cor(rbl.ev[valid_idx_mbc], mbc.ev[valid_idx_mbc], use="complete.obs") < 0) mbc.ev <- -mbc.ev

    valid_idx_nbc <- which(!is.na(rbl.ev) & !is.na(nbc.ev))
    if(length(valid_idx_nbc) > 10 && cor(rbl.ev[valid_idx_nbc], nbc.ev[valid_idx_nbc], use="complete.obs") < 0) nbc.ev <- -nbc.ev

    valid_idx_pc <- which(!is.na(rbl.ev) & !is.na(pc.ev))
    if(length(valid_idx_pc) > 10 && cor(rbl.ev[valid_idx_pc], pc.ev[valid_idx_pc], use="complete.obs") < 0) pc.ev <- -pc.ev

    EV.rbl[[seqlevels(rbl[[i]])]] <- rbl.ev
    EV.lcl[[seqlevels(lcl[[i]])]] <- lcl.ev
    EV.gcbc[[seqlevels(gcbc[[i]])]] <- gcbc.ev
    EV.mbc[[seqlevels(mbc[[i]])]] <- mbc.ev
    EV.nbc[[seqlevels(nbc[[i]])]] <- nbc.ev
    EV.pc[[seqlevels(pc[[i]])]] <- pc.ev

    # No Chip-seq correlation or plots in this version as per your instruction
    # The original template's correlation plots with covk4me3 are omitted.
}

# save all EV values
save(EV.rbl, EV.lcl, EV.gcbc, EV.mbc, EV.nbc, EV.pc, file = file.path(out_dir, 'ev.100k_multi_samples_v8.rda'))

# Original template loaded a different file name, keeping it for structural consistency
# but recommend loading the new file:
load(file.path(out_dir, 'ev.100k_multi_samples_v8.rda'))

# Collect all EV lists into a single named list for easier iteration
all_ev_lists <- list(
    RBL = EV.rbl,
    LCL = EV.lcl,
    GCBC = EV.gcbc,
    MBC = EV.mbc,
    NBC = EV.nbc,
    PC = EV.pc
)
sample_names <- names(all_ev_lists)


# generate compartment a and b figures
chroms <- paste0('chr',c(1:22,'X'))
for(chr in chroms){
    # Adjust height to accommodate all 6 samples
    pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.pdf')), width=10, height=12)
    par(mfrow=c(6,1), font.lab=2, cex.lab=1.2, lty = 0, mar=c(0,4,0,0), oma = c(0,0,0,0), mgp = c(2, 1, 0), xaxs='i', yaxs='i')

    barplot(EV.rbl[[chr]], col = ifelse(EV.rbl[[chr]] > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.12,0.12))
    barplot(EV.lcl[[chr]], col = ifelse(EV.lcl[[chr]] > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.12,0.12))
    barplot(EV.gcbc[[chr]], col = ifelse(EV.gcbc[[chr]] > 0,"red","blue"),space = 0,ylab='GCBC EV',ylim=c(-0.12,0.12))
    barplot(EV.mbc[[chr]], col = ifelse(EV.mbc[[chr]] > 0,"red","blue"),space = 0,ylab='MBC EV',ylim=c(-0.12,0.12))
    barplot(EV.nbc[[chr]], col = ifelse(EV.nbc[[chr]] > 0,"red","blue"),space = 0,ylab='NBC EV',ylim=c(-0.12,0.12))
    barplot(EV.pc[[chr]], col = ifelse(EV.pc[[chr]] > 0,"red","blue"),space = 0,ylab='PC EV',ylim=c(-0.12,0.12))

    dev.off()
}

# EV histogram (for each sample individually)
pdf(file = file.path(out_dir, 'ev.hist_all_samples.pdf'), width = 9, height = length(chroms) * 2) # Adjusted height
par(mfrow=c(length(chroms), 6), font.lab=2, cex.lab=1.2) # Adjusted columns to 6
for(i in seq_along(EV.rbl)){
    hist(EV.rbl[[i]],n=50,main=paste0(names(EV.rbl)[i],': RBL'))
    hist(EV.lcl[[i]],n=50,main=paste0(names(EV.lcl)[i],': LCL'))
    hist(EV.gcbc[[i]],n=50,main=paste0(names(EV.gcbc)[i],': GCBC'))
    hist(EV.mbc[[i]],n=50,main=paste0(names(EV.mbc)[i],': MBC'))
    hist(EV.nbc[[i]],n=50,main=paste0(names(EV.nbc)[i],': NBC'))
    hist(EV.pc[[i]],n=50,main=paste0(names(EV.pc)[i],': PC'))
}
dev.off()


# EV histogram: pairwise comparison across all samples (where signs differ)
# This will generate (N * (N-1) / 2) * 3 plots per chromosome for N samples
# For 6 samples, this is (6*5/2) * 3 = 15 * 3 = 45 plots per chromosome
# Total plots: 45 * 23 chromosomes = 1035 plots. This PDF will be VERY large.
pdf(file = file.path(out_dir, 'ev.hist_pairwise_diff_signs.pdf'), width = 12, height = length(chroms) * (length(sample_names) * (length(sample_names) - 1) / 2) * 1.5)
par(mfrow=c(length(chroms), (length(sample_names) * (length(sample_names) - 1) / 2) * 3), font.lab=2, cex.lab=1.2) # Adjusting columns for pairwise plots
for(i in seq_along(chroms)){ # Loop through chromosomes
    chr <- chroms[i]
    # Pairwise comparison loops
    for (s1_idx in 1:(length(sample_names) - 1)) {
        for (s2_idx in (s1_idx + 1):length(sample_names)) {
            sample1_name <- sample_names[s1_idx]
            sample2_name <- sample_names[s2_idx]

            ev1 <- all_ev_lists[[sample1_name]][[chr]]
            ev2 <- all_ev_lists[[sample2_name]][[chr]]

            # Find indices where signs differ between the two samples
            idx_diff_sign <- sign(ev1) != sign(ev2)

            # Check if there are enough points with differing signs to plot
            if (sum(idx_diff_sign, na.rm = TRUE) > 0) {
                hist(ev1[idx_diff_sign], n=50, main=paste0(chr, ': ', sample1_name, ' (', sample2_name, ' diff sign)'))
                hist(ev2[idx_diff_sign], n=50, main=paste0(chr, ': ', sample2_name, ' (', sample1_name, ' diff sign)'))
                hist(ev2[idx_diff_sign] - ev1[idx_diff_sign], n=50, main=paste0(chr, ': ', sample2_name, '-', sample1_name, ' (diff sign)'))
            } else {
                # Plot empty histograms or a placeholder if no differing signs for that pair/chr
                plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main=paste0(chr, ': ', sample1_name, '-', sample2_name, ' (No diff signs)'))
                plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main="")
                plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main="")
            }
        }
    }
}
dev.off()


# MA plot: pairwise comparison across all samples
# Total plots: (N * (N-1) / 2) * 23 chromosomes = 15 * 23 = 345 plots.
pdf(file = file.path(out_dir, 'ev.ma_pairwise_plots.pdf'), width = 12, height = length(chroms) * (length(sample_names) * (length(sample_names) - 1) / 2)) # Adjust height
par(mfrow=c(length(chroms), (length(sample_names) * (length(sample_names) - 1) / 2)), font.lab=2, cex.lab=1.2)
for(i in seq_along(chroms)){ # Loop through chromosomes
    chr <- chroms[i]
    for (s1_idx in 1:(length(sample_names) - 1)) {
        for (s2_idx in (s1_idx + 1):length(sample_names)) {
            sample1_name <- sample_names[s1_idx]
            sample2_name <- sample_names[s2_idx]

            ev1 <- all_ev_lists[[sample1_name]][[chr]]
            ev2 <- all_ev_lists[[sample2_name]][[chr]]

            # Ensure both EV vectors are of the same length and contain valid numbers
            common_idx <- which(!is.na(ev1) & !is.na(ev2))

            if (length(common_idx) > 0) {
                EV_diff <- ev2[common_idx] - ev1[common_idx]
                EV_mean <- (ev1[common_idx] + ev2[common_idx]) / 2
                plot(EV_mean, EV_diff, pch=16, cex=0.25,
                     main=paste0(chr, ": ", sample2_name, " vs ", sample1_name),
                     xlab = "Mean EV", ylab = "EV Difference")
                abline(h=0,lty="dashed")
            } else {
                plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="",
                     main=paste0(chr, ": ", sample2_name, " vs ", sample1_name, " (No common data)"))
            }
        }
    }
}
dev.off()

# sign flipped: all samples (density plot of differences for sign-flipped bins)
# This will now include differences for ALL pairwise sample combinations where signs flip.
EV.flip_all_pairwise <- c()
for (s1_idx in 1:(length(sample_names) - 1)) {
    for (s2_idx in (s1_idx + 1):length(sample_names)) {
        sample1_name <- sample_names[s1_idx]
        sample2_name <- sample_names[s2_idx]

        for(i in seq_along(chroms)){
            chr <- chroms[i]
            ev1 <- all_ev_lists[[sample1_name]][[chr]]
            ev2 <- all_ev_lists[[sample2_name]][[chr]]

            idx <- sign(ev1) != sign(ev2)
            EV_diff_flipped <- ev2[idx] - ev1[idx]
            EV.flip_all_pairwise <- c(EV.flip_all_pairwise, EV_diff_flipped)
        }
    }
}

# The density plot will aggregate all these pairwise differences
ev_flip_density_pairwise <- density(EV.flip_all_pairwise, na.rm = TRUE)
pdf(file.path(out_dir, 'ev.flip_density_pairwise.pdf'),6,6)
plot(ev_flip_density_pairwise, main = "Density of EV Differences for Pairwise Sign-Flipped Bins (All Samples)")
abline(v=0,lty="dashed")
dev.off()


# --- Specific Chromosome Compartment Examples (fig1a, fig1b) ---
# fig1b: chr17 compartments (now showing all samples)
chr = 'chr17'
idx = (3.35e+7/1e+5 + 1):(4.9e+7/1e+5)
# Accessing the specific subsets for each sample
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
EV.gcbc.chr = EV.gcbc[[chr]][idx]
EV.mbc.chr = EV.mbc[[chr]][idx]
EV.nbc.chr = EV.nbc[[chr]][idx]
EV.pc.chr = EV.pc[[chr]][idx]

pdf(file.path(out_dir, paste0(chr,'.tad.ab.examp.pdf')),width=6,height=9) # Adjusted height for 6 plots
par(mfrow=c(6,1),font.lab=2,cex.lab=1.2,lty = 0,mar=c(0,4,0,0),oma = c(0,0,0,0),mgp = c(2, 1, 0),xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col = ifelse(EV.rbl.chr > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.05,0.05), main = paste0(chr, " Compartments (Subset)"))
barplot(EV.lcl.chr, col = ifelse(EV.lcl.chr > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.05,0.05))
barplot(EV.gcbc.chr, col = ifelse(EV.gcbc.chr > 0,"red","blue"),space = 0,ylab='GCBC EV',ylim=c(-0.05,0.05))
barplot(EV.mbc.chr, col = ifelse(EV.mbc.chr > 0,"red","blue"),space = 0,ylab='MBC EV',ylim=c(-0.05,0.05))
barplot(EV.nbc.chr, col = ifelse(EV.nbc.chr > 0,"red","blue"),space = 0,ylab='NBC EV',ylim=c(-0.05,0.05))
barplot(EV.pc.chr, col = ifelse(EV.pc.chr > 0,"red","blue"),space = 0,ylab='PC EV',ylim=c(-0.05,0.05))
dev.off()

# fig1a: chr3 compartments (now showing all samples)
chr = 'chr3'
idx = (1.2e+8/1e+5 + 1):(1.44e+8/1e+5)
# Accessing the specific subsets for each sample
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
EV.gcbc.chr = EV.gcbc[[chr]][idx]
EV.mbc.chr = EV.mbc[[chr]][idx]
EV.nbc.chr = EV.nbc[[chr]][idx]
EV.pc.chr = EV.pc[[chr]][idx]

pdf(file.path(out_dir, paste0(chr,'.tad.ab.examp.pdf')),width=10,height=9) # Adjusted height for 6 plots
par(mfrow=c(6,1),font.lab=2,cex.lab=1.2,lty = 0,mar=c(0,4,0,0),oma = c(0,0,0,0),mgp = c(2, 1, 0),xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col = ifelse(EV.rbl.chr > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.05,0.05), main = paste0(chr, " Compartments (Subset)"))
barplot(EV.lcl.chr, col = ifelse(EV.lcl.chr > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.05,0.05))
barplot(EV.gcbc.chr, col = ifelse(EV.gcbc[[chr]] > 0,"red","blue"),space = 0,ylab='GCBC EV',ylim=c(-0.05,0.05))
barplot(EV.mbc.chr, col = ifelse(EV.mbc[[chr]] > 0,"red","blue"),space = 0,ylab='MBC EV',ylim=c(-0.05,0.05))
barplot(EV.nbc.chr, col = ifelse(EV.nbc[[chr]] > 0,"red","blue"),space = 0,ylab='NBC EV',ylim=c(-0.05,0.05))
barplot(EV.pc.chr, col = ifelse(EV.pc[[chr]] > 0,"red","blue"),space = 0,ylab='PC EV',ylim=c(-0.05,0.05))
dev.off()
