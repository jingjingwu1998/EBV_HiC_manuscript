#-----------------------------------
# Complete script with plot_ev_bar_debug for EV plotting diagnostics
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges)  # for runmean

# output directory
out_dir <- 'figures/ab_fig_2_multi_samples_EVsign_debug'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# input
rbl <- importC(
  'matrix/Day0/iced/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/Day0/raw/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/Day0/raw/100000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)
lcl <- importC(
  'matrix/Day28/iced/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/Day28/raw/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/Day28/raw/100000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)
gcbc <- importC(
  'matrix/GCBC/iced/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/GCBC/raw/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/GCBC/raw/100000/GCBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)
mbc <- importC(
  'matrix/MBC/iced/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/MBC/raw/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/MBC/raw/100000/MBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)
nbc <- importC(
  'matrix/NBC/iced/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/NBC/raw/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/NBC/raw/100000/NBC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)
pc <- importC(
  'matrix/PC/iced/100000/PC_merge_read_GRCh37.bwt2pairs_100000_iced.matrix',
  xgi.bed = "matrix/PC/raw/100000/PC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  ygi.bed = "matrix/PC/raw/100000/PC_merge_read_GRCh37.bwt2pairs_100000_abs.bed",
  rm.trans = TRUE
)

binsize <- 100000

# manual sign control for RBL
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

# lists to store EVs
EV.rbl <- list()
EV.lcl <- list()
EV.gcbc <- list()
EV.mbc <- list()
EV.nbc <- list()
EV.pc <- list()

for(i in 23:1){
    # RBL processing
    x <- normPerExpected(rbl[[i]], method="loess", stdev=TRUE)
    rbl.xdata <- as.matrix(intdata(forceSymmetric(x)))
    cat("RBL NA ratio chr", seqlevels(rbl[[i]]), ":", sum(is.na(rbl.xdata))/length(rbl.xdata), "\n")
    rbl.xdata[is.na(rbl.xdata)] <- 0
    rbl.ev <- as.numeric(runmean(Rle(pca(rbl.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))
    # manual sign flip for RBL
    chr <- seqlevels(rbl[[i]])
    sign_factor <- ifelse(!is.null(manual_sign_rbl[[chr]]), manual_sign_rbl[[chr]], 1)
    rbl.ev <- rbl.ev * sign_factor

    # LCL processing
    y <- normPerExpected(lcl[[i]], method="loess", stdev=TRUE)
    lcl.xdata <- as.matrix(intdata(forceSymmetric(y)))
    cat("LCL NA ratio chr", seqlevels(lcl[[i]]), ":", sum(is.na(lcl.xdata))/length(lcl.xdata), "\n")
    lcl.xdata[is.na(lcl.xdata)] <- 0
    lcl.ev <- as.numeric(runmean(Rle(pca(lcl.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # GCBC processing
    z <- normPerExpected(gcbc[[i]], method="loess", stdev=TRUE)
    gcbc.xdata <- as.matrix(intdata(forceSymmetric(z)))
    cat("GCBC NA ratio chr", seqlevels(gcbc[[i]]), ":", sum(is.na(gcbc.xdata))/length(gcbc.xdata), "\n")
    gcbc.xdata[is.na(gcbc.xdata)] <- 0
    gcbc.ev <- as.numeric(runmean(Rle(pca(gcbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # MBC processing
    m <- normPerExpected(mbc[[i]], method="loess", stdev=TRUE)
    mbc.xdata <- as.matrix(intdata(forceSymmetric(m)))
    cat("MBC NA ratio chr", seqlevels(mbc[[i]]), ":", sum(is.na(mbc.xdata)) / length(mbc.xdata), "\n")
    mbc.xdata[is.na(mbc.xdata)] <- 0
    mbc.ev <- as.numeric(runmean(Rle(pca(mbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # NBC processing
    n <- normPerExpected(nbc[[i]], method="loess", stdev=TRUE)
    nbc.xdata <- as.matrix(intdata(forceSymmetric(n)))
    cat("NBC NA ratio chr", seqlevels(nbc[[i]]), ":", sum(is.na(nbc.xdata)) / length(nbc.xdata), "\n")
    nbc.xdata[is.na(nbc.xdata)] <- 0
    nbc.ev <- as.numeric(runmean(Rle(pca(nbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # PC processing
    p <- normPerExpected(pc[[i]], method="loess", stdev=TRUE)
    pc.xdata <- as.matrix(intdata(forceSymmetric(p)))
    cat("PC NA ratio chr", seqlevels(pc[[i]]), ":", sum(is.na(pc.xdata)) / length(pc.xdata), "\n")
    pc.xdata[is.na(pc.xdata)] <- 0
    pc.ev <- as.numeric(runmean(Rle(pca(pc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # Align other samples to RBL EV sign
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
}

# Debug plotting function for EV with color diagnostics
plot_ev_bar_debug <- function(ev, label) {
  ev[is.na(ev)] <- 0
  cat(paste0(label, " EV summary:\n"))
  print(summary(ev))
  print(table(sign(ev)))
  
  ev_col <- rep("gray", length(ev))
  ev_col[ev > 0] <- "red"
  ev_col[ev < 0] <- "blue"
  print(table(ev_col))
  
  barplot(ev, col=ev_col, space=0, ylab=paste(label, "EV"), ylim=c(-0.12, 0.12))
}

# Plot and save PDFs
chroms <- paste0('chr', c(1:22, 'X'))
for(chr in chroms){
  pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.debug.pdf')), width=10, height=12)
  par(mfrow=c(6,1), font.lab=2, cex.lab=1.2, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')

  plot_ev_bar_debug(EV.rbl[[chr]], "RBL")
  plot_ev_bar_debug(EV.lcl[[chr]], "LCL")
  plot_ev_bar_debug(EV.gcbc[[chr]], "GCBC")
  plot_ev_bar_debug(EV.mbc[[chr]], "MBC")
  plot_ev_bar_debug(EV.nbc[[chr]], "NBC")
  plot_ev_bar_debug(EV.pc[[chr]], "PC")

  dev.off()
}
