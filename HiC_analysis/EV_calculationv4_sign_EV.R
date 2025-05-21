#-----------------------------------
# This script is used to calculate EV value and identify compartment A and B of HiC at 100k resolution.
# And the generation of the compartment tracks in fig1a, fig1b and fig1c
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges)  # runmean 需要

# output
out_dir <- 'figures/ab_fig_2_multi_samples_EVsign'
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

# 手动设定RBL每条染色体EV的符号
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

# 计算EV，符号对齐，NA过滤确保相关计算安全
EV.rbl <- list()
EV.lcl <- list()
EV.gcbc <- list()
EV.mbc <- list()
EV.nbc <- list()
EV.pc <- list()

for(i in 23:1){
    # RBL处理
    x <- normPerExpected(rbl[[i]], method="loess", stdev=TRUE)
    rbl.xdata <- as.matrix(intdata(forceSymmetric(x)))
    cat("RBL NA ratio chr", seqlevels(rbl[[i]]), ":", sum(is.na(rbl.xdata))/length(rbl.xdata), "\n")
    rbl.xdata[is.na(rbl.xdata)] <- 0
    rbl.ev <- as.numeric(runmean(Rle(pca(rbl.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))
    # 加入手动符号翻转
    chr <- seqlevels(rbl[[i]])
    sign_factor <- ifelse(!is.null(manual_sign_rbl[[chr]]), manual_sign_rbl[[chr]], 1)
    rbl.ev <- rbl.ev * sign_factor

    # LCL处理
    y <- normPerExpected(lcl[[i]], method="loess", stdev=TRUE)
    lcl.xdata <- as.matrix(intdata(forceSymmetric(y)))
    cat("LCL NA ratio chr", seqlevels(lcl[[i]]), ":", sum(is.na(lcl.xdata))/length(lcl.xdata), "\n")
    lcl.xdata[is.na(lcl.xdata)] <- 0
    lcl.ev <- as.numeric(runmean(Rle(pca(lcl.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # GCBC处理
    z <- normPerExpected(gcbc[[i]], method="loess", stdev=TRUE)
    gcbc.xdata <- as.matrix(intdata(forceSymmetric(z)))
    cat("GCBC NA ratio chr", seqlevels(gcbc[[i]]), ":", sum(is.na(gcbc.xdata))/length(gcbc.xdata), "\n")
    gcbc.xdata[is.na(gcbc.xdata)] <- 0
    gcbc.ev <- as.numeric(runmean(Rle(pca(gcbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # MBC处理
    m <- normPerExpected(mbc[[i]], method="loess", stdev=TRUE)
    mbc.xdata <- as.matrix(intdata(forceSymmetric(m)))
    cat("MBC NA ratio chr", seqlevels(mbc[[i]]), ":", sum(is.na(mbc.xdata)) / length(mbc.xdata), "\n")
    mbc.xdata[is.na(mbc.xdata)] <- 0
    mbc.ev <- as.numeric(runmean(Rle(pca(mbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # NBC处理
    n <- normPerExpected(nbc[[i]], method="loess", stdev=TRUE)
    nbc.xdata <- as.matrix(intdata(forceSymmetric(n)))
    cat("NBC NA ratio chr", seqlevels(nbc[[i]]), ":", sum(is.na(nbc.xdata)) / length(nbc.xdata), "\n")
    nbc.xdata[is.na(nbc.xdata)] <- 0
    nbc.ev <- as.numeric(runmean(Rle(pca(nbc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # PC处理
    p <- normPerExpected(pc[[i]], method="loess", stdev=TRUE)
    pc.xdata <- as.matrix(intdata(forceSymmetric(p)))
    cat("PC NA ratio chr", seqlevels(pc[[i]]), ":", sum(is.na(pc.xdata)) / length(pc.xdata), "\n")
    pc.xdata[is.na(pc.xdata)] <- 0
    pc.ev <- as.numeric(runmean(Rle(pca(pc.xdata, ncomp=2, center=TRUE, scale=FALSE, logratio='none')$rotation[,1]), 5, endrule='constant'))

    # 符号对齐，以RBL为基准，逐个翻转其他样本
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


# 保存结果
save(EV.rbl, EV.lcl, EV.gcbc, EV.mbc, EV.nbc, EV.pc, file=file.path(out_dir, 'ev.100k_multi_samples.rda'))

# 生成A/B区室柱状图
chroms <- paste0('chr', c(1:22, 'X'))
for(chr in chroms){
    pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.pdf')), width=10, height=12)  # 高度适当加大以容纳6行图
    par(mfrow=c(6,1), font.lab=2, cex.lab=1.2, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')

    barplot(EV.rbl[[chr]], col=ifelse(EV.rbl[[chr]] > 0, "red", "blue"), space=0, ylab='RBL EV', ylim=c(-0.12, 0.12))
    barplot(EV.lcl[[chr]], col=ifelse(EV.lcl[[chr]] > 0, "red", "blue"), space=0, ylab='LCL EV', ylim=c(-0.12, 0.12))
    barplot(EV.gcbc[[chr]], col=ifelse(EV.gcbc[[chr]] > 0, "red", "blue"), space=0, ylab='GCBC EV', ylim=c(-0.12, 0.12))
    barplot(EV.mbc[[chr]], col=ifelse(EV.mbc[[chr]] > 0, "red", "blue"), space=0, ylab='MBC EV', ylim=c(-0.12, 0.12))
    barplot(EV.nbc[[chr]], col=ifelse(EV.nbc[[chr]] > 0, "red", "blue"), space=0, ylab='NBC EV', ylim=c(-0.12, 0.12))
    barplot(EV.pc[[chr]], col=ifelse(EV.pc[[chr]] > 0, "red", "blue"), space=0, ylab='PC EV', ylim=c(-0.12, 0.12))

    dev.off()
}
