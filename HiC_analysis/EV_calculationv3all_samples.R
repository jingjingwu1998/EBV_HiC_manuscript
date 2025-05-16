#-----------------------------------
# This script is used to calculate EV value and identify compartment A and B of HiC at 100k resolution.
# And the generation of the compartment tracks in fig1a, fig1b and fig1c
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)
library(IRanges)  # runmean 需要

# output
out_dir <- 'figures/ab_fig_2_multi_samples'
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


# EV直方图，分染色体分别保存，防止拥挤
hist_dir <- file.path(out_dir, 'ev_hist_chrs')
dir.create(hist_dir, recursive = TRUE, showWarnings = FALSE)

for(i in seq_along(EV.rbl)){
  chr_name <- names(EV.rbl)[i]
  pdf(file = file.path(hist_dir, paste0(chr_name, '_ev_hist.pdf')), width=10, height=16)
  par(mfrow=c(4,3), font.lab=2, cex.lab=1.2)  # 12图，4行3列布局
  
  # 各样本 EV 分布
  hist(EV.rbl[[i]], n=50, main=paste0(chr_name, ': RBL'), xlab='EV')
  hist(EV.lcl[[i]], n=50, main=paste0(chr_name, ': LCL'), xlab='EV')
  hist(EV.gcbc[[i]], n=50, main=paste0(chr_name, ': GCBC'), xlab='EV')
  hist(EV.mbc[[i]], n=50, main=paste0(chr_name, ': MBC'), xlab='EV')
  hist(EV.nbc[[i]], n=50, main=paste0(chr_name, ': NBC'), xlab='EV')
  hist(EV.pc[[i]], n=50, main=paste0(chr_name, ': PC'), xlab='EV')

  # 差异直方图（样本两两比较）
  hist(EV.lcl[[i]] - EV.rbl[[i]], n=50, main=paste0(chr_name, ': LCL - RBL'), xlab='EV diff')
  hist(EV.gcbc[[i]] - EV.rbl[[i]], n=50, main=paste0(chr_name, ': GCBC - RBL'), xlab='EV diff')
  hist(EV.mbc[[i]] - EV.rbl[[i]], n=50, main=paste0(chr_name, ': MBC - RBL'), xlab='EV diff')
  hist(EV.nbc[[i]] - EV.rbl[[i]], n=50, main=paste0(chr_name, ': NBC - RBL'), xlab='EV diff')
  hist(EV.pc[[i]] - EV.rbl[[i]], n=50, main=paste0(chr_name, ': PC - RBL'), xlab='EV diff')
  hist(EV.gcbc[[i]] - EV.lcl[[i]], n=50, main=paste0(chr_name, ': GCBC - LCL'), xlab='EV diff')

  dev.off()
}


# EV符号不同区域的直方图，分三组对比
pdf(file = file.path(out_dir, 'ev.sign_diff_hist.pdf'), width=15, height=92)
par(mfrow = c(23, 18), font.lab = 2, cex.lab = 1.2)  # 23行，18列 (6样本两两比较，每组3张图)

for(i in seq_along(EV.rbl)){
  chr_name <- names(EV.rbl)[i]
  
  # 样本对列表，方便循环
  sample_pairs <- list(
    c("rbl", "lcl"),
    c("rbl", "gcbc"),
    c("rbl", "mbc"),
    c("rbl", "nbc"),
    c("rbl", "pc"),
    c("lcl", "gcbc"),
    c("lcl", "mbc"),
    c("lcl", "nbc"),
    c("lcl", "pc"),
    c("gcbc", "mbc"),
    c("gcbc", "nbc"),
    c("gcbc", "pc"),
    c("mbc", "nbc"),
    c("mbc", "pc"),
    c("nbc", "pc")
  )
  
  for(pair in sample_pairs){
    sample1 <- pair[1]
    sample2 <- pair[2]
    ev1 <- get(paste0("EV.", sample1))[[i]]
    ev2 <- get(paste0("EV.", sample2))[[i]]
    idx_diff <- sign(ev1) != sign(ev2)
    
    hist(ev1[idx_diff], n=50, main=paste0(chr_name, ": ", toupper(sample1), " (vs ", toupper(sample2), " sign diff)"), xlab="EV")
    hist(ev2[idx_diff], n=50, main=paste0(chr_name, ": ", toupper(sample2), " (vs ", toupper(sample1), " sign diff)"), xlab="EV")
    hist(ev2[idx_diff] - ev1[idx_diff], n=50, main=paste0(chr_name, ": ", toupper(sample2), " - ", toupper(sample1), " (sign diff)"), xlab="EV diff")
  }
}
dev.off()


# 全染色体所有bin EV差异核密度图，包含所有样本两两对比
EV_diff_all <- list(
  LCL_RBL = c(),
  GCBC_RBL = c(),
  GCBC_LCL = c(),
  MBC_RBL = c(),
  NBC_RBL = c(),
  PC_RBL = c(),
  MBC_LCL = c(),
  NBC_LCL = c(),
  PC_LCL = c(),
  NBC_GCBC = c(),
  PC_GCBC = c(),
  PC_MBC = c(),
  NBC_MBC = c(),
  PC_NBC = c()
)

for(i in seq_along(EV.rbl)){
  EV_diff_all$LCL_RBL <- c(EV_diff_all$LCL_RBL, EV.lcl[[i]] - EV.rbl[[i]])
  EV_diff_all$GCBC_RBL <- c(EV_diff_all$GCBC_RBL, EV.gcbc[[i]] - EV.rbl[[i]])
  EV_diff_all$GCBC_LCL <- c(EV_diff_all$GCBC_LCL, EV.gcbc[[i]] - EV.lcl[[i]])
  EV_diff_all$MBC_RBL <- c(EV_diff_all$MBC_RBL, EV.mbc[[i]] - EV.rbl[[i]])
  EV_diff_all$NBC_RBL <- c(EV_diff_all$NBC_RBL, EV.nbc[[i]] - EV.rbl[[i]])
  EV_diff_all$PC_RBL <- c(EV_diff_all$PC_RBL, EV.pc[[i]] - EV.rbl[[i]])
  EV_diff_all$MBC_LCL <- c(EV_diff_all$MBC_LCL, EV.mbc[[i]] - EV.lcl[[i]])
  EV_diff_all$NBC_LCL <- c(EV_diff_all$NBC_LCL, EV.nbc[[i]] - EV.lcl[[i]])
  EV_diff_all$PC_LCL <- c(EV_diff_all$PC_LCL, EV.pc[[i]] - EV.lcl[[i]])
  EV_diff_all$NBC_GCBC <- c(EV_diff_all$NBC_GCBC, EV.nbc[[i]] - EV.gcbc[[i]])
  EV_diff_all$PC_GCBC <- c(EV_diff_all$PC_GCBC, EV.pc[[i]] - EV.gcbc[[i]])
  EV_diff_all$PC_MBC <- c(EV_diff_all$PC_MBC, EV.pc[[i]] - EV.mbc[[i]])
  EV_diff_all$NBC_MBC <- c(EV_diff_all$NBC_MBC, EV.nbc[[i]] - EV.mbc[[i]])
  EV_diff_all$PC_NBC <- c(EV_diff_all$PC_NBC, EV.pc[[i]] - EV.nbc[[i]])
}

pdf(file = file.path(out_dir, 'ev.in_de.multi_samples.pdf'), width=15, height=10)
plot(density(EV_diff_all$LCL_RBL), main='Density of EV Difference: LCL - RBL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$GCBC_RBL), main='Density of EV Difference: GCBC - RBL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$GCBC_LCL), main='Density of EV Difference: GCBC - LCL', xlab='EV Difference')
abline(v=0, lty="dashed")

plot(density(EV_diff_all$MBC_RBL), main='Density of EV Difference: MBC - RBL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$NBC_RBL), main='Density of EV Difference: NBC - RBL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$PC_RBL), main='Density of EV Difference: PC - RBL', xlab='EV Difference')
abline(v=0, lty="dashed")

plot(density(EV_diff_all$MBC_LCL), main='Density of EV Difference: MBC - LCL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$NBC_LCL), main='Density of EV Difference: NBC - LCL', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$PC_LCL), main='Density of EV Difference: PC - LCL', xlab='EV Difference')
abline(v=0, lty="dashed")

plot(density(EV_diff_all$NBC_GCBC), main='Density of EV Difference: NBC - GCBC', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$PC_GCBC), main='Density of EV Difference: PC - GCBC', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$PC_MBC), main='Density of EV Difference: PC - MBC', xlab='EV Difference')
abline(v=0, lty="dashed")

plot(density(EV_diff_all$NBC_MBC), main='Density of EV Difference: NBC - MBC', xlab='EV Difference')
abline(v=0, lty="dashed")
plot(density(EV_diff_all$PC_NBC), main='Density of EV Difference: PC - NBC', xlab='EV Difference')
abline(v=0, lty="dashed")
dev.off()

# MA图，6个样本两两对比
pdf(file = file.path(out_dir, 'ev.ma.multi_samples.pdf'), width=18, height=138)  # 23染色体 × 15组对比 = 345图，3列布局
par(mfrow = c(23, 15), font.lab = 2, cex.lab = 1.2)

for(i in seq_along(EV.rbl)){
  chr_name <- names(EV.rbl)[i]

  sample_pairs <- list(
    c("lcl", "rbl"),
    c("gcbc", "rbl"),
    c("gcbc", "lcl"),
    c("mbc", "rbl"),
    c("nbc", "rbl"),
    c("pc", "rbl"),
    c("mbc", "lcl"),
    c("nbc", "lcl"),
    c("pc", "lcl"),
    c("nbc", "gcbc"),
    c("pc", "gcbc"),
    c("pc", "mbc"),
    c("nbc", "mbc"),
    c("pc", "nbc")
  )

  for(pair in sample_pairs){
    s1 <- pair[1]
    s2 <- pair[2]
    ev1 <- get(paste0("EV.", s1))[[i]]
    ev2 <- get(paste0("EV.", s2))[[i]]
    
    EV_diff <- ev1 - ev2
    EV_mean <- (ev1 + ev2) / 2
    
    # 计算x轴范围：中点 ± 1.5倍半宽度 = 3倍宽度
    x_mid <- mean(range(EV_mean, na.rm=TRUE))
    x_half_width <- diff(range(EV_mean, na.rm=TRUE)) / 2
    xlim_expanded <- c(x_mid - 1.5 * x_half_width, x_mid + 1.5 * x_half_width)
    
    plot(EV_mean, EV_diff, pch=16, cex=0.25,
         main=paste0(chr_name, ': ', toupper(s1), ' vs ', toupper(s2)),
         xlab='Mean EV', ylab='Diff EV', xlim=xlim_expanded)
    abline(h=0, lty="dashed")
  }
}
dev.off()

# chr17 染色体示例绘图
chr = 'chr17'
idx = as.integer((3.35e+7 / 1e+5 + 1):(4.9e+7 / 1e+5))
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
EV.gcbc.chr = EV.gcbc[[chr]][idx]
EV.mbc.chr = EV.mbc[[chr]][idx]
EV.nbc.chr = EV.nbc[[chr]][idx]
EV.pc.chr = EV.pc[[chr]][idx]

pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.examp.pdf')), width=6, height=12)  # 高度加大以容纳6个样本
par(mfrow = c(6,1), font.lab=2, cex.lab=1.2, lty=0, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col=ifelse(EV.rbl.chr > 0, "red", "blue"), space=0, ylab='RBL EV', ylim=c(-0.05, 0.05))
barplot(EV.lcl.chr, col=ifelse(EV.lcl.chr > 0, "red", "blue"), space=0, ylab='LCL EV', ylim=c(-0.05, 0.05))
barplot(EV.gcbc.chr, col=ifelse(EV.gcbc.chr > 0, "red", "blue"), space=0, ylab='GCBC EV', ylim=c(-0.05, 0.05))
barplot(EV.mbc.chr, col=ifelse(EV.mbc.chr > 0, "red", "blue"), space=0, ylab='MBC EV', ylim=c(-0.05, 0.05))
barplot(EV.nbc.chr, col=ifelse(EV.nbc.chr > 0, "red", "blue"), space=0, ylab='NBC EV', ylim=c(-0.05, 0.05))
barplot(EV.pc.chr,  col=ifelse(EV.pc.chr  > 0, "red", "blue"), space=0, ylab='PC EV',  ylim=c(-0.05, 0.05))
dev.off()

# chr3 染色体示例绘图
chr = 'chr3'
idx = as.integer((1.2e+8 / 1e+5 + 1):(1.44e+8 / 1e+5))
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
EV.gcbc.chr = EV.gcbc[[chr]][idx]
EV.mbc.chr = EV.mbc[[chr]][idx]
EV.nbc.chr = EV.nbc[[chr]][idx]
EV.pc.chr = EV.pc[[chr]][idx]

pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.examp.pdf')), width=10, height=12)
par(mfrow = c(6,1), font.lab=2, cex.lab=1.2, lty=0, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col=ifelse(EV.rbl.chr > 0, "red", "blue"), space=0, ylab='RBL EV', ylim=c(-0.05, 0.05))
barplot(EV.lcl.chr, col=ifelse(EV.lcl.chr > 0, "red", "blue"), space=0, ylab='LCL EV', ylim=c(-0.05, 0.05))
barplot(EV.gcbc.chr, col=ifelse(EV.gcbc.chr > 0, "red", "blue"), space=0, ylab='GCBC EV', ylim=c(-0.05, 0.05))
barplot(EV.mbc.chr, col=ifelse(EV.mbc.chr > 0, "red", "blue"), space=0, ylab='MBC EV', ylim=c(-0.05, 0.05))
barplot(EV.nbc.chr, col=ifelse(EV.nbc.chr > 0, "red", "blue"), space=0, ylab='NBC EV', ylim=c(-0.05, 0.05))
barplot(EV.pc.chr,  col=ifelse(EV.pc.chr  > 0, "red", "blue"), space=0, ylab='PC EV',  ylim=c(-0.05, 0.05))
dev.off()
