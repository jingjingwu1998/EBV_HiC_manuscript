#-----------------------------------
# This script is used to calculate EV value and identify compartment A and B of HiC at 100k resolution.
# And the generation of the compartment tracks in fig1a,fig1b and fig1c
#-----------------------------------
library(HiTC)
library(rtracklayer)
library(mixOmics)

# input 
rbl <- importC('matrix/Day0/iced/100000/Day0_100000_iced.matrix',xgi.bed="matrix/Day0/raw/100000/Day0_100000_abs.bed",ygi.bed="matrix/Day0/raw/100000/Day0_100000_abs.bed",rm.trans=T)
lcl <- importC('matrix/Day28/iced/100000/Day28_100000_iced.matrix',xgi.bed="matrix/Day28/raw/100000/Day28_100000_abs.bed",ygi.bed="matrix/Day28/raw/100000/Day28_100000_abs.bed",rm.trans=T)

binsize <- 100000

# get bigwig file
mem.maxVSize(vsize = Inf)
cov <- import('../mt/ENCFF818GNV.GM12878.H3K4Me3.bigWig',as = "Rle")

# calculate EV for each chromosome

EV.rbl <- EV.lcl <- list()
for(i in 23:1){
    ## 1st approach
    x <- normPerExpected(rbl[[i]], method="loess",stdev=TRUE)
    y <- normPerExpected(lcl[[i]], method="loess",stdev=TRUE)
    rbl.xdata <- as.matrix(intdata(forceSymmetric(x)))
    lcl.xdata <- as.matrix(intdata(forceSymmetric(y)))
    cat(sum(is.na(rbl.xdata))/length(rbl.xdata),sum(is.na(lcl.xdata))/length(lcl.xdata),'\n')
    rbl.xdata[is.na(rbl.xdata)] = 0
    lcl.xdata[is.na(lcl.xdata)] = 0
    rbl.ev <- as.numeric(runmean(Rle(pca(rbl.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))
    lcl.ev <- as.numeric(runmean(Rle(pca(lcl.xdata,ncomp = 2,center = T,scale = F,logratio = 'none')$rotation[,1]),5,endrule='constant'))
    
    covk4me3 <- log2(viewSums(Views(cov[[seqlevels(rbl[[i]])]],IRanges(seq(1,length(cov[[seqlevels(rbl[[i]])]]),binsize),width=binsize)))+1)
    idx_lcl = covk4me3!=0 & lcl.ev!=0
    idx_rbl = covk4me3!=0 & rbl.ev!=0
    
    cat(i, cor(covk4me3,lcl.ev), cor(covk4me3[idx_lcl],lcl.ev[idx_lcl]),'\n')
    cat(i, cor(covk4me3,rbl.ev), cor(covk4me3[idx_rbl],rbl.ev[idx_rbl]),'\n')
    
    if(cor(covk4me3[idx_lcl],lcl.ev[idx_lcl])<0) lcl.ev <- -lcl.ev
    if(cor(rbl.ev,lcl.ev)<0) rbl.ev <- -rbl.ev
    EV.rbl[[seqlevels(rbl[[i]])]] <- rbl.ev
    EV.lcl[[seqlevels(lcl[[i]])]] <- lcl.ev
    
    # correlation figure
    pdf(paste0('figures/correlation/',seqlevels(rbl[[i]]),"_correlation.pdf"),6,6)
    par(mfrow=c(2,2))
    plot(lcl.ev,covk4me3,pch=16,col=rgb(0,0,0,0.1),main="LCL raw correlation")
    plot(rbl.ev,covk4me3,pch=16,col=rgb(0,0,0,0.1),main="RBL raw correlation")
    plot(lcl.ev[idx_lcl],covk4me3[idx_lcl],pch=16,col=rgb(0,0,0,0.1),main="LCL remove both 0 correlation")
    plot(rbl.ev[idx_rbl],covk4me3[idx_rbl],pch=16,col=rgb(0,0,0,0.1),main="RBL remove both 0 correlation")
    #abline(v=0,lty="dashed")
    dev.off()
    
}

save(EV.rbl,EV.lcl,file='ev.100k_remove_0.rda')

load('ev.100k_remove_0.rda')

# generate compartment a and b figures
chroms <- paste0('chr',c(1:22,'X'))
for(chr in chroms){
    pdf(paste0('figures/ab_fig_2/',chr,'.tad.ab.pdf'),10,4)
    par(mfrow=c(2,1),font.lab=2,cex.lab=1.2,lty = 0,mar=c(0,4,0,0),oma = c(0,0,0,0),mgp = c(2, 1, 0),xaxs='i', yaxs='i')
    barplot(EV.rbl[[chr]], col = ifelse(EV.rbl[[chr]] > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.12,0.12))
    barplot(EV.lcl[[chr]], col = ifelse(EV.lcl[[chr]] > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.12,0.12))
    dev.off()
}

# EV histogram
pdf('figures/ab_fig_2/ev.hist.pdf',9,69)
par(mfrow=c(23,3),font.lab=2,cex.lab=1.2)
for(i in seq_along(EV.rbl)){
    hist(EV.rbl[[i]],n=50,main=paste0(names(EV.rbl)[i],': RBL'))
    hist(EV.lcl[[i]],n=50,main=paste0(names(EV.lcl)[i],': LCL'))
    hist(EV.lcl[[i]]-EV.rbl[[i]],n=50,main=paste0(names(EV.lcl)[i],': LCL-RBL'))
}
dev.off()

# EV histogram that rbl and lcl have different signs
pdf('figures/ab_fig_2/ev.hist1_2.pdf',9,69)
par(mfrow=c(23,3),font.lab=2,cex.lab=1.2)
for(i in seq_along(EV.rbl)){
    idx <- sign(EV.rbl[[i]]) != sign(EV.lcl[[i]])
    hist(EV.rbl[[i]][idx],n=50,main=paste0(names(EV.rbl)[i],': RBL'))
    hist(EV.lcl[[i]][idx],n=50,main=paste0(names(EV.lcl)[i],': LCL'))
    hist(EV.lcl[[i]][idx]-EV.rbl[[i]][idx],n=50,main=paste0(names(EV.lcl)[i],': LCL-RBL'))
}
dev.off()

# figure 1c all 100kb genomic bins increase and decrease
# increase and decrease
EV.all <- c()
for(i in seq_along(EV.rbl)){
    # if(i %in% c(9,13,19)) {
    #     EV.lcl[[i]] <- -EV.lcl[[i]]
    #     EV.rbl[[i]] <- -EV.rbl[[i]]
    # }
    EV_diff <- EV.lcl[[i]]-EV.rbl[[i]]
    EV.all <- c(EV.all,EV_diff)
}
ev_in_de_density <- density(EV.all)
pdf('ev.in_de.pdf',6,6)
plot(ev_in_de_density)
abline(v=0,lty="dashed")
dev.off()

# MA plot
pdf('ev.ma.pdf',9,69)
par(mfrow=c(23,3),font.lab=2,cex.lab=1.2)
for(i in seq_along(EV.rbl)){
    EV_diff <- EV.lcl[[i]]-EV.rbl[[i]]
    EV_mean <- (EV.rbl[[i]]+EV.lcl[[i]])/2
    plot(EV_mean,EV_diff,pch=16,main=paste0(names(EV.rbl)[i]),cex=0.25)
    abline(h=0,lty="dashed")
}
dev.off()

# sign flipped
EV.flip_all <- c()
for(i in seq_along(EV.rbl)){
    idx <- sign(EV.rbl[[i]]) != sign(EV.lcl[[i]])
    EV_diff<- EV.lcl[[i]][idx]-EV.rbl[[i]][idx]
    EV.flip_all <- c(EV.flip_all,EV_diff)
}

ev_flip_density <- density(EV.flip_all)
pdf('ev.flip.pdf',6,6)
plot(ev_flip_density)
abline(v=0,lty="dashed")
dev.off()


# fig1b: chr17 compartments
chr = 'chr17'
idx = (3.35e+7/1e+5 + 1):(4.9e+7/1e+5)
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
pdf(paste0('figures/ab_fig_2/',chr,'.tad.ab.examp.pdf'),6,3)
par(mfrow=c(2,1),font.lab=2,cex.lab=1.2,lty = 0,mar=c(0,4,0,0),oma = c(0,0,0,0),mgp = c(2, 1, 0),xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col = ifelse(EV.rbl.chr > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.05,0.05))
barplot(EV.lcl.chr, col = ifelse(EV.lcl.chr > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.05,0.05))
dev.off()

# fig1a: chr3 compartments
chr = 'chr3'
idx = (1.2e+8/1e+5 + 1):(1.44e+8/1e+5)
EV.rbl.chr = EV.rbl[[chr]][idx]
EV.lcl.chr = EV.lcl[[chr]][idx]
pdf(paste0('figures/ab_fig_2/',chr,'.tad.ab.examp.pdf'),10,3)
par(mfrow=c(2,1),font.lab=2,cex.lab=1.2,lty = 0,mar=c(0,4,0,0),oma = c(0,0,0,0),mgp = c(2, 1, 0),xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col = ifelse(EV.rbl.chr > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.05,0.05))
barplot(EV.lcl.chr, col = ifelse(EV.lcl.chr > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.05,0.05))
dev.off()
