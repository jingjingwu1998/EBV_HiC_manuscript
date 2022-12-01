#-----------------------------------
# This script is used to create fig2d HiC heatmap
#-----------------------------------
library(HiTC)
library(mixOmics)

rbl <- importC('matrix/Day0/iced/25000/Day0_25000_iced.matrix',xgi.bed="matrix/Day0/raw/25000/Day0_25000_abs.bed",
               ygi.bed="matrix/Day0/raw/25000/Day0_25000_abs.bed",rm.trans=T)
lcl <- importC('matrix/Day28/iced/25000/Day28_25000_iced.matrix',xgi.bed="matrix/Day28/raw/25000/Day28_25000_abs.bed",
               ygi.bed="matrix/Day28/raw/25000/Day28_25000_abs.bed",rm.trans=T)

binsize = 25000
mapC <- function(x, y, tracks=NULL,
                  minrange=NA, maxrange=NA, w=NA, trim.range=0.98,  show.zero=FALSE, show.na=FALSE, log.data=FALSE, value=FALSE, k=500,
                  col.pos=c("white","red"), col.neg=c("blue","white"), col.na="#CCCCCC",  col.zero="#FFFFFF", grid=FALSE, title=NULL){
    if (!isBinned(x) || !isBinned(y))
        stop("x and y have to be binned to plot them on the same scale")
    if (seqlevels(x) != seqlevels(y))
        stop("x and y have to come from the same chromosome")
    ## Remore long range contact if specified
    if (!is.na(w)){
        intdata(x)[!maskdiag(intdata(x), w=w)] <- NA
        intdata(y)[!maskdiag(intdata(y), w=w)] <- NA
    }
    ## Set Graphical Environment
    setEnvDisplay(x, y, tracks=tracks, view=2)
    ## Get data to map and plots
    xdata <- getData2Map(x, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)
    ydata <- getData2Map(y, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)
    ## Plots tracks and C map
    if (!is.null(tracks)){
        if (width(range(x))>=width(range(y))){
            addImageTracks(x, tracks, orientation="h")
        }else{
            addImageTracks(y, tracks, orientation="h")
        }
    }
    par(mar=c(0,5,5,5))
    triViewC(xdata, value=value, show.zero=show.zero, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, col.zero=col.zero, title=title[1], k=k, w=w)
    par(mar=c(5,5,0,5))
    triViewC(ydata, flip=TRUE, value=value, show.zero=show.zero, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, col.zero=col.zero, title=title[2], k=k, w=w)
}

refseq <- read.table('../mt/hg19.refGene.txt',sep='\t',stringsAsFactors=F,header=F)
refseq <- unique(GRanges(refseq[,3],IRanges(refseq[,5],refseq[,6]),strand=refseq[,4],name=refseq[,2],score=0,itemRgb=NA,thick=IRanges(refseq[,7],refseq[,8])))
refseq = refseq[substr(refseq$name,1,2)=='NM']

# gene level
# chr10:124550000-125500000 BUB3
# chr1:116000000-118000000 ATP1A1 
# chr2: 1.53e+8-1.59e+8
genes <- GRanges(c('chr2','chr10','chr1'),IRanges(c(1.53e+8,124550000,116000000),c(1.59e+8,125500000,118000000)))
chrs <- as.character(seqnames(genes))
currentchr <- ""
for(i in seq_along(genes)){
    if(currentchr!=chrs[i]){
        idxrbl <- match(paste0(chrs[i],chrs[i]),names(rbl))
        idxlcl<- match(paste0(chrs[i],chrs[i]),names(lcl))
        currentchr <- chrs[i]
    }
    rblgene <- extractRegion(rbl[[idxrbl]], c(1,2), chr=chrs[i], from=start(genes)[i], to=end(genes)[i])
    lclgene <- extractRegion(lcl[[idxlcl]], c(1,2), chr=chrs[i], from=start(genes)[i], to=end(genes)[i])
    if(i==1){
        pdf(paste0('chr2.153-159mb.pdf'),6,6)
        mapC(rblgene,lclgene,maxrange=40,minrange=0.5,col.pos=c("white", "orange",'red',"black"),title=c('RBL','LCL'))
        dev.off()
    }else if(i==2){
        pdf(paste0('figures/map.white.chr10.124.55-125.5mb.pdf'),6,6)
        mapC(rblgene,lclgene,maxrange=40,col.pos=c("white", "orange", "red", "black"),title=c('RBL','LCL'))
    }else if(i==3){
        pdf(paste0('figures/map.white.chr1.116-118mb.pdf'),6,6)
        mapC(rblgene,lclgene,maxrange=60,col.pos=c("white", "orange", "red", "black"),title=c('RBL','LCL'))
    }
    dev.off()
}

