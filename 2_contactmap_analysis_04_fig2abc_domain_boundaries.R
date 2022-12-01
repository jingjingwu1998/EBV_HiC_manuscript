#-----------------------------------
# This script is used to analysis overlaps of domain boundaries (LCL & RBL) and generate heatmap figures (fig2a and fig2b)
#-----------------------------------
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VennDiagram)
library(rtracklayer)
library(matrixStats)
library(ggplot2)
library(gplots)
library(ChIPseeker)
library(RColorBrewer)
library(pheatmap)

#overlaps of domain boundaries: LCL & RBL
chrs <- paste0('chr',c(1:22,'X'))
lclgr = rblgr = GRanges()
for(chr in chrs){
  lcl = read.table(paste0('day28_tmp/sysmat.',chr,'.is500001.ids200001.insulation.boundaries.bed'),sep='\t',skip=1)
  rbl = read.table(paste0('day0_tmp/sysmat.',chr,'.is500001.ids200001.insulation.boundaries.bed'),sep='\t',skip=1)
  lclsplit = t(sapply(strsplit(as.character(lcl$V4), split="(:|-)"),function(x) as.integer(x[2:3])))
  rblsplit = t(sapply(strsplit(as.character(rbl$V4), split="(:|-)"),function(x) as.integer(x[2:3])))
  ## m2: srounding 3 bins
  lclgr = c(lclgr,GRanges(lcl$V1,IRanges(lcl$V2,lcl$V3)))
  rblgr = c(rblgr,GRanges(rbl$V1,IRanges(rbl$V2,rbl$V3)))
  cat(chr,'\n')
}
seqlengths(lclgr) <- seqlengths(rblgr) <- seqlengths(Hsapiens)[seqlevels(lclgr)]
fo = findOverlaps(lclgr,rblgr)
cat(length(rblgr),length(lclgr),length(unique(subjectHits(fo))),length(unique(queryHits(fo))),'\n')

# fig2a
pdf('figures/fig2A.pdf',4,4)
plot.new()
draw.pairwise.venn(length(rblgr),length(lclgr),min(length(unique(subjectHits(fo))),
                                                   length(unique(queryHits(fo)))),c('RBL','LCL'),main='Contact Domain Boundaries')
title('Contact Domain Boundaries')
dev.off()

# fig2b up (contact domains & domain bounrdaries: LCL)
lclcd = gaps(lclgr)
lclcd = lclcd[strand(lclcd)=='*' & start(lclcd)!=1 & !(paste0(seqnames(lclcd),'_',end(lclcd)) %in% paste0(seqlevels(lclcd),'_',seqlengths(lclcd)))]

meta = read.csv('../mt/meta.csv',stringsAsFactors = F,header=F)
lclcd.left = lclcd.right = lclcd
end(lclcd.left) = start(lclcd.left)
start(lclcd.left) = end(lclcd.left) -75000
start(lclcd.right) = end(lclcd.right)
end(lclcd.right) = start(lclcd.right) + 75000
ctcf = read.table(gzfile('../mt/ENCFF551KNM.bed.gz'),sep='\t',stringsAsFactors=F)

ctcf = sort(sortSeqlevels(GRanges(ctcf[,1],IRanges(ctcf[,2],ctcf[,3]))))
fo.ctcf.left = findOverlaps(ctcf,lclcd.left)
fo.ctcf.right = findOverlaps(ctcf,lclcd.right)
size = 4000
for(i in seq_along(lclcd.left)){
  if(i %in% subjectHits(fo.ctcf.left)){
    y <- max(queryHits(fo.ctcf.left)[!is.na(match(subjectHits(fo.ctcf.left),i))])
    tmp.ctcf <- ctcf[y]
    start(tmp.ctcf) = start(tmp.ctcf) + round(width(tmp.ctcf)/2) - size/2
    end(tmp.ctcf) <- start(tmp.ctcf) + size
    lclcd.left[i] <- tmp.ctcf
  }else{
    start(lclcd.left[i]) =  end(lclcd.left[i]) - size
  }
}
for(i in seq_along(lclcd.right)){
  if(i %in% subjectHits(fo.ctcf.right)){
    y <- min(queryHits(fo.ctcf.right)[!is.na(match(subjectHits(fo.ctcf.right),i))])
    tmp.ctcf <- ctcf[y]
    start(tmp.ctcf) = start(tmp.ctcf) + round(width(tmp.ctcf)/2) - size/2
    end(tmp.ctcf) <- start(tmp.ctcf) + size
    lclcd.right[i] <- tmp.ctcf
  }else{
    end(lclcd.right[i]) =  start(lclcd.right[i]) + size
  }
}
binnum = 40

for(i in c(3,4,6,8)){
  wig = paste0(meta[i,2],'/',meta[i,4])
  lclcd.left.mat <- summary(BigWigFile(wig),lclcd.left,size=binnum,type="mean",defaultValue=0,as='matrix')
  lclcd.right.mat = summary(BigWigFile(wig),lclcd.right,size=binnum,type="mean",defaultValue=0,as='matrix')
  breaks = quantile(cbind(lclcd.left.mat,lclcd.right.mat),seq(0.5,0.9,0.1))
  pheatmap(lclcd.left.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
           breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'left',sep='.'),width=1,height=6,
           filename=paste('LCL.CD',meta[i,1],meta[i,2],'left.pdf',sep='.'),fontsize=3)
  pheatmap(lclcd.right.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
           breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'right',sep='.'),width=1,height=6,
           filename=paste('LCL.CD',meta[i,1],meta[i,2],'right.pdf',sep='.'),fontsize=3)
}


# fig2b bottom ( contact domains & domain bounrdaries: rbl)
rblcd = gaps(rblgr)
rblcd = rblcd[strand(rblcd)=='*' & start(rblcd)!=1 & !(paste0(seqnames(rblcd),'_',end(rblcd)) %in% paste0(seqlevels(rblcd),'_',seqlengths(rblcd)))]

meta = read.csv('../mt/meta.csv',stringsAsFactors = F,header=F)
rblcd.left = rblcd.right = rblcd
end(rblcd.left) = start(rblcd.left)
start(rblcd.left) = end(rblcd.left) -75000
start(rblcd.right) = end(rblcd.right)
end(rblcd.right) = start(rblcd.right) + 75000
ctcf = read.table(gzfile('../mt/ENCFF209GQK.bed.gz'),sep='\t',stringsAsFactors=F)

ctcf = sort(sortSeqlevels(GRanges(ctcf[,1],IRanges(ctcf[,2],ctcf[,3]))))
fo.ctcf.left = findOverlaps(ctcf,rblcd.left)
fo.ctcf.right = findOverlaps(ctcf,rblcd.right)
size = 4000
for(i in seq_along(rblcd.left)){
  if(i %in% subjectHits(fo.ctcf.left)){
    y <- max(queryHits(fo.ctcf.left)[!is.na(match(subjectHits(fo.ctcf.left),i))])
    tmp.ctcf <- ctcf[y]
    start(tmp.ctcf) = start(tmp.ctcf) + round(width(tmp.ctcf)/2) - size/2
    end(tmp.ctcf) <- start(tmp.ctcf) + size
    rblcd.left[i] <- tmp.ctcf
  }else{
    start(rblcd.left[i]) =  end(rblcd.left[i]) - size
  }
}
for(i in seq_along(rblcd.right)){
  if(i %in% subjectHits(fo.ctcf.right)){
    y <- min(queryHits(fo.ctcf.right)[!is.na(match(subjectHits(fo.ctcf.right),i))])
    tmp.ctcf <- ctcf[y]
    start(tmp.ctcf) = start(tmp.ctcf) + round(width(tmp.ctcf)/2) - size/2
    end(tmp.ctcf) <- start(tmp.ctcf) + size
    rblcd.right[i] <- tmp.ctcf
  }else{
    end(rblcd.right[i]) =  start(rblcd.right[i]) + size
  }
}
binnum = 40
i=3
wig = paste0(meta[i,2],'/',meta[i,4])
rblcd.left.mat <- summary(BigWigFile(wig),rblcd.left,size=binnum,type="mean",defaultValue=0,as='matrix')
rblcd.right.mat = summary(BigWigFile(wig),rblcd.right,size=binnum,type="mean",defaultValue=0,as='matrix')
index <- order(rowSums(cbind(rblcd.left.mat,rblcd.right.mat)),decreasing=T)
breaks = quantile(cbind(rblcd.left.mat,rblcd.right.mat),seq(0.5,0.9,0.1))
pheatmap(rblcd.left.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
         breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'left',sep='.'),
         width=1,height=6,filename=paste('RBL.CD',meta[i,1],meta[i,2],'left.pdf',sep='.'),
         fontsize=3)
pheatmap(rblcd.right.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
         breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'right',sep='.'),width=1,height=6,filename=paste('RBL.CD',meta[i,1],meta[i,2],'right.pdf',sep='.'),fontsize=3)
i=4
wig = paste0(meta[i,2],'/',meta[i,4])
rblcd.left.mat <- summary(BigWigFile(wig),rblcd.left,size=binnum,type="mean",defaultValue=0,as='matrix')
rblcd.right.mat = summary(BigWigFile(wig),rblcd.right,size=binnum,type="mean",defaultValue=0,as='matrix')
breaks = quantile(cbind(rblcd.left.mat,rblcd.right.mat),seq(0.5,0.9,0.1))
pheatmap(rblcd.left.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
         breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'left',sep='.'),width=1,height=6,filename=paste('RBL.CD',meta[i,1],meta[i,2],'left.pdf',sep='.'),fontsize=3)
pheatmap(rblcd.right.mat[index,], cluster_rows = FALSE, cluster_cols=F,col=c('white',brewer.pal(3, "Oranges")),
         breaks=breaks,legend=F,main=paste(meta[i,1],meta[i,2],'right',sep='.'),width=1,height=6,filename=paste('RBL.CD',meta[i,1],meta[i,2],'right.pdf',sep='.'),fontsize=3)


