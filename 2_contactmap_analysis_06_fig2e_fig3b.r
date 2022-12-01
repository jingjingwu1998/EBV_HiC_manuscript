#-----------------------------------
# This script is used to generate fig2e and fig3b
#-----------------------------------
library(matrixStats)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# lcl
bed <- read.table('matrix/Day0/raw/25000/Day0_25000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/Day0/iced/25000/Day0_25000_iced.matrix',sep='\t',stringsAsFactors=F)
chrs <- paste0('chr',c(1:22,'X'))
dlr.lcl <- list()
for(chr in chrs){
    bedtmp <- bed[bed[,1]==chr,]
    mat <- matrix(NA,nrow(bedtmp),nrow(bedtmp))
    rownames(mat) <- colnames(mat) <- paste0('bin',bedtmp[,4],'|hg19|',bedtmp[,1],":",bedtmp[,2]+1,"-",bedtmp[,3]+1) 
    icetmp <- ice[ice[,1]>=min(bedtmp[,4]) & ice[,1]<=max(bedtmp[,4]) & ice[,2]>=min(bedtmp[,4]) & ice[,2]<=max(bedtmp[,4]),]
    icetmp[,3] <- round(icetmp[,3],2)
    icetmptmp <- icetmp
    icetmptmp[,1] <- icetmp[,2]
    icetmptmp[,2] <- icetmp[,1]
    icetmp <- rbind(icetmp,icetmptmp)
    index <- as.matrix(icetmp[,1:2])
    index <- index - min(bedtmp[,4]) + 1
    mat[index] <- icetmp[,3]
    dlrtmp = array(0,dim=ncol(mat))
    for(i in seq_len(ncol(mat))){ #3000000/25000 = 120
        selfi = sum(mat[i,i],na.rm=T)
        locali = sum(mat[max(1,i-60):min(i+60,ncol(mat)),i],na.rm=T) - selfi
        disti = sum(mat[,i],na.rm=T)-locali-selfi
        dlrtmp[i] = log2((disti+1)/(locali+1))
        if(i%%500==0) cat(chr,i,'\n')
    }
    dlr.lcl[[chr]] <- dlrtmp
}

# rbl
bed <- read.table('matrix/Day28/raw/25000/Day28_25000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/Day28/iced/25000/Day28_25000_iced.matrix',sep='\t',stringsAsFactors=F)
chrs <- paste0('chr',c(1:22,'X'))
dlr.rbl <- list()
for(chr in chrs){
    bedtmp <- bed[bed[,1]==chr,]
    mat <- matrix(NA,nrow(bedtmp),nrow(bedtmp))
    rownames(mat) <- colnames(mat) <- paste0('bin',bedtmp[,4],'|hg19|',bedtmp[,1],":",bedtmp[,2]+1,"-",bedtmp[,3]+1) 
    icetmp <- ice[ice[,1]>=min(bedtmp[,4]) & ice[,1]<=max(bedtmp[,4]) & ice[,2]>=min(bedtmp[,4]) & ice[,2]<=max(bedtmp[,4]),]
    icetmp[,3] <- round(icetmp[,3],2)
    icetmptmp <- icetmp
    icetmptmp[,1] <- icetmp[,2]
    icetmptmp[,2] <- icetmp[,1]
    icetmp <- rbind(icetmp,icetmptmp)
    index <- as.matrix(icetmp[,1:2])
    index <- index - min(bedtmp[,4]) + 1
    mat[index] <- icetmp[,3]
    dlrtmp = array(0,dim=ncol(mat))
    for(i in seq_len(ncol(mat))){ #3000000/25000 = 120
        selfi = sum(mat[i,i],na.rm=T)
        locali = sum(mat[max(1,i-60):min(i+60,ncol(mat)),i],na.rm=T) - selfi
        disti = sum(mat[,i],na.rm=T)-locali-selfi
        dlrtmp[i] = log2((disti+1)/(locali+1))
        if(i%%500==0) cat(chr,i,'\n')
    }
    dlr.rbl[[chr]] <- dlrtmp
}
bin = 25000
save(dlr.lcl,dlr.rbl,bin,file='dlr.25k.rda')

## dlr in gene groups
load('dlr.25k.rda')
load('genegroup.rda')
load('genes.rda')

# gene groups
avana <- genes[genes$gene_name %in% genegroup[[1]],2:7]
ese <- genes[genes$gene_name %in% genegroup[[2]],2:7]
ese <- ese[!(ese$gene_name %in% c("SNORA9",  "SNORD61", "SNORD74", "SNORD75", "SNORD77", "SNORD78", "SNORD81")),]
avana <- makeGRangesFromDataFrame(avana,keep=T)
ese <- makeGRangesFromDataFrame(ese,keep=T)
gene <- makeGRangesFromDataFrame(genes,keep=T)
proml <- 5000
avana.prom <- promoters(avana,upstream=proml, downstream=0)
ese.prom <- promoters(ese,upstream=proml, downstream=0)
seqs <- seqlengths(Hsapiens)[1:23]
prom <- promoters(gene,upstream=proml, downstream=0)
seqlevels(prom,pruning.mode='coarse') = names(seqs)

# genome DLR
bins <- sapply(seqs, function(i) length(seq(1, i, bin)))
starts <- unlist(lapply(seqs, function(i) seq(1, i, bin)))
ends <- unlist(lapply(seqs, function(i) c(seq(bin, i, bin),i)))
chrs <- rep(names(seqs), times=bins)
regions <- GRanges(chrs,IRanges(start=starts,end=ends),lcl=round(unlist(dlr.lcl),3),rbl=round(unlist(dlr.rbl),3))
fo.avana <- findOverlaps(avana.prom,regions,minoverlap=proml/2)
fo.ese <- findOverlaps(ese.prom,regions,minoverlap=proml/2)
fo.gene <- findOverlaps(prom,regions,minoverlap=proml/2)
up <- 5
down <- 10

avan.idx <- subjectHits(fo.avana)[match(seq_along(avana.prom),queryHits(fo.avana))]
avan.idx.mat <- t(sapply(seq_along(avan.idx),function(i) if(as.logical(strand(avana.prom)[i]=='+')) (avan.idx[i]-up):(avan.idx[i]+down)
                         else (avan.idx[i]+up):(avan.idx[i]-down)))
avan.idx.mat <- avan.idx.mat[rowMins(avan.idx.mat)>=1 & rowMaxs(avan.idx.mat)<=length(regions),]
filter.avan <- apply(matrix(seqnames(regions)[as.integer(avan.idx.mat)],nrow(avan.idx.mat),ncol(avan.idx.mat)),1,function(x) length(table(x))==1)
avan.idx.mat <- avan.idx.mat[filter.avan,]
avana.lcl <- matrix(regions$lcl[as.integer(avan.idx.mat)],nrow(avan.idx.mat),ncol(avan.idx.mat))
avana.rbl <- matrix(regions$rbl[as.integer(avan.idx.mat)],nrow(avan.idx.mat),ncol(avan.idx.mat))
avana.lcl <- avana.lcl-rowMedians(avana.lcl)
avana.rbl <- avana.rbl-rowMedians(avana.rbl)

ese.idx <- subjectHits(fo.ese)[match(seq_along(ese.prom),queryHits(fo.ese))]
ese.idx.mat <- t(sapply(seq_along(ese.idx),function(i) if(as.logical(strand(ese.prom)[i]=='+')) (ese.idx[i]-up):(ese.idx[i]+down)
                        else (ese.idx[i]+up):(ese.idx[i]-down)))
ese.idx.mat <- ese.idx.mat[rowMins(ese.idx.mat)>=1 & rowMaxs(ese.idx.mat)<=length(regions),]
filter.ese <- apply(matrix(seqnames(regions)[as.integer(ese.idx.mat)],nrow(ese.idx.mat),ncol(ese.idx.mat)),1,function(x) length(table(x))==1)
ese.idx.mat <- ese.idx.mat[filter.ese,]
ese.lcl <- matrix(regions$lcl[as.integer(ese.idx.mat)],nrow(ese.idx.mat),ncol(ese.idx.mat))
ese.rbl <- matrix(regions$rbl[as.integer(ese.idx.mat)],nrow(ese.idx.mat),ncol(ese.idx.mat))
ese.lcl <- ese.lcl-rowMedians(ese.lcl)
ese.rbl <- ese.rbl-rowMedians(ese.rbl)

gene.idx <- subjectHits(fo.gene)[match(seq_along(prom),queryHits(fo.gene))]
gene.idx.mat <- t(sapply(seq_along(gene.idx),function(i) if(as.logical(strand(prom)[i]=='+')) (gene.idx[i]-up):(gene.idx[i]+down)
                         else (gene.idx[i]+up):(gene.idx[i]-down)))
gene.idx.mat <- gene.idx.mat[rowMins(gene.idx.mat)>=1 & rowMaxs(gene.idx.mat)<=length(regions),]
filter.gene <- apply(matrix(seqnames(regions)[as.integer(gene.idx.mat)],nrow(gene.idx.mat),ncol(gene.idx.mat)),1,function(x) length(table(x))==1)
gene.idx.mat <- gene.idx.mat[filter.gene,]
gene.lcl <- matrix(regions$lcl[as.integer(gene.idx.mat)],nrow(gene.idx.mat),ncol(gene.idx.mat))
gene.rbl <- matrix(regions$rbl[as.integer(gene.idx.mat)],nrow(gene.idx.mat),ncol(gene.idx.mat))
gene.lcl <- gene.lcl-rowMedians(gene.lcl)
gene.rbl <- gene.rbl-rowMedians(gene.rbl)

# figure 2e
pdf("figures/fig2e.pdf",width = 5,height=5)
plot(colMedians(gene.lcl),ylim=c(-0.1,0.1),type="l",lty='dashed',lwd=1,
     col='#5B9BD5',main='LCL Essential Genes',
     ylab='Normalized DLR',xaxt = 'n',xlab='Distance to TSS')
lines(colMedians(gene.rbl),col='orange',lty='dashed',lwd=1)
axis(1, at=c(1,6,11,16), labels=c(-1.25e+5,0,1.25e+5,2.5e+5))
legend('bottomright',c('LCL All Genes','RBL All Genes'),col=c('#5B9BD5','orange'),lty='dashed',lwd=2,bty='n')
lines(colMedians(avana.lcl),col='#5B9BD5',type='o',lwd=1)
lines(colMedians(avana.rbl),col='orange',type='o',lwd=1)
legend('topleft',c('LCL','RBL'),col=c('#5B9BD5','orange'),lwd=1,bty='n')
dev.off()


# figure 3b
pdf("figures/fig3b.pdf",width = 5,height=5)
plot(colMedians(gene.lcl),ylim=c(-0.2,0.2),type="l",lty='dashed',lwd=1,
     col='#5B9BD5',main='EBV Super Enhancers',
     ylab='Normalized DLR',xaxt = 'n',xlab='Distance to TSS')
lines(colMedians(gene.rbl),col='orange',lty='dashed',lwd=1)
axis(1, at=c(1,6,11,16), labels=c(-1.25e+5,0,1.25e+5,2.5e+5))
legend('bottomright',c('LCL All Genes','RBL All Genes'),col=c('#5B9BD5','orange'),lty='dashed',lwd=2,bty='n')
lines(colMedians(ese.lcl),col='#5B9BD5',type='o',pch=15,lwd=1)
lines(colMedians(ese.rbl),col='orange',type='o',pch=15,lwd=1)
legend('topleft',c('LCL','RBL'),col=c('#5B9BD5','orange'),lwd=1,bty='n')
dev.off()

