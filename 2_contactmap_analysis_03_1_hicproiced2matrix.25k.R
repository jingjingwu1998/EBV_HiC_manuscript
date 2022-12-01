#-----------------------------------
# This script is used to convert HiC-Pro contact maps to matrix file for matrix2insulation.pl to use.
#-----------------------------------

# lcl,day28
bed <- read.table('matrix/Day28/raw/25000/Day28_25000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/Day28/iced/25000/Day28_25000_iced.matrix',sep='\t',stringsAsFactors=F)
chrs <- paste0('chr',c(1:22,'X','Y'))
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
    mat[index] <- round(icetmp[,3],2)
    idx2 <- bedtmp[,3]-bedtmp[,2] == 25000
    mat2 <- mat[idx2,idx2]
    write.table(mat2,file=paste0('day28_tmp/',chr,'.tmp.matrix'),sep='\t',quote=F)
    cat(chr,'\n')
}
system(paste0("cat <(echo -ne '\t') day28_tmp/",chr,'.tmp.matrix > day28_tmp/','sysmat.',chr,".matrix"))

# rbl day0
bed <- read.table('matrix/Day0/raw/25000/Day0_25000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/Day0/iced/25000/Day0_25000_iced.matrix',sep='\t',stringsAsFactors=F)
chrs <- paste0('chr',c(1:22,'X','Y'))
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
    mat[index] <- round(icetmp[,3],2)
    idx2 <- bedtmp[,3]-bedtmp[,2] == 25000
    mat2 <- mat[idx2,idx2]
    write.table(mat2,file=paste0('day0_tmp/',chr,'.tmp.matrix'),sep='\t',quote=F)
    cat(chr,'\n')
}
system(paste0("cat <(echo -ne '\t') day0_tmp/",chr,'.tmp.matrix > day0_tmp/','sysmat.',chr,".matrix"))
