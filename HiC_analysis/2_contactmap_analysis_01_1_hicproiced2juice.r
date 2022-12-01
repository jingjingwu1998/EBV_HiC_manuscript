#-----------------------------------
# This script is used to convert HiC-Pro contact maps results to Juicebox format (*.pairs)
#-----------------------------------

# day0 RBL
setwd('/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/contact_maps')
bed <- read.table('hic_results/matrix/Day0/raw/5000/Day0_5000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('hic_results/matrix/Day0/iced/5000/Day0_5000_iced.matrix',sep='\t',stringsAsFactors=F)
out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
                  0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
write.table(out,file='hicpro2juice_out/Day0_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)


# Day28 LCL
bed <- read.table('hic_results/matrix/Day28/raw/5000/Day28_5000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('hic_results/matrix/Day28/iced/5000/Day28_5000_iced.matrix',sep='\t',stringsAsFactors=F)
out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
                  0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
write.table(out,file='hicpro2juice_out/Day28_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)
