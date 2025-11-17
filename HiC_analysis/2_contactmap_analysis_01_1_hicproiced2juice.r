#-----------------------------------
# This script is used to convert HiC-Pro contact maps results to Juicebox format (*.pairs)
# by Xiang Liu
#-----------------------------------

# day0 RBL
# setwd('/share/lab_teng/xiangliu/hic_analysis/hicpro_mapping_2/contact_maps')
# bed <- read.table('hic_results/matrix/Day0/raw/5000/Day0_5000_abs.bed',sep='\t',stringsAsFactors=F)
# ice <- read.table('hic_results/matrix/Day0/iced/5000/Day0_5000_iced.matrix',sep='\t',stringsAsFactors=F)
# out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
#                  0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
# write.table(out,file='hicpro2juice_out/Day0_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)


# Day28 LCL
# bed <- read.table('hic_results/matrix/Day28/raw/5000/Day28_5000_abs.bed',sep='\t',stringsAsFactors=F)
# ice <- read.table('hic_results/matrix/Day28/iced/5000/Day28_5000_iced.matrix',sep='\t',stringsAsFactors=F)
# out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
#                   0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
# write.table(out,file='hicpro2juice_out/Day28_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)

#-----------------------------------
# This script is used to convert HiC-Pro contact maps results to Juicebox format (*.pairs)
# by Jingjing Wu
#-----------------------------------

setwd('/share/lab_teng/trainee/JingjingWu/EBV')

# Day0 RBL
bed <- read.table('matrix/Day0/raw/5000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/Day0/iced/5000/Day0_CKDL220009103-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/RBL_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


# Day28 LCL
bed <- read.table('matrix/Day28/raw/5000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/Day28/iced/5000/Day28_CKDL220009104-1a_HN2CMDSX3_L3_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/LCL_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


# GCBC
bed <- read.table('matrix/GCBC/raw/5000/GCBC_merge_read_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/GCBC/iced/5000/GCBC_merge_read_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/GCBC_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


# MBC
bed <- read.table('matrix/MBC/raw/5000/MBC_merge_read_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/MBC/iced/5000/MBC_merge_read_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/MBC_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


# NBC
bed <- read.table('matrix/NBC/raw/5000/NBC_merge_read_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/NBC/iced/5000/NBC_merge_read_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/NBC_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


# PC
bed <- read.table('matrix/PC/raw/5000/PC_merge_read_GRCh37.bwt2pairs_5000_abs.bed',
                  sep='\t', stringsAsFactors=FALSE)
ice <- read.table('matrix/PC/iced/5000/PC_merge_read_GRCh37.bwt2pairs_5000_iced.matrix',
                  sep='\t', stringsAsFactors=FALSE)

out <- data.frame(0L,
                  bed[ice[,1],1], bed[ice[,1],2] + 1, 0L,
                  0L,
                  bed[ice[,2],1], bed[ice[,2],2] + 1, 1L,
                  round(ice[,3], 2))

write.table(out, file='hicpro2juice_out/PC_5000_iced.pairs',
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

# cMCL
bed <- read.table('matrix/cMCL/raw/5000/cMCL_5000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/cMCL/iced/5000/cMCL_5000_iced.matrix',sep='\t',stringsAsFactors=F)
out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
                  0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
write.table(out,file='hicpro2juice_out/cMCL_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)

# nnMCL
bed <- read.table('matrix/nnMCL/raw/5000/nnMCL_5000_abs.bed',sep='\t',stringsAsFactors=F)
ice <- read.table('matrix/nnMCL/iced/5000/nnMCL_5000_iced.matrix',sep='\t',stringsAsFactors=F)
out <- data.frame(0L,bed[ice[,1],1],bed[ice[,1],2]+1,0L,
                  0L,bed[ice[,2],1],bed[ice[,2],2]+1,1L,round(ice[,3],2))
write.table(out,file='hicpro2juice_out/nnMCL_5000_iced.pairs',sep='\t',row.names=F,col.names=F,quote=F)
