#------------------------------------
# This script is used to convert HiC-Pro contact map to bedpe format
#------------------------------------
library(HiCcompare)

mat <- read.table("../analysis/contact_map/matrix/Day28/iced/100000/Day28_100000_iced.matrix")
bed <- read.table("../analysis/contact_map/matrix/Day28/raw/100000/Day28_100000_abs.bed")

out <- hicpro2bedpe(mat,bed)

# merge chrom
out_final <- do.call("rbind", out$cis)

# remove 0 value chrY,chrM
out_final_2 <- out_final[which(!out_final$chr1 %in% c("chrY","chrM") & out_final$IF !=0),]

write.table(out_final,"data/day28_100k_bedpe.bed",sep="\t",row.names = F,quote = F,col.names = F)
