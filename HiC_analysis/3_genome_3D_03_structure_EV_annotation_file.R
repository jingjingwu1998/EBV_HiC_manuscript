#------------------------------------
# This script is used to generate 3D structure color annotation files,
# and also used to convert structure file to bed format for g3dtools
#------------------------------------

library(HiTC)
library(mixOmics)

# Positive eigenvector represents active chromatin (A compartment). 
# Negative eigenvector represents repressed chromatin (B compartment)

load('data/ev.100k_remove_0.rda')
file_list <- read.table("data/bedpe_list",sep="\t",header=F)

rbl_list <- file_list$V1[grep("day0",file_list$V1)]
lcl_list <- file_list$V1[grep("day28",file_list$V1)]

# format structure file
loop_list <- rbl_list
loop_list <- lcl_list

g3d_input <- data.frame()
ev_paint <- data.frame()
ev_final <- data.frame()
for (i in 1:length(loop_list)) {
  file_name_tmp <- paste0("miniMDS/structure_out/",loop_list[i],"_str.tsv")
  in_tmp <- file(file_name_tmp,"r")
  header_tmp <- readLines(in_tmp,n=3)
  in_file_tmp <- read.table(file_name_tmp,sep="\t",header=F,skip =3)
  
  out_put <- data.frame(V1=rep(header_tmp[1],nrow(in_file_tmp)),
                        V2=as.integer(as.numeric(header_tmp[3])+in_file_tmp$V1*100000),
                        #V3=as.integer(as.numeric(header_tmp[3])+in_file_tmp$V1*100000+100000),
                        V4=in_file_tmp$V2,
                        V5=in_file_tmp$V3,
                        V6=in_file_tmp$V4)

  # remove NaN
  out_put_final <- out_put[complete.cases(out_put),]
  #---------------------------------------------
  # format annotation file
  
  ev_rbl_idx <- which(names(EV.rbl) == header_tmp[1])
  ev_rbl_tmp <- as.data.frame(EV.rbl[[ev_rbl_idx]])
  
  colnames(ev_rbl_tmp) <- "EV"
  # map rows of ev_rbl and out_put_final
  ev_rbl_tmp$start <- as.integer(header_tmp[3])+100000*(as.integer(row.names(ev_rbl_tmp))-1)
  ev_rbl_tmp_out <- ev_rbl_tmp[which(ev_rbl_tmp$start %in% out_put_final$V2),]
  
  # formate final output
  ev_rbl_final <- data.frame(v1=rep(header_tmp[1],nrow(out_put_final)),
                             v2=as.integer(ev_rbl_tmp_out$start),
                             v3=as.integer(ev_rbl_tmp_out$start+100000),
                             v4=ifelse(ev_rbl_tmp_out$EV > 0 ,"#0080FF","#CCCC00"))  
  ev_rbl <- data.frame(v1=rep(header_tmp[1],nrow(out_put_final)),
                       v2=as.integer(ev_rbl_tmp_out$start),
                       v3=as.integer(ev_rbl_tmp_out$start+100000),
                       v4=paste0(header_tmp[1],"_",as.integer(ev_rbl_tmp_out$start),"_",as.integer(ev_rbl_tmp_out$start+100000)),
                       v5=ev_rbl_tmp_out$EV)
  
  #---------------------------------------------
  # combine all chr
  g3d_input <- rbind(g3d_input,out_put_final)
  ev_paint <- rbind(ev_paint,ev_rbl_final)
  ev_final <- rbind(ev_final,ev_rbl)
}

write.table(ev_final,"washU_input/day28_ev.bed",sep="\t",quote = F,
            row.names = F,col.names = F)
write.table(g3d_input,"g3d_input/day0_100k_g3d_input.bed",sep="\t",quote = F,
            row.names = F,col.names = F)
write.table(ev_paint,"washU_input/day0_100k_ev_paint.bed",sep="\t",quote = F,
            row.names = F,col.names = F)

