#the following 4 lines need to be adjusted to your folders

pathGT<-"/path/to/M3/M3.GT.chr"
pathCol<-"/path/to/columns.txt"
pathLong<-"/path/to/M3/M3.long.chr"
pathOut<-"/path/to/M3/finalMask/"


library(tidyverse)

#funtion to merge the genotypes
mergeGT<-function(GT){
  nGT<-length(GT)
  if("1/1" %in% GT){
    return("1/1")
  }
  if("1/0" %in% GT | "0/1" %in% GT | "1/." %in% GT | "./1" %in% GT){
    return("0/1")
  }
  return("0/0")
}

columns<-scan(pathCol, what=character())

for(i in 1:22){
 
  M3.GT <- read.table(paste0(pathGT,i,".txt"), quote="\"", comment.char="")
  colnames(M3.GT)[2:ncol(M3.GT)]<-columns[10:length(columns)]
  colnames(M3.GT)[1]<-"variants"
  
  M3.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M3.long)<-c("gene", "variants")
  
  #tmp M3 file to work on to merge genotypes
  M3.tmp<-merge(M3.long, M3.GT)
  M3.tmp$gene<-droplevels(M3.tmp$gene)
  nVar<-which(colnames(M3.tmp)=="variants")
  M3.tmp<-M3.tmp[,-nVar]

  M3.tmp<-M3.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M3.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M3.tmp$ID))
  M3<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M3.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M3)[1]<-"#CHROM"
  M3<-merge(M3,M3.tmp)
  M3<-M3 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M3, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M3<-M3[which(rows.to.keep),]
  
  write.table(M3, paste0(pathOut,"M3.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

}
