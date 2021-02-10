#the following 4 lines need to be adjusted to your folders

pathGT<-"/path/to/M2/M2.GT.chr"
pathGTExclusion<-"/path/to/M2.2/M2.2.GT.chr"
pathCol<-"/path/to/columns.txt"
pathLong<-"/path/to/M2/M2.long.chr"
pathOut<-"/path/to/M2/finalMask/"
pathOutExclusion<-"/path/to/M2.2/finalMask/"



library(tidyverse)

#funtion to merge the genotypes
mergeGT<-function(GT){
  nGT<-length(GT)
  if("1/1" %in% GT | "1|1" %in% GT){
    return("1/1")
  }
  if("1/0" %in% GT | "0/1" %in% GT | "1/." %in% GT | "./1" %in% GT | "1|0" %in% GT | "0|1" %in% GT | "1|." %in% GT | ".|1" %in% GT){
    return("0/1")
  }
  return("0/0")
}

columns<-scan(pathCol, what=character())

for(i in 1:22){
  #the masks with gnomad/esp only
  M2.GT <- read.table(paste0(pathGT,i,".txt"), quote="\"", comment.char="")
  colnames(M2.GT)[2:ncol(M2.GT)]<-columns[10:length(columns)]
  colnames(M2.GT)[1]<-"variants"
  
  M2.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M2.long)<-c("gene", "variants")
  
  #tmp M2 file to work on to merge genotypes
  M2.tmp<-merge(M2.long, M2.GT)
  M2.tmp$gene<-droplevels(M2.tmp$gene)
  nVar<-which(colnames(M2.tmp)=="variants")
  M2.tmp<-M2.tmp[,-nVar]

  M2.tmp<-M2.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M2.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M2.tmp$ID))
  M2<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M2.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M2)[1]<-"#CHROM"
  M2<-merge(M2,M2.tmp)
  M2<-M2 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M2, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M2<-M2[which(rows.to.keep),]
  
  write.table(M2, paste0(pathOut,"M2.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

  
  #the masks with the merged exclusion list too
  
  M2.2.GT <- read.table(paste0(pathGTExclusion,i,".txt"), quote="\"", comment.char="")
  colnames(M2.2.GT)[2:ncol(M2.2.GT)]<-columns[10:length(columns)]
  colnames(M2.2.GT)[1]<-"variants"
  
  M2.2.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M2.2.long)<-c("gene", "variants")
  
  #tmp M2.2 file to work on to merge genotypes
  M2.2.tmp<-merge(M2.2.long, M2.2.GT)
  M2.2.tmp$gene<-droplevels(M2.2.tmp$gene)
  nVar<-which(colnames(M2.2.tmp)=="variants")
  M2.2.tmp<-M2.2.tmp[,-nVar]

  M2.2.tmp<-M2.2.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M2.2.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M2.2.tmp$ID))
  M2.2<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M2.2.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M2.2)[1]<-"#CHROM"
  M2.2<-merge(M2.2,M2.2.tmp)
  M2.2<-M2.2 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M2.2, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M2.2<-M2.2[which(rows.to.keep),]
  
  write.table(M2.2, paste0(pathOutExclusion,"M2.2.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

}
