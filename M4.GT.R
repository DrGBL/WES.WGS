#the following 4 lines need to be adjusted to your folders

pathGT<-"/path/to/M4/M4.GT.chr"
pathGTExclusion<-"/path/to/M4.2/M4.2.GT.chr"
pathCol<-"/path/to/columns.txt"
pathLong<-"/path/to/M4/M4.long.chr"
pathOut<-"/path/to/M4/finalMask/"
pathOutExclusion<-"/path/to/M4.2/finalMask/"



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
  #the masks with gnomad/esp only
  M4.GT <- read.table(paste0(pathGT,i,".txt"), quote="\"", comment.char="")
  colnames(M4.GT)[2:ncol(M4.GT)]<-columns[10:length(columns)]
  colnames(M4.GT)[1]<-"variants"
  
  M4.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M4.long)<-c("gene", "variants")
  
  #tmp M4 file to work on to merge genotypes
  M4.tmp<-merge(M4.long, M4.GT)
  M4.tmp$gene<-droplevels(M4.tmp$gene)
  nVar<-which(colnames(M4.tmp)=="variants")
  M4.tmp<-M4.tmp[,-nVar]

  M4.tmp<-M4.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M4.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M4.tmp$ID))
  M4<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M4.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M4)[1]<-"#CHROM"
  M4<-merge(M4,M4.tmp)
  M4<-M4 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M4, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M4<-M4[which(rows.to.keep),]
  
  write.table(M4, paste0(pathOut,"M4.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

  
  #the masks with the merged exclusion list too
  
  M4.2.GT <- read.table(paste0(pathGTExclusion,i,".txt"), quote="\"", comment.char="")
  colnames(M4.2.GT)[2:ncol(M4.2.GT)]<-columns[10:length(columns)]
  colnames(M4.2.GT)[1]<-"variants"
  
  M4.2.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M4.2.long)<-c("gene", "variants")
  
  #tmp M4.2 file to work on to merge genotypes
  M4.2.tmp<-merge(M4.2.long, M4.2.GT)
  M4.2.tmp$gene<-droplevels(M4.2.tmp$gene)
  nVar<-which(colnames(M4.2.tmp)=="variants")
  M4.2.tmp<-M4.2.tmp[,-nVar]

  M4.2.tmp<-M4.2.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M4.2.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M4.2.tmp$ID))
  M4.2<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M4.2.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M4.2)[1]<-"#CHROM"
  M4.2<-merge(M4.2,M4.2.tmp)
  M4.2<-M4.2 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M4.2, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M4.2<-M4.2[which(rows.to.keep),]
  
  write.table(M4.2, paste0(pathOutExclusion,"M4.2.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

}
