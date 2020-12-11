#the following 4 lines need to be adjusted to your folders

pathGT<-"/path/to/M1/M1.GT.chr"
pathGTExclusion<-"/path/to/M1.2/M1.2.GT.chr"
pathCol<-"/path/to/columns.txt"
pathLong<-"/path/to/M1/M1.long.chr"
pathOut<-"/path/to/M1/finalMask/"
pathOutExclusion<-"/path/to/M1.2/finalMask/"



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
  M1.GT <- read.table(paste0(pathGT,i,".txt"), quote="\"", comment.char="")
  colnames(M1.GT)[2:ncol(M1.GT)]<-columns[10:length(columns)]
  colnames(M1.GT)[1]<-"variants"
  
  M1.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M1.long)<-c("gene", "variants")
  
  #tmp M1 file to work on to merge genotypes
  M1.tmp<-merge(M1.long, M1.GT)
  M1.tmp$gene<-droplevels(M1.tmp$gene)
  nVar<-which(colnames(M1.tmp)=="variants")
  M1.tmp<-M1.tmp[,-nVar]

  M1.tmp<-M1.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M1.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M1.tmp$ID))
  M1<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M1.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M1)[1]<-"#CHROM"
  M1<-merge(M1,M1.tmp)
  M1<-M1 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M1, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M1<-M1[which(rows.to.keep),]
  
  write.table(M1, paste0(pathOut,"M1.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

  
  #the masks with the merged exclusion list too
  
  M1.2.GT <- read.table(paste0(pathGTExclusion,i,".txt"), quote="\"", comment.char="")
  colnames(M1.2.GT)[2:ncol(M1.2.GT)]<-columns[10:length(columns)]
  colnames(M1.2.GT)[1]<-"variants"
  
  M1.2.long<-read.table(paste0(pathLong,i,".txt"), quote="\"", comment.char="")
  colnames(M1.2.long)<-c("gene", "variants")
  
  #tmp M1.2 file to work on to merge genotypes
  M1.2.tmp<-merge(M1.2.long, M1.2.GT)
  M1.2.tmp$gene<-droplevels(M1.2.tmp$gene)
  nVar<-which(colnames(M1.2.tmp)=="variants")
  M1.2.tmp<-M1.2.tmp[,-nVar]

  M1.2.tmp<-M1.2.tmp %>% 
    group_by(gene) %>%
    summarise(across(columns[10:length(columns)],mergeGT))
  
  colnames(M1.2.tmp)[1]<-"ID"
  
  #this is the vcf
  nGene<-nrow(table(M1.2.tmp$ID))
  M1.2<-data.frame("#CHROM"=rep(paste0("chr",i), nGene), POS=c(1:nGene), ID=names(table(M1.2.tmp$ID)), REF=rep("C",nGene), ALT=rep("A",nGene), QUAL=rep(".", nGene), FILTER=rep(".", nGene), INFO=rep(".", nGene), FORMAT=rep("GT", nGene))
  colnames(M1.2)[1]<-"#CHROM"
  M1.2<-merge(M1.2,M1.2.tmp)
  M1.2<-M1.2 %>% relocate("#CHROM", "POS")
  
  #this removes rows that are all 0/0 or 1/1 since there are no variants here, this is an artefact of subsetting on ancestry
  rows.to.keep<-apply(M1.2, 1, function(r) {!(all(r == "1/1") | all(r == "0/0"))})
  M1.2<-M1.2[which(rows.to.keep),]
  
  write.table(M1.2, paste0(pathOutExclusion,"M1.2.chr",i,".txt"),row.names=FALSE,sep="\t", quote = FALSE)

}
