#this uses the output from the 1000G PCA and projected on your cohort, as done in 03.ancestryPCA.sh
#it assumes the same file names as the outputs from the previous step

#here set your working directory (the pathWork in the previous step)
setwd("/scratch/richards/guillaume.butler-laporte/WGS/1000G.PCA/")

#will use randomForest for the classifier. Need to install the package below if not done.
library(randomForest)


# read in the eigenvectors, produced in PLINK
eigenvec <- read.table('plink.eigenvec', header = FALSE, skip=0, sep = ' ')
eigenvec <- eigenvec[,2:ncol(eigenvec)]
colnames(eigenvec) <- c("Individual.ID",paste('PC', c(1:20), sep = ''))

# read in the PED data
PED <- read.table('20130606_g1k.ped', header = TRUE, skip = 0, sep = '\t')

#build data frame for random forest classifier
dataRF <- merge(eigenvec, PED[, c("Individual.ID", "Population")], all.x=TRUE)

#build plot
dataRF$Population <- factor(dataRF$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

dataRF$Continental <- rep(NA_character_, nrow(dataRF))
dataRF$Continental[which(dataRF$Population %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"))]<-"AFR"
dataRF$Continental[which(dataRF$Population %in% c("CLM","MXL","PEL","PUR"))]<-"AMR"
dataRF$Continental[which(dataRF$Population %in% c("CDX","CHB","CHS","JPT","KHV"))]<-"EAS"
dataRF$Continental[which(dataRF$Population %in% c("CEU","FIN","GBR","IBS","TSI"))]<-"EUR"
dataRF$Continental[which(dataRF$Population %in% c("BEB","GIH","ITU","PJL","STU"))]<-"SAS"
dataRF$Continental<-as.factor(dataRF$Continental)


col <- colorRampPalette(c(
  "yellow","forestgreen","grey","royalblue","black"))(length(unique(dataRF$Continental)))[factor(dataRF$Continental)]


#PCA plots, to see if you get ok results from the 1000G PCA, you should see nice clusters, if not, then the PCA didn't work
#colors are as follows
#yellow=AFR
#forestgreen=AMR
#grey=EAS
#royalblue=EUR
#black=SAS

par(mar = c(5,5,5,5), cex = 2.0,
    cex.main = 7, cex.axis = 2.75, cex.lab = 2.75, mfrow = c(1,2))

plot(dataRF[,2], dataRF[,3],
     type = 'n',
     main = 'A',
     adj = 0.5,
     xlab = 'First component',
     ylab = 'Second component',
     font = 2,
     font.lab = 2)
points(dataRF[,2], dataRF[,3], col = col, pch = 20, cex = 2.25)
legend('bottomright',
       bty = 'n',
       cex = 3.0,
       title = '',
       c('Population 1', 'Population 2', 'Population 3',
         'Population 4', 'Population 5'),
       fill = c('yellow', 'forestgreen', 'grey', 'royalblue', 'black'))

plot(dataRF[,2], dataRF[,4],
     type="n",
     main="B",
     adj=0.5,
     xlab="First component",
     ylab="Third component",
     font=2,
     font.lab=2)
points(dataRF[,2], dataRF[,4], col=col, pch=20, cex=2.25)


# here we train the random forest classifier on the 1000G reference using the first 6 PC as features
rf_classifier = randomForest(Continental ~ ., data=dataRF[which(!is.na(dataRF$Continental)),c("PC1","PC2","PC3","PC4","PC5","PC6", "Continental")], ntree=10000, importance=TRUE)

#predict ancestries in your cohort
dataPred<-dataRF[which(is.na(dataRF$Continental)),c("Individual.ID", "PC1", "PC2", "PC3","PC4","PC5","PC6")]
dataPred$Prediction<-rep(NA, nrow(dataPred))
dataPred$Prediction<-predict(rf_classifier,dataPred[,c("PC1","PC2","PC3","PC4","PC5","PC6")])

#if you want to see the results, here I provide an example with AFR ancestry
#from here, if the population is very admixted, other algorithms can be applied, or at least visual inspection can be done to remove clear outliers.
#the blank dots should cluster with the yellow dots (for AFR)
dataAfr<-dataPred %>% filter(Prediction=="AFR")
plot(dataRF[,2], dataRF[,3],
     type = 'n',
     main = 'A',
     adj = 0.5,
     xlab = 'First component',
     ylab = 'Second component',
     font = 2,
     font.lab = 2)
points(dataRF[,2], dataRF[,3], col = col, pch = 20, cex = 2.25)
points(dataAfr[,2], dataAfr[,3])

plot(dataRF[,2], dataRF[,4],
     type="n",
     main="B",
     adj=0.5,
     xlab="First component",
     ylab="Third component",
     font=2,
     font.lab=2)
points(dataRF[,2], dataRF[,4], col=col, pch=20, cex=2.25)
points(dataAfr[,2], dataAfr[,4])

# lastly write IDs per ancestry in files to be used for further analyses
write(as.character(subset(dataPred, Prediction=="AFR")$Individual.ID), "afrIDsPCA.txt")
write(as.character(subset(dataPred, Prediction=="AMR")$Individual.ID), "amrIDsPCA.txt")
write(as.character(subset(dataPred, Prediction=="EAS")$Individual.ID), "easIDsPCA.txt")
write(as.character(subset(dataPred, Prediction=="EUR")$Individual.ID), "eurIDsPCA.txt")
write(as.character(subset(dataPred, Prediction=="SAS")$Individual.ID), "sasIDsPCA.txt")

