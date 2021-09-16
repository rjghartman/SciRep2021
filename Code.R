
# R code associated with Sex-dependent co-expression in the human body
# by R.J.G. Hartman, M. Mokry, G. Pasterkamp and H. M. den Ruijter
# published in Scientific Reports (2021)

##### Required packages #####

library(WGCNA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DESeq2)
library(pheatmap)
library(VennDiagram)


##### Loading, cleaning, and separating the GTEx data in to separate tissues #####


GTex.TPM <- read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", header = T, skip = 2)
GTex.genes <- as.character(GTex.TPM$Name)
Sample.Attr <- read.delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = T, row.names =1)
rownames(Sample.Attr) <- gsub("-","\\.",rownames(Sample.Attr))
GTex.TPM <- GTex.TPM[,-c(1,2)]
GTex.vars <- colnames(GTex.TPM)
colData <- Sample.Attr[GTex.vars,]
rm(Sample.Attr)
colData <- colData[order(colData$SMTSD),]
GTex.TPM <- GTex.TPM[,rownames(colData)]

## Adding sex, age and hardy scale to the colData frame ##

pheno <- read.delim("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", header = T, row.names = 1)
rownames(pheno) <- sub("-","\\.", rownames(pheno))
colData$phenoID <- sub("^([^\\.]*\\.[^\\.]*)\\..*$","\\1",rownames(colData))
colData$sex <- pheno[colData$phenoID,]$SEX
colData$age <- pheno[colData$phenoID,]$AGE
colData$death <- pheno[colData$phenoID,]$DTHHRDY

## We only kept samples where the RNA integrity number was higher than 6 ##

colData.rin <- subset(colData, SMRIN > 6)

GTex.TPM.rin <- GTex.TPM[,rownames(colData.rin)]

## Lastly, we only kept tissues where we at least had 80 samples of both sexes ##

tissue.selection <- names(table(colData.rin$SMTSD))[which(table(colData.rin$SMTSD, colData.rin$sex)[,2]>80&
                                                            table(colData.rin$SMTSD, colData.rin$sex)[,1]>80)]

##### Writing the cleaned GTEx TPM and column data into separate tissue subsets #####

for(i in 1:length(tissue.selection)){
  write.table(GTex.TPM.rin[,colData.rin$SMTSD == tissue.selection[i]], file = paste0(paste0("TPM ",tissue.selection[i]),".txt"))
  write.table(colData.rin[colData.rin$SMTSD == tissue.selection[i],], file = paste0(paste0("colData ",tissue.selection[i]),".txt"))
}
saveRDS(GTex.genes,"gtex genes.RDS")
rm(colData)
rm(GTex.TPM)
rm(colData.rin)
rm(GTex.TPM.rin)

##### Function for reloading the TPM and colData #####

LoadTissue <- function(tissue){
  message(paste0("Loading count and coldata for ", tissue))
  countData <- read.table(paste0("TPM ",paste0(tissue,".txt")), sep = " ")
  GTex.genes <- read.table("gtex genes.RDS")
  
  rownames(countData) <- as.character(unlist(GTex.genes))
  colData <- read.delim(paste0("colData ",paste0(tissue,".txt")), sep = " ", row.names = 1)
  colData <- colData[colnames(countData),]
  
  TissueData <- list(countData = countData,
                     colData = colData)
  return(TissueData)
}




### Note:
### had to change the Esophagus - Gastroesophageal Junction, Heart Atrial Appendage, Pancreas coldata manually, as there were extra quote signs






##### Selecting genes based on high enough expression and variance values #####

## These genes will contain the most information for co-expression analyses ##


selectingGenes <- function(tissue, geneselection, gsCutoffMean, gsCutoffVar){
  require(org.Hs.eg.db)
  require(genefilter)
  
  message("Loading count and coldata...")
  countData <- read.table(paste0("TPM ",paste0(tissue,".txt")), sep = " ", row.names = 1)
  
  GTex.genes <- read.table("gtex genes.RDS")
  
  rownames(countData) <- as.character(unlist(GTex.genes))
  colData <- read.delim(paste0("colData ",paste0(tissue,".txt")), sep = " ", row.names = 1)
  colData <- colData[colnames(countData),]
  message("Selecting genes...")
  if(geneselection == "both"){
    
    genes.M <- rownames(countData)[which(rowMeans(log2(countData[,colData$sex ==1]+0.001))>gsCutoffMean &
                                           genefilter::rowVars(log2(countData[,colData$sex == 1]+0.001))>gsCutoffVar)]
    genes.F <- rownames(countData)[which(rowMeans(log2(countData[,colData$sex ==2]+0.001))>gsCutoffMean &
                                           genefilter::rowVars(log2(countData[,colData$sex == 2]+0.001))>gsCutoffVar)]
  }
  
  else if(geneselection == "var"){
    genes.M <- rownames(countData)[which(rowVars(log2(countData[,colData$sex == 1]+0.001))>gsCutoffVar)]
    genes.F <- rownames(countData)[which(rowVars(log2(countData[,colData$sex == 2]+0.001))>gsCutoffVar)]
  }
  
  else if(geneselection == "mean"){
    genes.M <- rownames(countData)[which(rowMeans(log2(countData[,colData$sex == 1]+0.001))>gsCutoffMean)]
    genes.F <- rownames(countData)[which(rowMeans(log2(countData[,colData$sex == 2]+0.001))>gsCutoffMean)]
  }
  
  genes.intersect <- intersect(genes.M, genes.F)
  return(genes.intersect)
}

genes.intersectList <- list()
for(i in 1:24){
  message(paste0("Working on ", tissue.selection[i]))
  genes.intersectList[[i]] <- selectingGenes(tissue.selection[i],geneselection = "both",gsCutoffMean = log2(1.001), gsCutoffVar = 1)
}



##### Calculating sex-bias in 100 permutations #####

connShiftSexes_nonLog <- function(power = 5, colData, countData, genes.intersect){
  require(org.Hs.eg.db)
  require(genefilter)
  require(WGCNA)
  
  samplesID.female <- rownames(colData)[colData$sex ==2]
  
  countData.log <- log2(countData+0.001)
  countData.log.DFconn <- countData.log[genes.intersect,]
  df.conn <- NULL
  
  message(paste0("Calculating connectivities over ", paste0(floor(length(samplesID.female)/2), paste0(" samples and ", paste0(length(genes.intersect), " genes...")))))
  for(i in 1:floor(length(samplesID.female)/2)){
    
    samplesID <- rownames(colData)
    
    samplesID.male <- rownames(colData)[colData$sex ==1]
    
    samplesID.male.smaller <- sample(samplesID.male, length(samplesID.female),replace=F)
    
    SIZE <- floor(length(samplesID.female)/2)
    
    
    samplesID1 <- sample(samplesID.female, SIZE, replace = F)
    samplesID2 <- sample(samplesID.male.smaller, SIZE, replace = F)
    
    samples.conn <- c(sample(samplesID1,((length(samplesID1)+1)-i),replace = F),
                      sample(samplesID2,(i-1),replace = F))
    d.conn <- t(countData.log.DFconn[,samples.conn])
    conn.genes <- softConnectivity(d.conn, power = power, verbose = 0)
    
    df.conn <- cbind(df.conn, conn.genes)
  }
  
  return(df.conn)
}

TissueData <- list()

for(i in 1:length(tissue.selection)){
  TissueData[[i]] <- LoadTissue(tissue.selection[i])
}

megalist <- list()
for(i in 1:100){
  message(paste0("Working on iteration ", i))
  set.seed(i)
  Conns <- list()
  
  for(j in 1:24){
    #message(paste0("Working on ", tissue.selection[i]))
    Conns[[j]] <- connShiftSexes_nonLog(power = 5, TissueData[[j]][[2]],TissueData[[j]][[1]], genes.intersectList[[j]])
  }
  Conns_Quant <- lapply(Conns, function(x)apply(x, 1, function(x){
    x <- unlist(x)
    quant <- stats::quantile(1:length(x), probs = seq(0, 1, 0.2))
    return(tapply(x, findInterval(1:length(x), quant, all.inside = T), median))
  }))
  
  
  Conns_QuantFold <- lapply(Conns_Quant, function(x)apply(t(x), 1, function(x){
    return(x[5]/x[1])
  }))
  
  Conns_Genes_male <- lapply(Conns_QuantFold, function(x){
    (log2(x)>1)
  })
  Conns_Genes_female <- lapply(Conns_QuantFold, function(x){
    (log2(x)<(-1))
  })
  megalist[[i]] <- list(male = Conns_Genes_male, female = Conns_Genes_female)
}

### We wrangle the list to a form in which we have the sex-bias call per permutation per tissue

### Male-biased list 

megalist.male <- list()

for(i in 1:100){
  megalist.male[[i]] <- megalist[[i]]$male
}

#make list per tissue
tissue.list.male <- list()

for(j in 1:24){
  tissue.list.male[[j]] <- list()
}


for(j in 1:24){
  for(i in 1:100){
    tissue.list.male[[j]][[i]] <- genes.intersectList[[j]][megalist.male[[i]][[j]]]
  }
}

saveRDS(tissue.list.male,"Genelist per permutation by tissue male.RDS")

### Female-biased list

megalist.female <- list()

for(i in 1:100){
  megalist.female[[i]] <- megalist[[i]]$female
}

tissue.list.female <- list()

for(j in 1:24){
  tissue.list.female[[j]] <- list()
}


for(j in 1:24){
  for(i in 1:100){
    tissue.list.female[[j]][[i]] <- genes.intersectList[[j]][megalist.female[[i]][[j]]]
  }
}
saveRDS(tissue.list.female,"Genelist per permutation by tissue female.RDS")

## These sex-biased gene lists are added as an RDS file to the GitHub (github.com/rjghartman)

### Lastly, we only took the genes that were called sex-biased in at least 20 permutations

genes.female.20 <- list()
for(i in 1:24){
  genes.female.20[[i]] <- names(which(table(unlist(tissue.list.female[[i]]))>20))
}
names(genes.female.20) <- tissue.selection
genes.female.20 <- lapply(genes.female.20, function(x)sub("\\..*","",x))
genes.female.20 <- lapply(genes.female.20, function(x)x[!is.na(x)])
genes.female.20.symbol <- lapply(genes.female.20, function(x)mapIds(org.Hs.eg.db, x,"SYMBOL","ENSEMBL"))

genes.male.20 <- list()
for(i in 1:24){
  genes.male.20[[i]] <- names(which(table(unlist(tissue.list.male[[i]]))>20))
}
names(genes.male.20) <- tissue.selection
genes.male.20 <- lapply(genes.male.20, function(x)sub("\\..*","",x))
genes.male.20 <- lapply(genes.male.20, function(x)x[!is.na(x)])
genes.male.20.symbol <- lapply(genes.male.20, function(x)mapIds(org.Hs.eg.db, x,"SYMBOL","ENSEMBL"))


## These lists were used for subsequent enrichment analyses
