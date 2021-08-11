library(readr)
#importing data
args <- commandArgs(TRUE)
filename <-args[1]
filename1 <-args[2]
filename2 <-args[3]
filename3 <-args[4]
annovar <- read.delim(filename,header=TRUE)
eyeVarP <- read.delim(filename1,header=TRUE)
geneP <- read_csv(filename2)
annovar$combine=paste(annovar$Gene.refGene,annovar$Chr,annovar$Otherinfo5,annovar$Otherinfo7,annovar$Otherinfo8)
eyeVarP$combine=paste(eyeVarP$GENE_NAME,eyeVarP$X.CHROM,eyeVarP$POS,eyeVarP$REF,eyeVarP$ALT)
blend=merge(annovar,eyeVarP, by="combine")
blend$DamagePredCount=as.numeric(as.character(blend$DamagePredCount))
blend$Priority_Score=as.numeric(as.character(blend$Priority_Score))
blend[is.na(blend)]=0
blend$eyeVarP=blend$DamagePredCount+blend$Priority_Score
data=blend[c(2,3,4,5,6,7,8,9,10,11,156,181,13,12,14,16,24,37,38,39,40,41,42,43,44,46,48,50,52,54,56,58,60,62,63,64,65,66,67,68,79,83,84,85,86,87,88,89,90,91,92,95,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,157,158,159,160,161,162,163,164,165,166,167,168)]
data$X1000g2015aug_all=as.numeric(as.character(data$X1000g2015aug_all))
data$AF=as.numeric(as.character(data$AF))
data[is.na(data)]=0
data_fiter=data[data$'X1000g2015aug_all' <=0.01,]
data_fiter_AF=data_fiter[data_fiter$'AF' <=0.01,]
geneP$gene=paste(geneP$Gene.refGene)
data_fiter_AF$gene=paste(data_fiter_AF$Gene.refGene)
blend_gene=merge(data_fiter_AF,geneP, by="gene")
data_filter_gene=blend_gene[c(2,3,4,5,6,7,8,9,10,11,12,13,125,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123)]
write.csv(data_filter_gene,paste0(filename3,".csv"),row.names = F)
#library(readr)
#importing data
#args <- commandArgs(TRUE)
#filename <-args[1]
#filename1 <-args[2]
#filename2 <-args[3]
#filename3 <-args[4]
#annovar <- read.delim(filename,header=TRUE)
#eyeVarP <- read.delim(filename1,header=TRUE)
#hc <- read.delim(filename2,header=TRUE)
#annovar$combine=paste(annovar$Gene.refGene,annovar$Chr,annovar$Otherinfo5,annovar$Otherinfo7,annovar$Otherinfo8)
#eyeVarP$combine=paste(eyeVarP$GENE_NAME,eyeVarP$X.CHROM,eyeVarP$POS,eyeVarP$REF,eyeVarP$ALT)
#blend=merge(annovar,eyeVarP, by="combine")
#blend$DamagePredCount=as.numeric(as.character(blend$DamagePredCount))
#blend$Priority_Score=as.numeric(as.character(blend$Priority_Score))
#blend[is.na(blend)]=0
#blend$eyeVarP=blend$DamagePredCount+blend$Priority_Score
#data=blend[c(2,3,4,5,6,7,8,9,10,11,156,181,13,12,14,16,24,37,38,39,40,41,42,43,44,46,48,50,52,54,56,58,60,62,63,64,65,66,67,68,79,83,84,85,86,87,88,89,90,91,92,95,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,157,158,159,160,161,162,163,164,165,166,167,168)]
#data$X1000g2015aug_all=as.numeric(as.character(data$X1000g2015aug_all))
#data[is.na(data)]=0
#data_fiter=data[data$'X1000g2015aug_all' <=0.01,]
#data_fiter$gene=paste(data_fiter$Gene.refGene)
#hc$gene=paste(hc$Gene.refGene)
#Gene_filter=merge(data_fiter,hc, by="gene")
#eyeVarP_Gene=Gene_filter[c(2,3,4,5,6,7,8,9,10,11,12,13,125,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123)]
#write.csv(Gene_filter,paste0(filename3,".csv"),row.names = F)

