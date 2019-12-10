### Using R for file conversion
df <- read.delim("GAPIT.Genotype.Numerical.txt", header=T)
dim(df)
df[1:6,1:6]
class(df)
df2 <- df[,-1]
rownames(df2) <- as.vector(as.matrix(df[,1]))
df2 <- as.matrix(df2)
df2-1
df2[1:6,1:6]


load("Inf_meanBLUPs.rda")

Inf_meanBLUPs <- Inf_meanBLUPs[,c(1,2,6)]
Inf_meanBLUPs.copy <- Inf_meanBLUPs
nam <- data.frame(rownames(df2))
names(nam) <- "Taxa"
CrossPheno <- merge(Inf_meanBLUPs, nam, by="Taxa")
df2.copy <- df2 
df2 <- df2[match(CrossPheno$Taxa, rownames(df2)),]

ph <- CrossPheno
data <- ph
for(i in 3:ncol(data)){
  tab <-data.frame(data[,c(1,2)], data[,i])
  colnames(tab) <- c("Taxa", "fam", colnames(data)[i])
  cat("<Phenotype>",file=paste("tassel_",colnames(data)[i],".txt", sep=""), sep="\n")
  write.table(t(data.frame(c("taxa", "factor", "data"))), paste("tassel_", colnames(data)[i],".txt", sep=""), sep="\t", append=TRUE, quote=F, row.names=F, col.names=F)
  write.table(t(data.frame(c("Taxa","fam", colnames(data)[i]))), paste("tassel_", colnames(data)[i],".txt", sep=""), append=T, sep="\t",quote=F,row.names=F,col.names=F)
  write.table(tab, paste("tassel_", colnames(data)[i],".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
}
