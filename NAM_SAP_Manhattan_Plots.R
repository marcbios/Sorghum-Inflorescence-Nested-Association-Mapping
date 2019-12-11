rm(list=ls())
get_genome_pos <- function( chr_num, position, buffer=12000000){
  # chr_num: vector with chromosome numbers
  # position: vector with positions, corresponding to each chr_num
  # buffer: space between chromosomes in bp
  # Chromosome lengths for sorghum:
  chr_length <- c(
    chromosome_1  = 80884392,
    chromosome_2  = 77742459,
    chromosome_3  = 74386277,  
    chromosome_4  = 68658214,
    chromosome_5  = 71854669,
    chromosome_6  = 61277060,
    chromosome_7  = 65505356,
    chromosome_8  = 62686529,
    chromosome_9  = 59416394,
    chromosome_10 = 61233695
  )
  chr_length <- chr_length + buffer
  
  position + sapply( chr_num, function(x) sum(c(0,chr_length[-length(chr_length)])[ 1:x ]))
}


chromo <- c(1,2,3,4,5,6,7,8,9,10)
posi <- c(80884392, 77742459, 74386277, 68658214, 71854669, 61277060, 65505356, 62686529, 59416394, 61233695)
chrendpos <- data.frame(chromo, posi)
head(chrendpos)
names(chrendpos) <- c("chr", "pos")
chrendpos$genome_pos <- get_genome_pos(chrendpos$chr, chrendpos$pos)

SAP_GLM_LPBL <- read.table("~/Desktop/GAPIT..Branch_length.GWAS.Results_GLM.csv", sep=",", header=T, stringsAsFactor=F)
names(SAP_GLM_LPBL)
SAP_GLM_LPBL <- SAP_GLM_LPBL[, c(1:5)]
names(SAP_GLM_LPBL) <- c("snp", "chr", "pos", "pvalue", "maf")
head(SAP_GLM_LPBL)
str(SAP_GLM_LPBL)
SAP_GLM_LPBL <- SAP_GLM_LPBL[order(SAP_GLM_LPBL$chr, SAP_GLM_LPBL$pos),]
head(SAP_GLM_LPBL)
rownames(SAP_GLM_LPBL) <- NULL
SAP_GLM_LPBL$genome_pos <- get_genome_pos(SAP_GLM_LPBL$chr, SAP_GLM_LPBL$pos)

mean_genome_pos <- tapply(SAP_GLM_LPBL$genome_pos, SAP_GLM_LPBL$chr, function(x) (min(x)+max(x))/2 )
vek<-as.numeric(SAP_GLM_LPBL$chr)%%2


SAP_LPBL_CMLM <- read.table("~/Desktop/GAPIT..Branch_length.GWAS.Results_CMLM.csv", sep=",", header=T, stringsAsFactor=F)
head(SAP_LPBL_CMLM)
names(SAP_LPBL_CMLM)
SAP_LPBL_CMLM <- SAP_LPBL_CMLM[,c(1:5)]
names(SAP_LPBL_CMLM) <- c("snp", "chr", "pos", "pvalue", "maf")
SAP_LPBL_CMLM <- SAP_LPBL_CMLM[order(SAP_LPBL_CMLM$chr, SAP_LPBL_CMLM$pos),]
head(SAP_LPBL_CMLM)
SAP_LPBL_CMLM$genome_pos <- get_genome_pos(SAP_LPBL_CMLM$chr, SAP_LPBL_CMLM$pos)
head(SAP_LPBL_CMLM)
tail(SAP_LPBL_CMLM)

# color position
vek1<-as.numeric(SAP_LPBL_CMLM$chr)%%2


LBL_JL1 <- read.table("~/Desktop/NAM/tassel_LPBL_BLUP_stepwisePerm100nMnwF1.txt", skip=2, header=T, sep="\t", stringsAsFactors=F)
LBL_JL1 <- LBL_JL1[,c(2:4,9)]
names(LBL_JL1)<- c("snp", "chr", "pos", "pvalue")
LBL_JL1 <- LBL_JL1[-nrow(LBL_JL1),]

LBL_JL2 <- read.table("~/Desktop/NAM/tassel_LPBL4BLUP_stepwisePerm1001.txt", skip=2, header=T, sep="\t", stringsAsFactors=F)
LBL_JL2 <- LBL_JL2[,c(2:4,9)]
names(LBL_JL2)<- c("snp", "chr", "pos", "pvalue")
LBL_JL2 <- LBL_JL2[-nrow(LBL_JL2),]

LBL_JL <- rbind(LBL_JL1, LBL_JL2)

head(LBL_JL)
tail(LBL_JL)

for(i in 2:ncol(LBL_JL)){
  LBL_JL[,i] <- as.numeric(LBL_JL[,i])
  
}

LBL_JL <- LBL_JL[order(LBL_JL$chr, LBL_JL$pos),]
head(LBL_JL)
rownames(LBL_JL) <- NULL

LBL_JL <- LBL_JL[order(LBL_JL$chr, LBL_JL$pos),]

LBL_JL$genome_pos <- get_genome_pos(LBL_JL$chr, LBL_JL$pos)
head(LBL_JL)
LBL_JL$logpvalue <- -log10(LBL_JL$pvalue)
tail(LBL_JL)
LBL_JL$logpvalue[which(LBL_JL$logpvalue==Inf)] <- 0

vek2<-as.numeric(LBL_JL$chr)%%2

QTL_SNPs.GLM.match <- read.csv("/Users/omo/Desktop/NAM/GLM_NAM_JL.LBL.top5percent.07172018.csv", header=T)
QTL_SNPs.GLM.match$Chr <- gsub("_.+$", "", QTL_SNPs.GLM.match$NAM_SNP)
QTL_SNPs.GLM.match$Chr <- gsub("S", "", QTL_SNPs.GLM.match$Chr)
QTL_SNPs.GLM.match$Pos <- gsub("^.+_", "", QTL_SNPs.GLM.match$NAM_SNP)
QTL_SNPs.GLM.match$Chr <- as.numeric(QTL_SNPs.GLM.match$Chr)
QTL_SNPs.GLM.match$Pos <- as.numeric(QTL_SNPs.GLM.match$Pos)
QTL_SNPs.GLM.match$genome_pos <- get_genome_pos(QTL_SNPs.GLM.match$Chr, QTL_SNPs.GLM.match$Pos)

QTL_SNPs.CLM.match <- read.csv("/Users/omo/Desktop/NAM/CLM.NAM_LPBL.top5percent.07172018.csv", header=T)
QTL_SNPs.CLM.match$Chr <- gsub("_.+$", "", QTL_SNPs.CLM.match$NAM_SNP)
QTL_SNPs.CLM.match$Chr <- gsub("S", "", QTL_SNPs.CLM.match$Chr)
QTL_SNPs.CLM.match$Pos <- gsub("^.+_", "", QTL_SNPs.CLM.match$NAM_SNP)
QTL_SNPs.CLM.match$Chr <- as.numeric(QTL_SNPs.CLM.match$Chr)
QTL_SNPs.CLM.match$Pos <- as.numeric(QTL_SNPs.CLM.match$Pos)
QTL_SNPs.CLM.match$genome_pos <- get_genome_pos(QTL_SNPs.CLM.match$Chr, QTL_SNPs.CLM.match$Pos)

#pdf("/Users/omo/Desktop/NAM/threemanhatty08232018_new2.pdf",width = 15, height = 6)
pdf("/Users/omo/Desktop/NAM/threemanhatty.LBL.11242019.2.pdf",width = 15, height = 6)

par(mar=c(6,6,2,2))
cexp <- 1.5
cexp2 <- 1.5
cexps <- 1.7
cexmain <- 2
plot( I(-log10(SAP_GLM_LPBL$pvalue)) ~ SAP_GLM_LPBL$genome_pos, cex=cexp, col=c("gray45","gray45")[as.factor(vek)], 
      ylab=expression(paste('-log'[10],(italic(p)))), cex.lab=cexps, cex.main=cexmain, cex.axis=cexps, pch=20,
      xlab='Chromosome', xaxt='n', main="", frame.plot =F, ylim=c(0,26))

#"chartreuse","chartreuse"
points( I(-log10(SAP_LPBL_CMLM$pvalue)) ~ SAP_LPBL_CMLM$genome_pos, col=c("darkgoldenrod2","darkgoldenrod2")[as.factor(vek1)], 
        pch=20, cex=cexp)#"chartreuse4","chartreuse3"

abline(v=c(QTL_SNPs.GLM.match$genome_pos), lty=2, lwd=3.0, col="purple")#adjustcolor("brown", alpha=0.3)
abline(v=c(QTL_SNPs.CLM.match$genome_pos), lty=3, lwd=3.0, col="blue")#adjustcolor("darkmagenta", alpha=0.3)

points( I(-log10(LBL_JL$pvalue)) ~ LBL_JL$genome_pos, col=c("black", "black")[as.factor(vek2)], 
        pch=21, cex=cexp2, bg=c("coral2", "coral2")[as.factor(vek2)])#firebrick4,"magenta","mediumvioletred"

axis(side=1, at=mean_genome_pos, labels=1:10, tick=F, cex.axis=cexps)




#S3_70190339. S3_70190339
#274817190

dev.off()

###########################


SAP_GLM_RL <- read.table("~/Desktop/GAPIT..Rachis_length.GWAS.Results_GLM.csv", sep=",", header=T, stringsAsFactor=F)
names(SAP_GLM_RL)
SAP_GLM_RL <- SAP_GLM_RL[, c(1:5)]
names(SAP_GLM_RL) <- c("snp", "chr", "pos", "pvalue", "maf")
head(SAP_GLM_RL)
str(SAP_GLM_RL)
SAP_GLM_RL <- SAP_GLM_RL[order(SAP_GLM_RL$chr, SAP_GLM_RL$pos),]
head(SAP_GLM_RL)
rownames(SAP_GLM_RL) <- NULL
SAP_GLM_RL$genome_pos <- get_genome_pos(SAP_GLM_RL$chr, SAP_GLM_RL$pos)

mean_genome_pos <- tapply(SAP_GLM_RL$genome_pos, SAP_GLM_RL$chr, function(x) (min(x)+max(x))/2 )
vek<-as.numeric(SAP_GLM_RL$chr)%%2


SAP_RL_CMLM <- read.table("~/Desktop/GAPIT..Rachis_length.GWAS.Results_CMLM.csv", sep=",", header=T, stringsAsFactor=F)
head(SAP_RL_CMLM)
names(SAP_RL_CMLM)
SAP_RL_CMLM <- SAP_RL_CMLM[,c(1:5)]
names(SAP_RL_CMLM) <- c("snp", "chr", "pos", "pvalue", "maf")
SAP_RL_CMLM <- SAP_RL_CMLM[order(SAP_RL_CMLM$chr, SAP_RL_CMLM$pos),]
head(SAP_RL_CMLM)
SAP_RL_CMLM$genome_pos <- get_genome_pos(SAP_RL_CMLM$chr, SAP_RL_CMLM$pos)
head(SAP_RL_CMLM)
tail(SAP_RL_CMLM)

# color position
vek2<-as.numeric(SAP_RL_CMLM$chr)%%2


RL_JL1 <- read.table("~/Desktop/NAM/tassel_RL_BLUP_stepwisePerm100nMnwF1.txt", skip=2, header=T, sep="\t", stringsAsFactors=F)
RL_JL1 <- RL_JL1[,c(2:4,9)]
names(RL_JL1)<- c("snp", "chr", "pos", "pvalue")
RL_JL1 <- RL_JL1[-nrow(RL_JL1),]

RL_JL2 <- read.table("~/Desktop/NAM/tassel_RL_BLUP_stepwisePerm1001.txt", skip=2, header=T, sep="\t", stringsAsFactors=F)
RL_JL2 <- RL_JL2[,c(2:4,9)]
names(RL_JL2)<- c("snp", "chr", "pos", "pvalue")
RL_JL2 <- RL_JL2[-nrow(RL_JL2),]

RL_JL <- rbind(RL_JL1, RL_JL2)

head(RL_JL)
tail(RL_JL)

for(i in 2:ncol(RL_JL)){
  RL_JL[,i] <- as.numeric(RL_JL[,i])
  
}

RL_JL <- RL_JL[order(RL_JL$chr, RL_JL$pos),]
head(RL_JL)
rownames(RL_JL) <- NULL

RL_JL <- RL_JL[order(RL_JL$chr, RL_JL$pos),]

RL_JL$genome_pos <- get_genome_pos(RL_JL$chr, RL_JL$pos)
head(RL_JL)
RL_JL$logpvalue <- -log10(RL_JL$pvalue)
tail(RL_JL)
RL_JL$logpvalue[which(RL_JL$logpvalue==Inf)] <- 0


vek3<-as.numeric(RL_JL$chr)%%2

QTL_SNPs.GLM.match <- read.csv("/Users/omo/Desktop/NAM/GLM.NAM.RL.top5percent.csv", header=T)
QTL_SNPs.GLM.match$Chr <- gsub("_.+$", "", QTL_SNPs.GLM.match$NAM_SNP)
QTL_SNPs.GLM.match$Chr <- gsub("S", "", QTL_SNPs.GLM.match$Chr)
QTL_SNPs.GLM.match$Pos <- gsub("^.+_", "", QTL_SNPs.GLM.match$NAM_SNP)
QTL_SNPs.GLM.match$Chr <- as.numeric(QTL_SNPs.GLM.match$Chr)
QTL_SNPs.GLM.match$Pos <- as.numeric(QTL_SNPs.GLM.match$Pos)
QTL_SNPs.GLM.match$genome_pos <- get_genome_pos(QTL_SNPs.GLM.match$Chr, QTL_SNPs.GLM.match$Pos)

QTL_SNPs.CLM.match <- read.csv("/Users/omo/Desktop/NAM/CLM.NAM.RL.top5percent.csv", header=T)
QTL_SNPs.CLM.match$Chr <- gsub("_.+$", "", QTL_SNPs.CLM.match$NAM_SNP)
QTL_SNPs.CLM.match$Chr <- gsub("S", "", QTL_SNPs.CLM.match$Chr)
QTL_SNPs.CLM.match$Pos <- gsub("^.+_", "", QTL_SNPs.CLM.match$NAM_SNP)
QTL_SNPs.CLM.match$Chr <- as.numeric(QTL_SNPs.CLM.match$Chr)
QTL_SNPs.CLM.match$Pos <- as.numeric(QTL_SNPs.CLM.match$Pos)
QTL_SNPs.CLM.match$genome_pos <- get_genome_pos(QTL_SNPs.CLM.match$Chr, QTL_SNPs.CLM.match$Pos)

pdf("/Users/omo/Desktop/NAM/threemanhatty.RL.11242019.2.pdf",width = 15, height = 6)
par(mar=c(6,6,2,2))
cexp <- 1.5
cexp2 <- 1.5
cexps <- 1.7
cexmain <- 2
plot( I(-log10(SAP_GLM_RL$pvalue)) ~ SAP_GLM_RL$genome_pos, cex=cexp, col=c("gray45","gray45")[as.factor(vek)], 
      ylab=expression(paste('-log'[10],(italic(p)))), cex.lab=cexps, cex.main=cexmain, cex.axis=cexps, pch=20,
      xlab='Chromosome', xaxt='n', main="", frame.plot =F, ylim=c(0,32))


points( I(-log10(SAP_RL_CMLM$pvalue)) ~ SAP_RL_CMLM$genome_pos, col=c("darkgoldenrod2","darkgoldenrod2")[as.factor(vek2)], 
        pch=20, cex=cexp)#"chartreuse4","chartreuse3", "gold3","gold"

abline(v=c(QTL_SNPs.GLM.match$genome_pos), lty=2, lwd=3.0, col="purple")#adjustcolor("brown", alpha=0.3)#deepskyblue1
abline(v=c(QTL_SNPs.CLM.match$genome_pos), lty=3, lwd=3.0, col="blue")#adjustcolor("darkmagenta", alpha=0.3)


points( I(-log10(RL_JL$pvalue)) ~ RL_JL$genome_pos, col=c("black", "black")[as.factor(vek2)], 
        pch=21, cex=cexp2, bg=c("coral2", "coral2")[as.factor(vek2)])#firebrick4,"magenta","mediumvioletred"

axis(side=1, at=mean_genome_pos, labels=1:10, tick=F, cex.axis=cexps)



dev.off()

#S3_70190156






