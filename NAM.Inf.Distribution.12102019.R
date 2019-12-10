Inf_mean.All <- read.csv("Sorghum.Inf.NAM.Distribution.Data.csv", header=T)

Inf_mean.All$Taxa <- as.character(Inf_mean.All$Taxa)
families= unique(as.character(as.vector(Inf_meanBLUPs$family_ID)))
families <- families[order(families)]
Inf_mean.All2 <- Inf_mean.All[Inf_mean.All$family_ID%in%families,]

phenames=as.character(as.vector(colnames(Inf_mean.All2[,-c(1,2)])))

#axis1 <- 'axis(side=1, at=seq(0,150,50), labels=seq(0,150,50), tick=F, cex.axis=2)'
#axis2 <- axis(side=1, at=seq(0,150,50), labels=seq(0,150,50), tick=F, cex.axis=2)
#axis3 <- axis(side=1, at=seq(100,500,100), labels=seq(100,500,100), tick=F, cex.axis=2)
#axis4 <- axis(side=1, at=seq(0,20,5), labels=seq(0,20,5), tick=F, cex.axis=2)
axis.val <- matrix(c(0,150,50,0,150,50,100,500,100,0,20,5),4,3, byrow = T)


for(i in 1:length(phenames)){
  min.val <- min(Inf_mean.All[,phenames[i]], na.rm=T)
  max.val <- max(Inf_mean.All[,phenames[i]], na.rm=T)
  pheno.trt <- data.frame(Inf_mean.All[,c(1,2)], Inf_mean.All[,phenames[i]])
  colnames(pheno.trt)[3] <- phenames[i]
  Lata.val <- pheno.trt[which(pheno.trt$Taxa=="RTx430"),3]
  
  
  pdf(paste("NAM.Inf.Distribution.Plot", phenames[i], "pdf", sep="."), 8, 11)
  par(mfcol=c(11,1))
  
  
  for(j in 1:length(families)){
    
    #, mai = c(1, 0.1, 0.1, 0.1, 0.1, 0.1)
    #par(mar=c(0,0,1,0))
    par(mar=c(0,7,4,7))
    
    pheno.trt.fam <- pheno.trt[which(pheno.trt$family_ID==families[j]),]
    
    donor.val <- pheno.trt[which(pheno.trt$Taxa==families[j]),3]
    
    dh <- density(pheno.trt.fam[,3], na.rm=T)
    y.vec <- c(dh$y)
    
    plot(dh, yaxt='n', xaxt='n', ann=FALSE, frame.plot=FALSE, xlim=c(round(min.val-10, digits=-1),round(max.val+10, digits=-1)), ylim=c(0,max(y.vec)), col="blue")
    polygon(dh, col="white", border="black", lwd=2.5)
    
    
    mn.hp <- mean(pheno.trt.fam[,3], na.rm=T)
    abline(v=Lata.val, col="blue", lwd=3)
    abline(v=donor.val, col="forestgreen", lwd=3)
    #text(x=round(min.val, digits=-1), y=(max(y.vec)/2), label=families[j], cex=2.5)
    #points(x=donor.val, y=(max(y.vec)/2), col="blue", pch=16, cex=2)
    points(x=mn.hp, y=(max(y.vec)), col="red", bg="red", pch=25, cex=3)
    
  }
  #par(mar=c(1,0,1,0))
  #plot(1, xlim=c(10,200))
  #axis(side=1, 10:200, labels = , col.axis = "black")
  #axis(side=1, at=seq(round(min.val-10, digits=-1),round(max.val+10, digits=-1),10), labels=seq(round(min.val-10, digits=-1),round(max.val+10, digits=-1),10), tick=F, cex.axis=2)
  axis(side=1, at=seq(axis.val[i,1], axis.val[i,2], axis.val[i,3]), labels=seq(axis.val[i,1], axis.val[i,2], axis.val[i,3]), tick=F, cex.axis=2.5)
  
  dev.off()
}
