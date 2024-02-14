setwd("C://Users/nm18d357/Dropbox/Boulot/CMPG_EvolEdge/SimuAdmixture/")
library(scales)
library(ggplot2)
library(reshape2)
library(viridis)


################################################
##General Info
Name="SimAdm"
deffile=read.table("SimAdm.def",header=TRUE)
NDef=dim(deffile)[1]
listSS=unique(deffile[,"SS1"])
listNS=unique(deffile[,"NS1"])
listTDiv=unique(deffile[,"TDiv"])
listTAdm=unique(deffile[,"TAdm"])
LevelsTDivTAdm=c("500_1","500_10","500_100","500_300","1000_1","1000_10","1000_100","1000_300")
LevelsSSNS=c("4_2000","4_5000","4_10000","8_2000","8_5000","8_10000","16_2000","16_5000","16_10000")
NRun=1
Length=1e8

#Different number of samples in the source populations depending on the def (2, 4 or 8 diploids)
for (nS in c(2,4,8)){
  assign(paste("Labels_nS",nS,sep=""),c(rep("Source1",nS),rep("Source2",nS),rep("Adm0%",10),rep("Adm5%",10),rep("Adm10%",10),rep("Adm20%",10),rep("Adm30%",10)))
  assign(paste("Col_nS",nS,sep=""),c(rep("darkblue",nS),rep("chartreuse2",nS),rep("gold",10),rep("darkorange1",10),rep("red",10),rep("darkred",10),rep("black",10)))
  assign(paste("Pch_nS",nS,sep=""),c(rep(8,2*nS),rep(20,50)))
  assign(paste("Cex_nS",nS,sep=""),c(rep(1.5,2*nS),rep(1,50)))
}

################################################
##Compute Pxy distance between samples##
# ##!!DONE ONCE NO NEED TO REDO!!
# MyFunction=function(a,b){
#   return((a/2+b/2-2*(a/2)*(b/2)))} #dxy
# MySum=function(Mat,Indi){
#   print(dim(Mat))
#   return(sweep(Mat,MARGIN=1,STATS=Mat[,Indi],FUN=MyFunction))}
# 
# for (def in c(1:NDef)){
#   print(def)
#   if (def %in% c(1:24)){
#     Labels=Labels_nS2
#   } else if (def %in% c(25:48)){
#     Labels=Labels_nS4
#   } else {Labels=Labels_nS8}
# 
#   Gen=read.table(paste(Name,"/",Name,"_",def,"_1.gen",sep=""), header = TRUE)
#   Gen=Gen[,5:ncol(Gen)]
#   NSamples=ncol(Gen)
#   #Convert the hapliod genotypes into diploid
#   Diplo=Gen[,seq(1,NSamples,by=2)]+Gen[,seq(2,NSamples,by=2)]
# 
#   Mat_Pxy=matrix(nrow=ncol(Diplo), ncol=ncol(Diplo), dimnames=list(colnames(Diplo),colnames(Diplo)), NA)
#   for (i in c(1:dim(Mat_Pxy)[1])){
#     print(i)
#     Mat_Pxy[i,]=colSums(MySum(Diplo,i),na.rm=TRUE) #Sum of dxy
#     Mat_Pxy[i,i]=0
#   }
#   Mat_Pxy=Mat_Pxy/Length
#   colnames(Mat_Pxy)=Labels
#   rownames(Mat_Pxy)=Labels
#   write.table(Mat_Pxy, file = paste("Pxy/Pxy_",Name,"_def",def,".txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"))
# }

################################################
##Plot MDS from the Pxy distances##
#I plot on the same page the 4 TAdm for the same SS, NS and TDiv 
for (SS in listSS){
  for (NS in listNS){
    for (TDiv in listTDiv){
      W=which((deffile[,"SS1"]==SS)&(deffile[,"NS1"]==NS)&(deffile[,"TDiv"]==TDiv))
      pdf(file=paste("MDS/Pxy_",Name,"_SS",SS,"_NS",NS,"_NT2000_TDiv",TDiv,".pdf",sep=""))
      par(mfrow = c(2, 2))
      for (def in W){
        print(def)
        if (def %in% c(1:24)){
          Col=Col_nS2; Pch=Pch_nS2; Cex=Cex_nS2
        } else if (def %in% c(25:48)){
          Col=Col_nS4; Pch=Pch_nS4; Cex=Cex_nS4
        } else {Col=Col_nS8; Pch=Pch_nS8; Cex=Cex_nS8}
        
        Mat_Pxy=read.table(paste("Pxy/Pxy_",Name,"_def",def,".txt",sep=""), header=FALSE, skip = 1)
        Labels=Mat_Pxy[,1]
        Mat_Pxy=Mat_Pxy[,-1]
        data.mds<-cmdscale(Mat_Pxy, eig=TRUE)
        plot(data.mds$points[,1:2],col=as.character(Col),pch=as.numeric(Pch),cex=as.numeric(Cex),cex.axis=0.7,cex.lab=0.9,
             main=NA,
             xlab=paste("Coordinate 1 (",round(abs(data.mds$eig[1])*100/sum(abs(data.mds$eig)),2),"%)",sep=""),
             ylab=paste("Coordinate 2 (",round(abs(data.mds$eig[2])*100/sum(abs(data.mds$eig)),2),"%)",sep=""))
        title(paste("TAdm=",deffile[def,"TAdm"],"g",sep=""), line = 0.5)
        if (def==W[4]){legend("topleft",legend=as.character(Labels[!duplicated(Labels)]),ncol=1,text.col=as.character(Col[!duplicated(Labels)]),pch=NA,cex=0.7, inset=0.01,adj=0.1, box.col="lightgrey")}
      }
      dev.off()
    }
  }
}


################################################
##FST between sources from the Pxy distances##
#(Hudson et al. 1992 Genetics, eq. 3)
#FST= 1 - Hw/Hb =  (Hb-Hw)/Hb
#where Hw =(Pxy_intraPop1 + Pxy_intraPop2)/2 and Hb is Pxy_interPop1-2
Pop=c("Source1","Source2")
Fst=matrix(ncol=1,nrow=NDef,dimnames = list(c(1:NDef),"Source1vs2"))
for (def in c(1:NDef)){
  print(def)
  Mat_Pxy=read.table(paste("Pxy/Pxy_",Name,"_def",def,".txt",sep=""), header=FALSE, skip = 1)
  Labels=Mat_Pxy[,1]
  Mat_Pxy=Mat_Pxy[,-1]
  
  Mean_Pxy=matrix(ncol=length(Pop), nrow=length(Pop), dimnames=list(Pop,Pop))
  for (i in Pop){
    for (j in Pop){
      Sub=Mat_Pxy[which(Labels==i),which(Labels==j)]
      if (i==j){
        Sub=Sub[lower.tri(Sub, diag = FALSE)]
        Mean_Pxy[i,j]=mean(Sub)
      } else {Mean_Pxy[i,j]=mean(unlist(Sub))}
    }
  }
  Fst[def,1]=(100*(1-((Mean_Pxy[Pop[1],Pop[1]]+Mean_Pxy[Pop[2],Pop[2]])/2)/Mean_Pxy[Pop[1],Pop[2]]))
}
write.table(Fst, file = "Fst/FstbtwSources.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

#PLOT
MyFst=as.data.frame(cbind(paste(deffile[,"TDiv"],deffile[,"TAdm"],sep="_"),paste(deffile[,"SS1"],deffile[,"NS1"],sep="_"),round(Fst,0)))
colnames(MyFst)=c("TDiv_TAdm","SS_NS","Fst")
MyFst$TDiv_TAdm<-factor(MyFst$TDiv_TAdm, levels=LevelsTDivTAdm)
MyFst$SS_NS<-factor(MyFst$SS_NS, levels=LevelsSSNS)
MyFst$Fst<-as.numeric(as.character(MyFst$Fst))

pdf("Fst/FstbtwSources.pdf")
ggplot(data = MyFst, aes(x=TDiv_TAdm, y=SS_NS, fill=Fst)) + 
  geom_tile() +
  scale_fill_viridis(option="plasma",limits = c(min(round(Fst,0),na.rm=TRUE), max(round(Fst,0),na.rm=TRUE)), oob = scales::squish,na.value="white") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=3.5, color="white",lty=1,lwd=2) +
  geom_hline(yintercept=6.5, color="white",lty=1,lwd=2) +
  geom_hline(yintercept=6.5, color="white",lty=1,lwd=2) +
  geom_vline(xintercept=4.5, color="white",lty=1,lwd=2)
dev.off()


