argsFull <- commandArgs()
Rscript <- argsFull[1]
scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
RepositoryPath = sub("/scripts","",scriptPath)
setwd(RepositoryPath)
message("Import data...")
load(paste(RepositoryPath,"/data/ensemblBP.RData",sep=""))

library(ROCR)
library(gplots)
message("Make ID...")
ID = c(1:nrow(poolData))

######################################
# get ROC analysis
######################################

message("Make table for ROC analysis...")
labels = poolData$UsedSite
score=data.frame(MES = poolData$MES,labels,
    SSF = poolData$SSF,labels,
    SVM = poolData$zscSVM,labels,
    BPP = poolData$zscBPP,labels,
    LB = poolData$scoreLB,labels,
    RNABPS = poolData$scoreRNABPS,labels
)

message("Make Predrediction...")
pred=prediction(predictions=score[,c(5,7,9,11)],labels=score[,c(6,8,10,12)])
perf=performance(pred,"tpr","fpr")
auc<- performance(pred,"auc")
auc<- unlist(slot(auc, "y.values"))
auc<- round(auc,3)

rang = order(auc,decreasing = TRUE)
annot = c("SVM-BPfinder, AUC =","BPP, AUC =","LaBranchoR, AUC =","RNABPS, AUC =")
color =  c("red","blue","black","orange")
lty = c(1:4)

message("Save PDF file...")
pdf(file="RocMergeBP.pdf",width =10, height = 10)
plot(perf, col= as.list(color),lty= as.list(lty),cex.lab=1.5, lwd=4,cex.main=2,xlab="1-Specificity",
	ylab="Sensitivity",cex.axis=5,text.adj=c(1.2,0.1),xlim=c(0,1),ylim=c(0,1),colorize=F,type="l")
legend(x=0.45, y=0.3, paste(annot[rang],auc[rang]), col = color[rang], lty = lty[rang],cex=1.5,text.col = "black",lwd=3)
segments(x0=0,y0=0,x1=1,y1=1,lty=2,col='grey')
dev.off()

getROCtable <- function(pred,i){
	cutoffs = unlist(slot(pred, "cutoffs")[i])
	tp = unlist(slot(pred, "tp")[i])
	fp = unlist(slot(pred, "fp")[i])
	tn = unlist(slot(pred, "tn")[i])
	fn = unlist(slot(pred, "fn")[i])
	Se =tp / (tp + fn)
	Sp =tn / (tn + fp)
	diff_Se_Sp = abs(Se - Sp)
	ACC = (tp+tn)/(tp+tn+fn+fp)
	VPP = tp / (tp + fp)
	VPN = tn / (tn + fn)
	data_ROC = data.frame(cutoffs, tp, fp, tn, fn, Se, Sp, diff_Se_Sp, ACC, VPP, VPN)
	return(data_ROC)
}

message("Make ROC table Output...")
rocTable_SVM = getROCtable(pred,1)
rocTable_BPP = getROCtable(pred,2)
rocTable_LB = getROCtable(pred,3)
rocTable_RNABPS = getROCtable(pred,4)
message("Save ROC table Output...")
write.table(rocTable_SVM, "rocTableSVM.txt", sep="\t",dec=",",row.names=F)
write.table(rocTable_BPP, "rocTableBPP.txt", sep="\t",dec=",",row.names=F)
write.table(rocTable_LB, "rocTableLB.txt", sep="\t",dec=",",row.names=F)
write.table(rocTable_RNABPS, "rocTableRNABPS.txt", sep="\t",dec=",",row.names=F)

message('score: cutoff; TP; FP; TN; FN; Se; Sp; diff; VPP; VPN')
message(paste("SVM-BPfinder: ",paste(rocTable_SVM[rocTable_SVM$diff_Se_Sp == min(rocTable_SVM$diff_Se_Sp),],collapse="; ")))
message(paste("BPP: ",paste(rocTable_BPP[rocTable_BPP$diff_Se_Sp == min(rocTable_BPP$diff_Se_Sp),],collapse="; ")))
message(paste("LaBranchoR: ",paste(rocTable_LB[rocTable_LB$diff_Se_Sp == min(rocTable_LB$diff_Se_Sp),],collapse="; ")))
message(paste("RNABPS: ",paste(rocTable_RNABPS[rocTable_RNABPS$diff_Se_Sp == min(rocTable_RNABPS$diff_Se_Sp),],collapse="; ")))

######################################
# get VENN diagramm
######################################

#Cutoffs from ROC analysis
cutOff_SVM = 0.706
cutOff_BPP = 5.384
cutOff_LB = 0.653
cutOff_RNABPS = 0.653

message("Remove NAs...")
poolData$zscSVM[is.na(poolData$zscSVM)] = 0
poolData$zscBPP[is.na(poolData$zscBPP)] = 0
poolData$scoreLB[is.na(poolData$scoreLB)] = 0
poolData$scoreRNABPS[is.na(poolData$scoreRNABPS)] = 0

message("Make vector...")
SVM = ID[poolData$zscSVM>=cutOff_SVM]
BPP = ID[poolData$zscBPP>=cutOff_BPP]
Branch = ID[poolData$Branch_class=="Yes"]
LB = ID[poolData$scoreLB>=cutOff_LB]
RNABPS = ID[poolData$scoreRNABPS>=cutOff_RNABPS]

message("Print first PDF...")
pdf(file="VennMergeAll.pdf",width =10, height = 10)
tmp = venn(list(SVM,BPP,Branch,LB,RNABPS),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"))
print(summary(attr(tmp,"intersections")))
dev.off()

message("Print second PDF...")
SVM_TP = ID[poolData$zscSVM>=cutOff_SVM & poolData$UsedSite==1]
BPP_TP = ID[poolData$zscBPP>=cutOff_BPP & poolData$UsedSite==1]
Branch_TP = ID[poolData$Branch_class=="Yes" & poolData$UsedSite==1]
LB_TP = ID[poolData$scoreLB>=cutOff_LB & poolData$UsedSite==1]
RNABPS_TP = ID[poolData$scoreRNABPS>=cutOff_RNABPS & poolData$UsedSite==1]

pdf(file="VennMergeTP.pdf",width =10, height = 10)
tmp = venn(list(SVM_TP,BPP_TP,Branch_TP,LB_TP,RNABPS_TP),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"))
print(summary(attr(tmp,"intersections")))
dev.off()

message("Print third PDF...")
SVM_FP = ID[poolData$zscSVM>=cutOff_SVM & poolData$UsedSite==0]
BPP_FP = ID[poolData$zscBPP>=cutOff_BPP & poolData$UsedSite==0]
Branch_FP = ID[poolData$Branch_class=="Yes" & poolData$UsedSite==0]
LB_FP = ID[poolData$scoreLB>=cutOff_LB & poolData$UsedSite==0]
RNABPS_FP = ID[poolData$scoreRNABPS>=cutOff_RNABPS & poolData$UsedSite==0]

pdf(file="VennMergeFP.pdf",width =10, height = 10)
tmp = venn(list(SVM_FP,BPP_FP,Branch_FP,LB_FP,RNABPS_FP),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"))
print(summary(attr(tmp,"intersections")))
dev.off()
