argsFull <- commandArgs()
Rscript <- argsFull[1]
scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
RepositoryPath = sub("/scripts","",scriptPath)
setwd(RepositoryPath)
dataToStudy = read.table(paste(RepositoryPath,"/data/RNAseqBP.txt",sep=""),header=T,sep="\t",dec=".")

dataToStudy$classBranch = "No"
dataToStudy$classBranch[dataToStudy$branchpoint_prob>0] = "Yes"

table(dataToStudy$classBranch, dataToStudy$UsedSite)
######################################
# get ROC analysis
######################################

library(ROCR)

labels = dataToStudy$UsedSite
BPP = dataToStudy$bp_zsc
SVM = dataToStudy$svm_scr
LB = dataToStudy$LB_score
RNABPS = dataToStudy$RNABPS_score

score=data.frame(BPP,labels,SVM,labels,LB,labels,RNABPS,labels)
pred=prediction(predictions=score[,c(1,3,5,7)],labels=score[,c(2,4,6,8)])
perf=performance(pred,"tpr","fpr")
auc<- performance(pred,"auc")
auc<- unlist(slot(auc, "y.values"))
auc=round(auc,3)

rang = order(auc,decreasing = TRUE)
annot = c("BPP, AUC =","SVMBP-finder, AUC =","LaBranchoR, AUC =","RNABPS, AUC =")
color =  c("blue","red","black","orange")
lty = c(2, 1, 3, 4)

pdf(file="RocScoreRNAseq.pdf",width =10, height = 10)
plot(perf,col= as.list(color),cex.lab=1.5,lty=as.list(lty), lwd=4,cex.main=2,xlab="1-Specificity",
	ylab="Sensitivity",cex.axis=5,text.adj=c(1.2,0.1),xlim=c(0,1),ylim=c(0,1),colorize=F)
legend(x=0.45, y=0.3, paste(annot[rang],auc[rang]), col = color[rang],cex=1.5,text.col = "black", lty = lty[rang],lwd=3)
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

rocTableMAPP_BPP = getROCtable(pred,1)
rocTablePON_SVM = getROCtable(pred,2)
rocTableMAPP_LB = getROCtable(pred,3)
rocTableMAPP_RNABPS = getROCtable(pred,4)

write.table(rocTableMAPP_BPP, "rocTable_RNAseq_BPP.txt", sep="\t",dec=",",row.names=F)
write.table(rocTablePON_SVM, "rocTable_RNAseq_SVM.txt", sep="\t",dec=",",row.names=F)
write.table(rocTableMAPP_LB, "rocTable_RNAseq_LaBranchoR.txt", sep="\t",dec=",",row.names=F)
write.table(rocTableMAPP_RNABPS, "rocTable_RNAseq_RNABPS.txt", sep="\t",dec=",",row.names=F)

######################################
# expression and presence or not of Branch point
######################################

dataToTest<-dataToStudy[!is.na(dataToStudy$Expression_Average_...),]
dataToTest<-dataToTest[dataToTest$Expression_Average_...>0,]

meanBranchYes <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$branchpoint_class=="Yes" & dataToTest$Expression_Average_... < 8000])),2)
N_BranchYes <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$branchpoint_class=="Yes"]))
meanBranchNo <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$branchpoint_class=="No" & dataToTest$Expression_Average_... < 8000])),2)
N_BranchNo <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$branchpoint_class=="No"]))

meanBPPYes <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$BPP_class=="Yes" & dataToTest$Expression_Average_... < 8000])),2)
N_BPPYes <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$BPP_class=="Yes"]))
meanBPPNo <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$BPP_class=="No" & dataToTest$Expression_Average_... < 8000])),2)
N_BPPNo <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$BPP_class=="No"]))

meanSVMYes <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$SVM_class=="Yes" & dataToTest$Expression_Average_... < 8000])),2)
N_SVMYes <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$SVM_class=="Yes"]))
meanSVMNo <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$SVM_class=="No" & dataToTest$Expression_Average_... < 8000])),2)
N_SVMNo <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$SVM_class=="No"]))

meanLBYes <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$LB_class=="Yes" & dataToTest$Expression_Average_... < 8000])),2)
N_LBYes <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$LB_class=="Yes"]))
meanLBNo <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$LB_class=="No" & dataToTest$Expression_Average_... < 8000])),2)
N_LBNo <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$LB_class=="No"]))

meanRNABPSYes <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$RNABPS_class=="Yes" & dataToTest$Expression_Average_... < 8000])),2)
N_RNABPSYes <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$RNABPS_class=="Yes"]))
meanRNABPSNo <- round(mean(na.omit(dataToTest$Expression_Average_...[dataToTest$RNABPS_class=="No" & dataToTest$Expression_Average_... < 8000])),2)
N_RNABPSNo <- length(na.omit(dataToTest$Expression_Average_...[dataToTest$RNABPS_class=="No"]))

bornMax = 10000
pdf(file="BoxplotPresenceRNAseq.pdf",width =20, height = 6)
par(yaxt="n",xaxt = "s", mfrow=c(1,5), mai = c(1, 1, 0, 0),oma = c(0.5, 0, 2, 0), font.lab=2, font.axis=2,bty="n",cex.lab=2,cex.axis=1.5, mex = 2)
boxplot(Expression_Average_...~SVM_class, data = dataToTest[dataToTest$UsedSite==1,],log="y",names=c(paste("No, N = ",N_SVMNo,sep=""),paste("Yes, N = ",N_SVMYes,sep="")),
	ylog=T, xlab = paste("SVM-BPfinder\n(No: ",meanSVMNo," %; Yes: ",meanSVMYes," %)",sep=""),yaxt="s",ylab = "Expression (log scale)")
segments(x0 = 1, x1 = 2, y0 = bornMax+1000, y1 = bornMax+1000)
text(x= 1.5,y = bornMax+8000,"***",cex=1.5)
par(mai = c(1, 0, 0, 0))
boxplot(Expression_Average_...~BPP_class, data = dataToTest[dataToTest$UsedSite==1,],log="y",names=c(paste("No, N = ",N_BPPNo,sep=""),paste("Yes, N = ",N_BPPYes,sep="")),
	ylog=T, xlab = paste("BPP\n(No: ",meanBPPNo," %; Yes: ",meanBPPYes," %)",sep=""))
segments(x0 = 1, x1 = 2, y0 = bornMax+1000, y1 = bornMax+1000)
text(x= 1.5,y = bornMax+8000,"***",cex=1.5)
boxplot(Expression_Average_...~branchpoint_class, data = dataToTest[dataToTest$UsedSite==1,],log="y",names=c(paste("No, N = ",N_BranchNo,sep=""),paste("Yes, N = ",N_BranchYes,sep="")),
	ylog=T, xlab =  paste("Branchpointer\n(No: ",meanBranchNo," %; Yes: ",meanBranchYes," %)",sep=""))
segments(x0 = 1, x1 = 2, y0 = bornMax+1000, y1 = bornMax+1000)
text(x= 1.5,y = bornMax+8000,"***",cex=1.5)
boxplot(Expression_Average_...~LB_class, data = dataToTest[dataToTest$UsedSite==1,],log="y",names=c(paste("No, N = ",N_LBNo,sep=""),paste("Yes, N = ",N_LBYes,sep="")),
	ylog=T, xlab =  paste("LaBranchoR\n(No: ",meanLBNo," %; Yes: ",meanLBYes," %)",sep=""))
segments(x0 = 1, x1 = 2, y0 = bornMax+1000, y1 = bornMax+1000)
text(x= 1.5,y = bornMax+8000,"***",cex=1.5)
boxplot(Expression_Average_...~RNABPS_class, data = dataToTest[dataToTest$UsedSite==1,],log="y",names=c(paste("No, N = ",N_RNABPSNo,sep=""),paste("Yes, N = ",N_RNABPSYes,sep="")),
	ylog=T, xlab =  paste("RNABPS\n(No: ",meanRNABPSNo," %; Yes: ",meanRNABPSYes," %)",sep=""))
segments(x0 = 1, x1 = 2, y0 = bornMax+1000, y1 = bornMax+1000)
text(x= 1.5,y = bornMax+8000,"***",cex=1.5)
dev.off()

######################################
# expression and branch point score
######################################

#correlation expression
library(ggplot2)
dataExpr = dataToStudy[dataToStudy$UsedSite==1 & !is.na(dataToStudy$Expression_Average_...),]
dataExpr = dataExpr[dataExpr$Expression_Average_...!=100,]
dataExpr = dataExpr[dataExpr$Expression_Average_...!=0,]

pdf(file="CorrScoreExprRNAseq.pdf",width =6, height = 5)
par (mfrow=c(1,5))
l = lm(log(Expression_Average_...)~svm_scr, data = dataExpr[!is.na(dataExpr$svm_scr),])
ggplot(dataExpr[!is.na(dataExpr$svm_scr),]) +
	aes(x = svm_scr, y = Expression_Average_..., fill = ..level..) +
	stat_density_2d(geom = "polygon",n=200) +
	xlab(paste("SVM-BPfinder, R² =", round(summary(l)$r.squared,4),
		", p-value = ",cor.test(log(dataExpr$Expression_Average_...),dataExpr$svm_scr)$p.value,sep="")) +
	ylab("Expression (log scale)") +
	scale_y_log10() +
	theme_classic() +
	geom_abline(intercept = coef(l)[1], slope = coef(l)[2],color="red") +
	labs(fill = "Densité")
l = lm(log(Expression_Average_...)~bp_zsc, data = dataExpr[!is.na(dataExpr$bp_zsc),])
ggplot(dataExpr[!is.na(dataExpr$bp_zsc),]) +
	aes(x = bp_zsc, y = Expression_Average_..., fill = ..level..) +
	stat_density_2d(geom = "polygon",n=200) +
	xlab(paste("BPP, R² =", round(summary(l)$r.squared,4),
		", p-value = ",cor.test(log(dataExpr$Expression_Average_...),dataExpr$bp_zsc)$p.value,sep="")) +
	ylab("Expression (log scale)") +
	scale_y_log10() +
	theme_classic() +
	geom_abline(intercept = coef(l)[1], slope = coef(l)[2],color="red") +
	labs(fill = "Densité")
l = lm(log(Expression_Average_...)~branchpoint_prob, data = dataExpr[dataExpr$branchpoint_prob!=0,])
ggplot(dataExpr[dataExpr$branchpoint_prob!=0,]) +
	aes(x = branchpoint_prob, y = Expression_Average_..., fill = ..level..) +
	stat_density_2d(geom = "polygon",n=200) +
	xlab(paste("Branchpointer, R² =", round(summary(l)$r.squared,4),
		", p-value = ",cor.test(log(dataExpr$Expression_Average_...[dataExpr$branchpoint_prob!=0]),
		dataExpr$branchpoint_prob[dataExpr$branchpoint_prob!=0])$p.value,sep="")) +
	ylab("Expression (log scale)") +
	scale_y_log10() +
	theme_classic() +
	geom_abline(intercept = coef(l)[1], slope = coef(l)[2],color="red") +
	labs(fill = "Densité")
l = lm(log(Expression_Average_...)~LB_score, data = dataExpr[!is.na(dataExpr$LB_score),])
ggplot(dataExpr[!is.na(dataExpr$LB_score),]) +
	aes(x = LB_score, y = Expression_Average_..., fill = ..level..) +
	stat_density_2d(geom = "polygon",n=200) +
	xlab(paste("LaBranchoR, R² =", round(summary(l)$r.squared,4),
		", p-value = ",cor.test(log(dataExpr$Expression_Average_...),dataExpr$LB_score)$p.value,sep="")) +
	ylab("Expression (log scale)") +
	scale_y_log10() +
	theme_classic() +
	geom_abline(intercept = coef(l)[1], slope = coef(l)[2],color="red") +
	labs(fill = "Densité")
l = lm(log(Expression_Average_...)~RNABPS_score, data = dataExpr[!is.na(dataExpr$RNABPS_score),])
ggplot(dataExpr[!is.na(dataExpr$RNABPS_score),]) +
	aes(x = RNABPS_score, y = Expression_Average_..., fill = ..level..) +
	stat_density_2d(geom = "polygon",n=200) +
	xlab(paste("RNABPS, R² =", round(summary(l)$r.squared,4),
		", p-value = ",cor.test(log(dataExpr$Expression_Average_...),dataExpr$RNABPS_score)$p.value,sep="")) +
	ylab("Expression (log scale)") +
	scale_y_log10() +
	theme_classic() +
	geom_abline(intercept = coef(l)[1], slope = coef(l)[2],color="red") +
	labs(fill = "Densité")
dev.off()

######################################
# get VENN diagramm
######################################

#upload library
library(VennDiagram)
library(gplots)

dataToStudy$SVM_class[is.na(dataToStudy$SVM_class)] <- "No"
dataToStudy$BPP_class[is.na(dataToStudy$BPP_class)] = "No"
dataToStudy$branchpoint_class[is.na(dataToStudy$branchpoint_class)] = "No"
dataToStudy$LB_class[is.na(dataToStudy$LB_class)] = "No"
dataToStudy$RNABPS_class[is.na(dataToStudy$RNABPS_class)] = "No"

dataVenn = dataToStudy[!(dataToStudy$SVM_class=="No" & dataToStudy$BPP_class=="No" & dataToStudy$branchpoint_class=="No" &
				dataToStudy$LB_class=="No" & dataToStudy$RNABPS_class=="No"),]

dataVenn$id=c(1:nrow(dataVenn))
SVM = dataVenn[dataVenn$SVM_class=="Yes","id"]
BPP = dataVenn[dataVenn$BPP_class=="Yes","id"]
Branch = dataVenn[dataVenn$branchpoint_class=="Yes","id"]
LB = dataVenn[dataVenn$LB_class=="Yes","id"]
RNABPS = dataVenn[dataVenn$RNABPS_class=="Yes","id"]

SVM_TP = dataVenn[dataVenn$SVM_class=="Yes" & dataVenn$UsedSite==1,"id"]
BPP_TP = dataVenn[dataVenn$BPP_class=="Yes" & dataVenn$UsedSite==1,"id"]
Branch_TP = dataVenn[dataVenn$branchpoint_class=="Yes" & dataVenn$UsedSite==1,"id"]
LB_TP = dataVenn[dataVenn$LB_class=="Yes" & dataVenn$UsedSite==1,"id"]
RNABPS_TP = dataVenn[dataVenn$RNABPS_class=="Yes" & dataVenn$UsedSite==1,"id"]

SVM_FP = dataVenn[dataVenn$SVM_class=="Yes" & dataVenn$UsedSite==0,"id"]
BPP_FP = dataVenn[dataVenn$BPP_class=="Yes" & dataVenn$UsedSite==0,"id"]
Branch_FP = dataVenn[dataVenn$branchpoint_class=="Yes" & dataVenn$UsedSite==0,"id"]
LB_FP = dataVenn[dataVenn$LB_class=="Yes" & dataVenn$UsedSite==0,"id"]
RNABPS_FP = dataVenn[dataVenn$RNABPS_class=="Yes" & dataVenn$UsedSite==0,"id"]

tmp = venn(list(SVM,BPP,Branch,LB,RNABPS),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),show.plot=F)
summary(attr(tmp,'intersections'))
venn.diagram(
	x = list(SVM , BPP, Branch, LB, RNABPS),category.names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),
	filename = 'RNAseq_DiagVenn.tiff',
	height = 6500, width = 6500, resolution = 600, imagetype = "tiff", units = "px",
	print.mode = "raw",
	col = "black",
	fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	alpha = 0.50,
	cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
	 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
	cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	cat.cex = 1.5,
	cat.fontface = "bold",
	margin = 0.05
)

tmp = venn(list(SVM_TP,BPP_TP,Branch_TP,LB_TP,RNABPS_TP),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),show.plot=F)
summary(attr(tmp,'intersections'))
venn.diagram(
	x = list(SVM_TP , BPP_TP, Branch_TP, LB_TP, RNABPS_TP),category.names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),
	filename = 'RNAseq_DiagVenn_TP.tiff',
	height = 6500, width = 6500, resolution = 600, imagetype = "tiff", units = "px",
	print.mode = "raw",
	col = "black",
	fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	alpha = 0.50,
	cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
	 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
	cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	cat.cex = 1.5,
	cat.fontface = "bold",
	margin = 0.05
)

tmp = venn(list(SVM_FP,BPP_FP,Branch_FP,LB_FP,RNABPS_FP),names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),show.plot=F)
summary(attr(tmp,'intersections'))
venn.diagram(
	x = list(SVM_FP , BPP_FP, Branch_FP, LB_FP, RNABPS_FP),category.names = c("SVM-BPfinder", "BPP", "Branchpointer", "LaBranchoR", "RNABPS"),
	filename = 'RNAseq_DiagVenn_FP.tiff',
	height = 6500, width = 6500, resolution = 600, imagetype = "tiff", units = "px",
	print.mode = "raw",
	col = "black",
	fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	alpha = 0.50,
	cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
	 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
	cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
	cat.cex = 1.5,
	cat.fontface = "bold",
	margin = 0.05
)
