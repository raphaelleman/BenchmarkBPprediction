argsFull <- commandArgs()
Rscript <- argsFull[1]
scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
RepositoryPath = sub("/scripts","",scriptPath)
setwd(RepositoryPath)
dataResult = read.table(paste(RepositoryPath,"/data/variantBP.txt",sep=""), sep="\t", header=T, dec=",")

####################################
#repartition variant
####################################

dataBP = dataResult[,c("class_effect","distSS")]
Obs = dataBP$class_effect
distSS = dataBP$distSS
distSScomp = c(-44:-18)[-which(c(-44:-18)%in%distSS)]
distSS = c(distSScomp,distSS)
Obs = c(rep("NA",length(distSScomp)), Obs)
effVar=table(Obs,distSS)

pdf("ReparVar.pdf", width = 16, height = 9)

barplot(effVar[1:2,],ylim=c(0,16),col=c('gray55','black'),xlab="",ylab="",cex.main=2.5, main="Branchpoint area",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=15,x=0.2,labels="Number of variants for each position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(14,length(effVar[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar[1,])-1)),by=1.2),
	labels=paste("-",c(44:18),sep=""),col='black',cex=1.5)
text(y=rep(12,length(effVar[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar[1,])-1)),by=1.2),
	labels=as.numeric(effVar[1,]),col='gray55',cex=2)
text(y=rep(13,length(effVar[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar[2,])-1)),by=1.2),
	labels=as.numeric(effVar[2,]),col='black',cex=2)

dev.off()

######################################
# get ROC analysis
######################################

library(ROCR)

labels = dataResult$class_effect
Branch = dataResult$Delta_prob*(-1)
BPP = dataResult$Delta_BPP*(-1)
SVM = dataResult$Delta_SVM*(-1)
HSF = dataResult$delta_HSF*(-1)
LB = dataResult$Delta_LB*(-1)
RBPS = dataResult$Delta_RBPS*(-1)

score=data.frame(Branch,labels,BPP,labels,SVM,labels,HSF,labels,LB,labels,RBPS,labels)
pred_multi=prediction(predictions=score[,c(1,3,5,7,9,11)],labels=score[,c(2,4,6,8,10,12)])
perf_multi=performance(pred_multi,"tpr","fpr")
auc_multi<- performance(pred_multi,"auc")
auc_multi<- unlist(slot(auc_multi, "y.values"))
auc_multi=round(auc_multi,3)

rang = order(auc_multi,decreasing = TRUE)
annot = c("Branchpointer, AUC =","BPP, AUC =","SVMBP-finder, AUC =","HSF, AUC =","LaBranchoR, AUC =","RNABPS, AUC =")
color =  c("lightblue","blue","red","violet","black","orange")
lty = c(5, 2, 1, 5, 3, 4)

pdf(file="RocVariantBPscore.pdf",width =10, height = 10)
plot(perf_multi,col= as.list(color),cex.lab=1.5,lty = as.list(lty),lwd=5,cex.main=2,xlab="1-Specificity",
	ylab="Sensitivity",cex.axis=5,text.adj=c(1.2,0.1),xlim=c(0,1),ylim=c(0,1),colorize=F,type="l")
legend(x=0.45, y=0.3, paste(annot[rang],auc_multi[rang]),col =color[rang],lty =lty[rang],cex=1.5,text.col = "black",lwd=3)
segments(x0=0,y0=0,x1=1,y1=1,lty=2,col='grey')
dev.off()

getROCtable <- function(pred_multi,i){
	cutoffs = unlist(slot(pred_multi, "cutoffs")[i])
	tp = unlist(slot(pred_multi, "tp")[i])
	fp = unlist(slot(pred_multi, "fp")[i])
	tn = unlist(slot(pred_multi, "tn")[i])
	fn = unlist(slot(pred_multi, "fn")[i])
	Se =tp / (tp + fn)
	Sp =tn / (tn + fp)
	diff_Se_Sp = abs(Se - Sp)
	ACC = (tp+tn)/(tp+tn+fn+fp)
	VPP = tp / (tp + fp)
	VPN = tn / (tn + fn)
	data_ROC = data.frame(cutoffs, tp, fp, tn, fn, Se, Sp, diff_Se_Sp, ACC, VPP, VPN)
	return(data_ROC)
}

data_ROC_Branch = getROCtable(pred_multi,1)
data_ROC_BPP = getROCtable(pred_multi,2)
data_ROC_SVM = getROCtable(pred_multi,3)
data_ROC_HSF = getROCtable(pred_multi,4)
data_ROC_LB = getROCtable(pred_multi,5)
data_ROC_RBPS = getROCtable(pred_multi,6)

write.table(data_ROC_SVM,"data_ROC_variant_SVM.txt",sep="\t",dec=",",row.names=F,quote=F)
write.table(data_ROC_BPP,"data_ROC_variant_BPP.txt",sep="\t",dec=",",row.names=F,quote=F)
write.table(data_ROC_Branch,"data_ROC_variant_Branch.txt",sep="\t",dec=",",row.names=F,quote=F)
write.table(data_ROC_HSF,"data_ROC_variant_HSF.txt",sep="\t",dec=",",row.names=F,quote=F)
write.table(data_ROC_LB,"data_ROC_variant_LB.txt",sep="\t",dec=",",row.names=F,quote=F)
write.table(data_ROC_RBPS,"data_ROC_variant_RBPS.txt",sep="\t",dec=",",row.names=F,quote=F)

####################################
#repartition relatif variant
####################################

getPosBP <- function(intBP,SVM="NON"){
	intBPsplit = unlist(strsplit(as.character(intBP),"-",fixed = T))
	intBPsplitNum = as.numeric(intBPsplit[length(intBPsplit)])+1
	return (intBPsplitNum)
}

posBP_Branch = rep(0,nrow(dataResult))
posBP_BPP = rep(0,nrow(dataResult))
posBP_SVM = rep(0,nrow(dataResult))
posBP_LB = rep(0,nrow(dataResult))
posBP_RBPS = rep(0,nrow(dataResult))

for(i in 1:nrow(dataResult)){
	posBP_Branch[i] = getPosBP(dataResult[i,"ProbPBarea_Branch"])
	posBP_BPP[i] = getPosBP(dataResult[i,"ProbPBarea_BPP"])
	posBP_SVM[i] = getPosBP(dataResult[i,"ProbPBarea_SVM"])
	posBP_LB[i] = getPosBP(dataResult[i,"ProbPBarea_LB"])
	posBP_RBPS[i] = getPosBP(dataResult[i,"ProbPBarea_RBPS"])
}

#repartition variant
Obs = dataResult$class_effect
distSS = dataResult$distSS
distSS_Branch = distSS + abs(posBP_Branch)
distSS_BPP = distSS + abs(posBP_BPP)
distSS_SVM = distSS + abs(posBP_SVM)
distSS_LB = distSS + abs(posBP_LB)
distSS_RBPS = distSS + abs(posBP_RBPS)
distSS = c(distSS_Branch,distSS_BPP,distSS_SVM,distSS_LB,distSS_RBPS)
distSS = distSS[!duplicated(distSS)]
effVar_Branch = table(Obs,distSS_Branch)
compBranch = matrix(data = 0,nrow=2,ncol=length(as.character(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_Branch)))]))))
rownames(compBranch) <- c(0,1)
colnames(compBranch) <- c(as.numeric(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_Branch)))])))
effVar_Branch2 = cbind(effVar_Branch,compBranch)
effVar_Branch2 = effVar_Branch2[,order(as.numeric(colnames(effVar_Branch2)))]

effVar_BPP = table(Obs,distSS_BPP)
compBPP = matrix(data = 0,nrow=2,ncol=length(as.character(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_BPP)))]))))
rownames(compBPP) <- c(0,1)
colnames(compBPP) <- c(as.numeric(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_BPP)))])))
effVar_BPP2 = cbind(effVar_BPP,compBPP)
effVar_BPP2 = effVar_BPP2[,order(as.numeric(colnames(effVar_BPP2)))]

effVar_SVM = table(Obs,distSS_SVM)
compSVM = matrix(data = 0,nrow=2,ncol=length(as.character(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_SVM)))]))))
rownames(compSVM) <- c(0,1)
colnames(compSVM) <- c(as.numeric(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_SVM)))])))
effVar_SVM2 = cbind(effVar_SVM,compSVM)
effVar_SVM2 = effVar_SVM2[,order(as.numeric(colnames(effVar_SVM2)))]

effVar_LB = table(Obs,distSS_LB)
compLB = matrix(data = 0,nrow=2,ncol=length(as.character(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_LB)))]))))
rownames(compLB) <- c(0,1)
colnames(compLB) <- c(as.numeric(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_LB)))])))
effVar_LB2 = cbind(effVar_LB,compLB)
effVar_LB2 = effVar_LB2[,order(as.numeric(colnames(effVar_LB2)))]

effVar_RBPS = table(Obs,distSS_RBPS)
compRBPS = matrix(data = 0,nrow=2,ncol=length(as.character(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_RBPS)))]))))
rownames(compRBPS) <- c(0,1)
colnames(compRBPS) <- c(as.numeric(na.omit(distSS[-which(distSS%in%as.numeric(colnames(effVar_RBPS)))])))
effVar_RBPS2 = cbind(effVar_RBPS,compRBPS)
effVar_RBPS2 = effVar_RBPS2[,order(as.numeric(colnames(effVar_RBPS2)))]

pdf("ReparVarRelatif.pdf", width = 18, height = 25)

par(mfrow=c(5,1), mar = c(2,2,2,2))

barplot(effVar_Branch2[1:2,],ylim=c(0,24),col=c('gray55','black'),xlab="",ylab="",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=22,x=0.2,labels="Branchpointer, Number of variants for each relative position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(21,length(effVar_Branch2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_Branch2[1,])-1)),by=1.2),
	labels=as.numeric(colnames(effVar_Branch2)),col='black',cex=1.5)
text(y=rep(19,length(effVar_Branch2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_Branch2[1,])-1)),by=1.2),
	labels=as.numeric(effVar_Branch2[1,]),col='gray55',cex=2)
text(y=rep(20,length(effVar_Branch2[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_Branch2[2,])-1)),by=1.2),
	labels=as.numeric(effVar_Branch2[2,]),col='black',cex=2)

barplot(effVar_BPP2[1:2,],ylim=c(0,24),col=c('gray55','black'),xlab="",ylab="",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=22,x=0.2,labels="BPP, Number of variants for each relative position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(21,length(effVar_BPP2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_BPP2[1,])-1)),by=1.2),
	labels=as.numeric(colnames(effVar_BPP2)),col='black',cex=1.5)
text(y=rep(19,length(effVar_BPP2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_BPP2[1,])-1)),by=1.2),
	labels=as.numeric(effVar_BPP2[1,]),col='gray55',cex=2)
text(y=rep(20,length(effVar_BPP2[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_BPP2[2,])-1)),by=1.2),
	labels=as.numeric(effVar_BPP2[2,]),col='black',cex=2)

barplot(effVar_SVM2[1:2,],ylim=c(0,24),col=c('gray55','black'),xlab="",ylab="",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=22,x=0.2,labels="SVM-BPfinder, Number of variants for each relative position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(21,length(effVar_SVM2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_SVM2[1,])-1)),by=1.2),
	labels = as.numeric(colnames(effVar_SVM2)),col='black',cex=1.5)
text(y=rep(19,length(effVar_SVM2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_SVM2[1,])-1)),by=1.2),
	labels=as.numeric(effVar_SVM2[1,]),col='gray55',cex=2)
text(y=rep(20,length(effVar_SVM2[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_SVM2[2,])-1)),by=1.2),
	labels=as.numeric(effVar_SVM2[2,]),col='black',cex=2)

barplot(effVar_LB2[1:2,],ylim=c(0,24),col=c('gray55','black'),xlab="",ylab="",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=22,x=0.2,labels="LaBranchoR, Number of variants for each relative position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(21,length(effVar_LB2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_LB2[1,])-1)),by=1.2),
	labels = as.numeric(colnames(effVar_LB2)),col='black',cex=1.5)
text(y=rep(19,length(effVar_LB2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_LB2[1,])-1)),by=1.2),
	labels=as.numeric(effVar_LB2[1,]),col='gray55',cex=2)
text(y=rep(20,length(effVar_LB2[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_LB2[2,])-1)),by=1.2),
	labels=as.numeric(effVar_LB2[2,]),col='black',cex=2)

barplot(effVar_RBPS2[1:2,],ylim=c(0,24),col=c('gray55','black'),xlab="",ylab="",cex.lab=0.1, yaxt="n", xaxt="n")
text(y=22,x=0.2,labels="RNABPS, Number of variants for each relative position:",col='black',cex=2,font=2,adj=c(0,0))
text(y=rep(21,length(effVar_RBPS2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_RBPS2[1,])-1)),by=1.2),
	labels = as.numeric(colnames(effVar_RBPS2)),col='black',cex=1.5)
text(y=rep(19,length(effVar_RBPS2[1,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_RBPS2[1,])-1)),by=1.2),
	labels=as.numeric(effVar_RBPS2[1,]),col='gray55',cex=2)
text(y=rep(20,length(effVar_RBPS2[2,])),x=seq(from=0.7, to = (0.7+1.21*(length(effVar_RBPS2[2,])-1)),by=1.2),
	labels=as.numeric(effVar_RBPS2[2,]),col='black',cex=2)

dev.off()

####################################
#Score combination by regresion logistic
####################################

Formula = NULL
IndPos = which(dataResult$class_effect==1)
IndNeg = which(dataResult$class_effect==0)
AUC = NULL
ACC = NULL
LRT = NULL
AUC_Uni = NULL

for (i in 1:1000){
	RandIndPos = sample(IndPos, round(2*length(IndPos)/3,0),replace=F)
	RandIndNeg = sample(IndNeg, round(2*length(IndNeg)/3,0),replace=F)
	RandInd = c(RandIndPos,RandIndNeg)
	# dataTrain = dataResult[RandInd, c("class_effect", "MutInPBarea_Branch", "Delta_prob", "MutInPBarea_BPP", "Delta_BPP",
										# "MutInPBarea_SVM","Delta_SVM","MutInPBarea_LB", "Delta_LB", "delta_HSF", "MutInPBarea_RBPS" , "Delta_RBPS")]
	# dataTest = dataResult[-RandInd, c("class_effect", "MutInPBarea_Branch", "Delta_prob", "MutInPBarea_BPP", "Delta_BPP",
										# "MutInPBarea_SVM","Delta_SVM","MutInPBarea_LB", "Delta_LB", "delta_HSF", "MutInPBarea_RBPS" , "Delta_RBPS")]
	dataTrain = dataResult[RandInd, c("class_effect", "Delta_prob", "MutInPBarea_BPP", "Delta_BPP",
										"Delta_SVM", "Delta_LB", "delta_HSF", "Delta_RBPS")]
	dataTest = dataResult[-RandInd, c("class_effect", "Delta_prob", "MutInPBarea_BPP", "Delta_BPP",
										"Delta_SVM", "Delta_LB", "delta_HSF", "Delta_RBPS")]
	dataTrain = na.omit(dataTrain)
	dataTest = na.omit(dataTest)
	mod0<- glm(class_effect ~ 1, data = dataTrain, family = "binomial")

	# slm1 <- step(object = mod0, scope = class_effect ~ MutInPBarea_Branch + Delta_prob + MutInPBarea_BPP + Delta_BPP +
										# MutInPBarea_SVM + MutInPBarea_LB + Delta_SVM + delta_HSF + Delta_LB + MutInPBarea_RBPS + Delta_RBPS,
											# trace = 0, direction = "both")
	slm1 <- step(object = mod0, scope = class_effect ~ Delta_prob + MutInPBarea_BPP + Delta_BPP +
										Delta_SVM + delta_HSF + Delta_LB + Delta_RBPS,
											trace = 0, direction = "both")
	mod = as.character(slm1$formula)[3]
	modsplit = unlist(strsplit(mod, " + ", fixed=TRUE))
	modPaste = paste(modsplit[order(modsplit)],collapse=" + ")
	mod1 = glm(class_effect ~ MutInPBarea_BPP, data = dataTrain, family = "binomial")

	labels = dataTest$class_effect
	Proba = predict.glm(slm1, newdata = dataTest, type = "response")
	ProbaUni = predict.glm(mod1, newdata = dataTest, type = "response")
	classProba = rep(0,length(Proba))
	classProba[Proba>0.5] = 1
	t = table(classProba,dataTest$class_effect,dnn = c("ClassPred","class_effect"))
	acc = sum(diag(t))/sum(t)
	ACC = c(ACC,acc)

	score_GLM=data.frame(Proba,labels)
	pred_GLM=prediction(predictions=score_GLM[,1],labels=score_GLM[,2])
	auc_GLM<- performance(pred_GLM,"auc")
	auc_GLM<- unlist(slot(auc_GLM, "y.values"))
	auc_GLM=round(auc_GLM,3)
	AUC = c(AUC,auc_GLM)
	score_Uni=data.frame(ProbaUni,labels)
	pred_Uni=prediction(predictions=score_Uni[,1],labels=score_Uni[,2])
	auc_Uni<- performance(pred_Uni,"auc")
	auc_Uni<- unlist(slot(auc_Uni, "y.values"))
	auc_Uni=round(auc_Uni,3)
	AUC_Uni = c(AUC_Uni,auc_Uni)

	Formula = c(Formula,modPaste)

	testanova = anova(slm1,mod1,test="Chisq")
	LRT = c(LRT,as.numeric(unlist(testanova["Pr(>Chi)"]))[2])

}

t = table(as.factor(Formula))
tcount = as.numeric(t)
tnames = names(t[order(tcount,decreasing =T)])
tcount = tcount[order(tcount,decreasing =T)]
t.test(x = AUC[Formula == "MutInPBarea_BPP + MutInPBarea_Branch"], y = AUC_Uni[Formula == "MutInPBarea_BPP + MutInPBarea_Branch"])
tnames = tnames[1:10]
tcount = tcount[1:10]

pdf("ModelStability.pdf", width = 20, height = 16)

par(mar=c(35,20,1,1))

x <- barplot(tcount, xaxt="n",xlim=c(0,12), ylab="Model count", cex.lab=2, cex.axis=1.5)
text(cex=2, x=x, y=-2, tnames, xpd=TRUE, srt=45,pos=2)

Formula2 = factor(Formula,levels = tnames)

par(mar=c(35,10,1,1))
boxplot(ACC~Formula2, las=2,ylab ="Accuracy (Proba>0.5)",xlim=c(0,11), cex.lab=2, cex.axis=1.5)

par(mar=c(35,10,1,1),mgp=c(5,1,0))
boxplot(LRT~Formula2, las=2, ylab ="LRT with mutInPBarea_BPP as ref", xlim=c(0,11), cex.lab=2, cex.axis=1.5)
abline(h=0.01,col="red")

layout(matrix(c(1,2,2,2,2,2), nc=6))
par(mar=c(35,4,1,1))
boxplot(AUC_Uni,ylim = c(0.6,1),ylab ="AUC univarié (BPP)", cex.lab=2, cex.axis=1.5)
boxplot(AUC~Formula2, las=2,ylab ="AUC",xlim=c(0,10),ylim = c(0.6,1), cex.lab=2, cex.axis=1.5)
abline(h=median(AUC_Uni),col="red")

dev.off()
