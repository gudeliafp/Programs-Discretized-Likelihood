#####################################################################
#Program to obtain Bp and BCa Bootstrap confidence intervals, using #
#the Carbon Fibers Data presented in Kundu and Gupta (2006).        #
#The data represent the strength measured in GPA for simple carbon  #
#fibers that were tested under tension at gauge lengths 20 mm       #
#(Data Set I) and 10 mm  (Data Set II).                             #
#####################################################################

#Required libraries
library(alabama)
library(Rmpfr)
library(pracma)

#Required functions
source("FunctionReliabilityParameter.R",local=FALSE)
source("FunctionMeasuringInstrument.R",local=FALSE)
source("FunctionDiscretizedSample.R",local=FALSE)
source("FunctionProbabilityGEVXmin.R",local=FALSE)
source("FunctionMinusLogLikeGEV.R",local=FALSE)
source("FunctionConstraintsGEV.R",local=FALSE)
source("FunctionBootstrap95CIExample.R",local=FALSE)
source("FunctionProbabilityXminRmpfr.R",local=FALSE)
source("FunctionMinusLogLike.R",local=FALSE)
source("FunctionReliabilityParameterWeibull.R",local=FALSE)

#Discretized samples
DATASET1<-read.table("DataSet1.txt",header=FALSE,sep="")
DATASET2<-read.table("DataSet2.txt",header=FALSE,sep="")
x<-DATASET1$V1
y<-DATASET2$V1
Valh<-0.0005

#MLE of the parameters (mu1,gamma1,alpha1) of Weibull model (X)
Ui<-matrix(c(-1,0,0,
             0,1,0,
             0,0,1),ncol=3,byrow=TRUE)
CiX<-c(-min(x+Valh),0,0)
OptX<-constrOptim(c(1,1,1),FunctionMinusLogLike,NULL,ui=Ui,ci=CiX,mu=1e-04,control=list(),method="Nelder-Mead",
                  outer.iterations = 100, outer.eps = 1e-05,x,Valh)
EmvX<-OptX$par

#MLE of the parameters (mu2,gamma2,alpha2) of Weibull model (Y)
CiY<-c(-min(y+Valh),0,0)
OptY<-constrOptim(c(1,1,1),FunctionMinusLogLike,NULL,ui=Ui,ci=CiY,mu=1e-04,control=list(),method="Nelder-Mead",
                  outer.iterations = 100, outer.eps = 1e-05,y,Valh)
EmvY<-OptY$par

#MLE of reliability parameter
Emvdelta<-integrate(FunctionReliabilityParameterWeibull,lower=max(c(EmvX[1],EmvY[1])),upper=Inf,EmvX[1],EmvY[1],EmvX[2],EmvY[2],EmvX[3],EmvY[3])$value

#Summary of results: Point estimation
EMV<-matrix(c(EmvX,EmvY,Emvdelta),ncol=7)
colnames(EMV)<-c("EMVmu1","EMVgamma1","EMValpha1","EMVmu2","EMVgamma2","EMValpha2","EMVdelta")
#write.csv(EMV, file="EMV_Carbon_Fibers_Example.csv")

###############################################################
#An Alternative to Optimization: Likelihood function is based #
#on the von Mises's family of GEV distribution with parameters#
#(a,b,c); see Hirose, H., & Lai, T. L. (1997).                #
###############################################################

#MLE of the parameters (a,b,c) of the GEV distribution (X)
ValMleGEVX<-auglag(par=c(min(x+Valh),0.4,0.2),fn=FunctionMinusLogLikeGEV,hin=FunctionConstraintsGEV,control.outer=list(trace=FALSE),datos=x,precision=Valh)$par

#MLE of the parameters (a,b,c) of the GEV distribution (Y)
ValMleGEVY<-auglag(par=c(min(y+Valh),0.8,0.2),fn=FunctionMinusLogLikeGEV,hin=FunctionConstraintsGEV,control.outer=list(trace=FALSE),datos=y,precision=Valh)$par

#Reparameterization of Weibull distribution (X) in terms of the parameters
#(muX=mu1,gammaX,alphaX=alpha1), where gammaX=gamma1^(1/alpha1)
alphaX<-ValMleGEVX[3]^(-1)
gammaX<-ValMleGEVX[2]*alphaX
muX<-ValMleGEVX[1]-gammaX

#Reparameterization of Weibull distribution (Y) in terms of the parameters
#(muY=mu2,gammaY,alphaY=alpha2), where gammaY=gamma2^(1/alpha2)
alphaY<-ValMleGEVY[3]^(-1)
gammaY<-ValMleGEVY[2]*alphaY
muY<-ValMleGEVY[1]-gammaY

#MLE of reliability parameter
MLEdelta<-integrate(FunctionReliabilityParameter,lower=max(c(muX,muY)),upper=Inf,muX,muY,gammaX,gammaY,alphaX,alphaY)$value

################################################################

#Bootstrap confidence interval
M<-10000
ResultsBCI<-FunctionBootstrap95CIExample(M,x,y,ValMleGEVX,ValMleGEVY,length(x),length(y),Valh,Valh)

#Summary of results
SummaryBootstrapCI<-matrix(ResultsBCI,ncol=8)
colnames(SummaryBootstrapCI)<-c("LowerBp","UpperBp","z0hat","ahat","alpha1Jack","alpha2Jack","LowerBCa","UpperBCa")
#write.csv(SummaryBootstrapCI, file="SummaryBootstrapCI_Example.csv")

#Plot 900 x 700
BaseDatosDeltaStarExample<-read.csv("BootstrapDeltaStar_Example.csv")
DeltaStarExample<-BaseDatosDeltaStarExample[,2]
BaseDatosBootstrapCIExample<-read.csv("SummaryBootstrapCI_Example.csv")[,-1]
EjeX<-seq(0.10,0.40,0.05)
EjeY<-seq(0,2000,500)
hist(DeltaStarExample,xlab=expression(delta),ylab="",main="",cex.lab=1.23,cex.axis=1.20, axes = FALSE)
axis(1, at =EjeX,labels=EjeX, col.axis="black", las=1, cex.axis=1.2)
axis(2, at =EjeY,labels=EjeY, col.axis="black", las=0, cex.axis=1.2)
mtext(side=2,text="Frequency", line=2.5,cex=1.20)
points(x=c(BaseDatosBootstrapCIExample$LowerBp,BaseDatosBootstrapCIExample$UpperBp),
       y=c(0,0),pch="+",cex=1.85)
points(x=c(BaseDatosBootstrapCIExample$LowerBCa,BaseDatosBootstrapCIExample$UpperBCa),
       y=c(0,0),pch="*",cex=2.5)
