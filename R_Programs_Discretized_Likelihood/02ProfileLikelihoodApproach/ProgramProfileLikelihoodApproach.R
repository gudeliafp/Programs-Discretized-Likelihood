#####################################################################
#Program to obtain profile likelihood inferences using the Carbon   #
#Fibers Data presented in Kundu and Gupta (2006).                   #
#The data represent the strength measured in GPA for simple carbon  #
#fibers that were tested under tension at gauge lengths 20 mm       #
#(Data Set I) and 10 mm  (Data Set II).                             #
#####################################################################

#Required libraries
#library(alabama)
library(Rmpfr)
#library(pracma)

#Required functions
source("FunctionProbabilityXminRmpfr.R",local=FALSE)
source("FunctionMinusLogLike.R",local=FALSE)
source("FunctionMinusLogLikeXY.R",local=FALSE)
source("FunctionMinusLogLikeThresholdXY.R",local=FALSE)
source("FunctionMinusLogLikeReliabilityXY.R",local=FALSE)
source("FunctionMinusPseudoLogLikeReliabilityXY.R",local=FALSE)

#Discretized samples
DATASET1<-read.table("DataSet1.txt",header=FALSE,sep="")
DATASET2<-read.table("DataSet2.txt",header=FALSE,sep="")
x<-DATASET1$V1
y<-DATASET2$V1
Valh<-0.0005

#MLE of the parameters (mu,gamma1,gamma2,alpha) and delta=gamma1/(gamma1+gamma2)
UiXY<-matrix(c(-1,0,0,0,
               0,1,0,0,
               0,0,1,0,
               0,0,0,1),ncol=4,byrow=TRUE)
z<-min(c(min(x),min(y)))
CiXY<-c(-(z+Valh),0,0,0)
OptXY<-constrOptim(c(0.5,1,1,1),FunctionMinusLogLikeXY,NULL,ui=UiXY,ci=CiXY,mu=1e-04,control=list(),method="Nelder-Mead",
                   outer.iterations = 100, outer.eps = 1e-05,x,y,Valh)
EmvXY<-OptXY$par
EMVXYdelta<-EmvXY[2]/(EmvXY[2]+EmvXY[3])

#Summary of results: Point estimation
EMVprofile<-matrix(c(EmvXY,EMVXYdelta),ncol=5)
colnames(EMVprofile)<-c("EMVmu","EMVgamma1","EMVgamma2","EMValpha","EMVdelta")
#write.csv(EMVprofile, file="EMVprofile_Carbon_Fibers_Example.csv")

#Profile likelihood function of mu
ValThreshold<-c(seq(0,z-Valh,length.out=800),seq(z-Valh,z+Valh-2*Valh/10,length.out=20)[-1])
UiThresholdXY<-matrix(c(1,0,0,
                        0,1,0,
                        0,0,1),ncol=3,byrow=TRUE)
CiThresholdXY<-c(0,0,0)
LogVeroThresholdXY<-c()
v0ThresholdXY<-EmvXY[2:4]
for(i in 1:length(ValThreshold))
{
  OptThresholdXY<-constrOptim(v0ThresholdXY,FunctionMinusLogLikeThresholdXY,NULL,ui=UiThresholdXY,ci=CiThresholdXY,mu=1e-04,control=list(),
                              method="Nelder-Mead",outer.iterations = 100, outer.eps = 1e-05,x,y,Valh,ValThreshold[i])
  v0ThresholdXY<-OptThresholdXY$par
  LogVeroThresholdXY[i]<--OptThresholdXY$value
}
RelativaPerfilThresholdXY<-exp(LogVeroThresholdXY-max(LogVeroThresholdXY))

#Profile likelihood interval for mu
vc<-0.1465
indice<-1
while(RelativaPerfilThresholdXY[indice]<vc)
{
  indice<-indice+1
}
LimiteInferiorThresholdXY<-ValThreshold[indice]

indice<-length(ValThreshold)
while(RelativaPerfilThresholdXY[indice]<vc)
{
  indice<-indice-1
}
LimiteSuperiorThresholdXY<-ValThreshold[indice]
c(LimiteInferiorThresholdXY,LimiteSuperiorThresholdXY)

#plot 900 x 700 jpeg
WXY<-c(0,0.375,0.687,1,1.3125)
plot(ValThreshold,RelativaPerfilThresholdXY,type="l",lty=1,lwd=2,col=1,xlab=expression(mu),ylab="",
     cex.lab=1.5,cex.axis=1.25,yaxs="i",ylim=c(0,1),xlim=c(0,z+Valh),xaxt="n")
axis(1, at=WXY,labels=WXY, col.axis="black", las=1, cex.axis=1.25)
mtext(side=2,text="Relative profile likelihood", line=2.5,cex=1.4)
mtext(side=1,text=expression("                                                                                                                                                                   " * tilde(mu)+h), line=2.5,cex=1)
segments(LimiteInferiorThresholdXY,vc,LimiteSuperiorThresholdXY,vc,col=1,lty=1)
text(EmvXY[1],0.1481,"*",cex=2.5,col=1)
text(0.91,0.19,paste("95% CI: [",round(LimiteInferiorThresholdXY,4),",",round(LimiteSuperiorThresholdXY,4),"]"),cex=1.25)

#Pseudo-MLE of the parameters (gamma1,gamma2,alpha)
UiXYEstimada<-matrix(c(1,0,0,
                       0,1,0,
                       0,0,1),ncol=3,byrow=TRUE)
CiXYEstimada<-c(-(z+Valh),0,0)
OptXYEstimada<-constrOptim(c(1,1,1),FunctionMinusLogLikeThresholdXY,NULL,ui=UiXYEstimada,ci=CiXYEstimada,mu=1e-04,control=list(),method="Nelder-Mead",
                           outer.iterations = 100, outer.eps = 1e-05,x,y,Valh,z)
EmvXYEstimada<-OptXYEstimada$par

#Profile likelihood functions of delta=P(Y<X)=gamma1/(gamma1+gamma2)
ValReliability<-seq(0.1,0.45,length.out=150)
ValReliabilityEstimada<-seq(0.1,0.45,length.out=150)
UiReliabilityXY<-matrix(c(-1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE)
UiReliabilityXYEstimada<-matrix(c(1,0,0,1),ncol=2,byrow=TRUE)
CiReliabilityXY<-c(-(z+Valh),0,0)
CiReliabilityXYEstimada<-c(0,0)
LogVeroReliabilityXY<-c()
LogVeroReliabilityXYEstimada<-c()
v0ReliabilityXY<-c(EmvXY[1:2],EmvXY[4])
v0ReliabilityXYEstimada<-c(EmvXYEstimada[1],EmvXYEstimada[3])
for(i in 1:length(ValReliability))
{
  OptReliabilityXY<-constrOptim(v0ReliabilityXY,FunctionMinusLogLikeReliabilityXY,NULL,ui=UiReliabilityXY,ci=CiReliabilityXY,mu=1e-04,control=list(),
                                method="Nelder-Mead",outer.iterations = 100, outer.eps = 1e-05,x,y,Valh,ValReliability[i])
  v0ReliabilityXY<-OptReliabilityXY$par
  LogVeroReliabilityXY[i]<--OptReliabilityXY$value
}
for(i in 1:length(ValReliabilityEstimada))
{
  OptReliabilityXYEstimada<-constrOptim(v0ReliabilityXYEstimada,FunctionMinusPseudoLogLikeReliabilityXY,NULL,ui=UiReliabilityXYEstimada,ci=CiReliabilityXYEstimada,mu=1e-04,control=list(),
                                        method="Nelder-Mead",outer.iterations = 100, outer.eps = 1e-05,x,y,Valh,ValReliabilityEstimada[i],z)
  v0ReliabilityXYEstimada<-OptReliabilityXYEstimada$par
  LogVeroReliabilityXYEstimada[i]<--OptReliabilityXYEstimada$value
}
RelativaPerfilReliabilityXY<-exp(LogVeroReliabilityXY-max(LogVeroReliabilityXY))
RelativaPerfilReliabilityXYEstimada<-exp(LogVeroReliabilityXYEstimada-max(LogVeroReliabilityXYEstimada))

#plot 900 x 700 jpeg
LIBp<-0.1527
LSBp<-0.3085
LIBCa<-0.1598
LSBCa<-0.3177
plot(ValReliability,RelativaPerfilReliabilityXY,type="l",lty=1,lwd=2,col=1,xlab=expression(delta),ylab="",
     cex.lab=1.5,cex.axis=1.25,yaxs="i",ylim=c(0,1),xlim=c(0.1,0.45))
mtext(side=2,text="Relative profile likelihood", line=2.5,cex=1.4)
points(ValReliabilityEstimada,RelativaPerfilReliabilityXYEstimada,type="l",lty=2,lwd=2,col=1)
legend(x=0.33,y=1,c(expression(R[p](delta)),expression(R[ep](delta*";"~mu*"="*tilde(mu)))),
       lty=c(1,2),lwd=2,col=1,bty="n",cex=1.4,y.intersp=0.5)
points(x=c(LIBp,LSBp),y=c(vc,vc),pch="+",cex=1.85)
points(x=c(LIBCa,LSBCa),y=c(vc,vc),pch="*",cex=2.5)

#Profile likelihood interval for delta=P(Y<X)=gamma1/(gamma1+gamma2)
vc<-0.1465
indice<-1
while(RelativaPerfilReliabilityXY[indice]<vc)
{
  indice<-indice+1
}
LimiteInferiorReliabilityXY<-ValReliability[indice-1]

indice<-length(ValReliability)
while(RelativaPerfilReliabilityXY[indice]<vc)
{
  indice<-indice-1
}
LimiteSuperiorReliabilityXY<-ValReliability[indice+1]
c(LimiteInferiorReliabilityXY,LimiteSuperiorReliabilityXY)


