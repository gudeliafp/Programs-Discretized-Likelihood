####################################################################
#Program to obtain contour plots, Likelihood Ratio Test and Akaike #
#Information Criterion using the Carbon Fibers Data presented in   #
#Kundu and Gupta (2006).                                           #
#The data represent the strength measured in GPA for simple carbon #
#fibers that were tested under tension at gauge lengths 20 mm      #
#(Data Set I) and 10 mm  (Data Set II).                            #
####################################################################

#Required libraries
library(Rmpfr)

#Required functions
source("FunctionProbabilityXminRmpfr.R",local=FALSE)
source("FunctionMinusLogLike.R",local=FALSE)
source("FunctionMinusLogLikeThresholdShape.R",local=FALSE)
source("FunctionMinusLogLikeXY.R",local=FALSE)

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

#Likelihood contours (Threshold,Shape)
ValThresholdXY<-seq(0,min(y+Valh),length.out=200)
ValShapeXY<-seq(1,6,length.out=200)
v0ThresholdShapeX<-EmvX[2]
v0ThresholdShapeY<-EmvY[2]
LogVeroThresholdShapeX<-matrix(rep(0,length(ValThresholdXY)*length(ValShapeXY)),ncol=length(ValShapeXY))
LogVeroThresholdShapeY<-matrix(rep(0,length(ValThresholdXY)*length(ValShapeXY)),ncol=length(ValShapeXY))
for(i in 1:length(ValThresholdXY))
{
  for(j in 1:length(ValShapeXY))
  {
    if(ValThresholdXY[i]<=min(x+Valh))
    {
      OptThresholdShapeX<-nlm(FunctionMinusLogLikeThresholdShape,v0ThresholdShapeX,x,Valh,ValThresholdXY[i],ValShapeXY[j])
      v0ThresholdShapeX<-OptThresholdShapeX$estimate
      LogVeroThresholdShapeX[i,j]<--OptThresholdShapeX$minimum
    }
    else
    {
      LogVeroThresholdShapeX[i,j]<--10^100
    }
  }
}
for(i in 1:length(ValThresholdXY))
{
  for(j in 1:length(ValShapeXY))
  {
    if(ValThresholdXY[i]>1 & ValShapeXY[j]<=5)
    {
      OptThresholdShapeY<-nlm(FunctionMinusLogLikeThresholdShape,v0ThresholdShapeY,y,Valh,ValThresholdXY[i],ValShapeXY[j])
      v0ThresholdShapeY<-OptThresholdShapeY$estimate
      LogVeroThresholdShapeY[i,j]<--OptThresholdShapeY$minimum
    }
    else
    {
      LogVeroThresholdShapeY[i,j]<--10^100
    }
  }
}
RelativaThresholdShapeY<-exp(LogVeroThresholdShapeY-max(LogVeroThresholdShapeY))
RelativaThresholdShapeX<-exp(LogVeroThresholdShapeX-max(LogVeroThresholdShapeX))
W<-c(0,0.375,0.687,1,1.3125,1.607,1.9015)
#plot 900 x 700 jpeg
contour(ValThresholdXY,ValShapeXY,RelativaThresholdShapeX,
        levels=c(1,0.5,0.1,0.05),
        labels=c(1,0.5,0.1,0.05),
        labcex=1,lty=1,lwd=2,
        xlab=expression(mu),
        ylab=expression(alpha),
        cex.lab=1.5,yaxs="i",xaxs="i",
        cex.axis=1.25,xlim=c(0,2),ylim=c(1,6),xaxt="n")
axis(1, at=W,labels=W, col.axis="black", las=1, cex.axis=1.25)
points(EmvX[1],EmvX[3],pch=8,cex=0.75,col=1)
contour(ValThresholdXY,ValShapeXY,RelativaThresholdShapeY,
        levels=c(1,0.5,0.1,0.05),
        labels=c(1,0.5,0.1,0.05),
        labcex=1,lty=2,lwd=2,
        xlab=expression(mu),
        ylab=expression(alpha),
        cex.lab=1.5,yaxs="i",xaxs="i",
        cex.axis=1.25,xlim=c(0,2),ylim=c(1,6),add=TRUE)
points(EmvY[1],EmvY[3],pch=8,cex=0.75,col=1)
legend("topright",c("X","Y"),title="Data set",lty=c(1,2),lwd=2,col=1,bty="n",cex=1.25,y.intersp=0.75,adj=1)
mtext(side=1,text=expression("                                                      " * x[(1)]+h), line=2.5,cex=1)
mtext(side=1,text=expression("                                                                                                                                                              " * y[(1)]+h), line=2.5,cex=1)

#likelihood ratio test
UiXY<-matrix(c(-1,0,0,0,
               0,1,0,0,
               0,0,1,0,
               0,0,0,1),ncol=4,byrow=TRUE)
z<-min(c(min(x),min(y)))
CiXY<-c(-(z+Valh),0,0,0)
OptXY<-constrOptim(c(1,1,1,1),FunctionMinusLogLikeXY,NULL,ui=UiXY,ci=CiXY,mu=1e-04,control=list(),method="Nelder-Mead",
                   outer.iterations = 100, outer.eps = 1e-05,x,y,Valh)
EmvXY<-OptXY$par
MaxLogVero<--(OptX$value+OptY$value)
MaxLogVeroH0<--OptXY$value
D<--2*(MaxLogVeroH0-MaxLogVero)
Pvalor<-1-pchisq(D, df=2)
c(D,Pvalor)

#The Akaike information criterion value using correct likelihood
AICModel6Par<-2*6-2*MaxLogVero
AICModel4Par<-2*4-2*MaxLogVeroH0
c(AICModel6Par,AICModel4Par)

#The Akaike information criterion value using MLE's of pseudolikelihood
MaxLogVeroXPseudo<-FunctionMinusLogLike(c(1.312,248.3652,4.6344),x,Valh)
MaxLogVeroYPseudo<-FunctionMinusLogLike(c(1.312,86.9579,4.6344),y,Valh)
MaxLogVeroXYPseudo<--(MaxLogVeroXPseudo+MaxLogVeroYPseudo)
AICModelPseudo<-2*3-2*MaxLogVeroXYPseudo
AICModelPseudo


