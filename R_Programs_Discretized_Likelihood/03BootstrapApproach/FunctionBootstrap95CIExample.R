FunctionBootstrap95CIExample<-function(MBootstrap,ValXobs,ValYobs,EMVGEVX,EMVGEVY,nBootstrapX,mBootstrapY,ValhX,ValhY)
{
  InstrumentX<-FunctionMeasuringInstrument(ValhX,EMVGEVX[1],EMVGEVX[2],EMVGEVX[3])
  InstrumentY<-FunctionMeasuringInstrument(ValhY,EMVGEVY[1],EMVGEVY[2],EMVGEVY[3])
  MiddlePointX<-InstrumentX[[1]]
  ExtremesX<-InstrumentX[[2]]
  MiddlePointY<-InstrumentY[[1]]
  ExtremesY<-InstrumentY[[2]]
  BootSampleMatrixX<-c()
  BootSampleMatrixY<-c()
  deltaStar<-c()
  for(i in 1:MBootstrap)
  {
    XobsSimBoot<-FunctionDiscretizedSample(EMVGEVX[1],EMVGEVX[2],EMVGEVX[3],nBootstrapX,MiddlePointX,ExtremesX)
    YobsSimBoot<-FunctionDiscretizedSample(EMVGEVY[1],EMVGEVY[2],EMVGEVY[3],mBootstrapY,MiddlePointY,ExtremesY)
    BootSampleMatrixX<-rbind(BootSampleMatrixX,XobsSimBoot)
    BootSampleMatrixY<-rbind(BootSampleMatrixY,YobsSimBoot)
    MatrizXYobsSimBoot<-meshgrid(XobsSimBoot,YobsSimBoot)
    deltaStar[i]<-sum(MatrizXYobsSimBoot$Y<MatrizXYobsSimBoot$X)/(nBootstrapX*mBootstrapY)
  }
  #Save bootstrap samples
  write.csv(BootSampleMatrixX, file="BootstrapSampleX_Example.csv")
  write.csv(BootSampleMatrixY, file="BootstrapSampleY_Example.csv")
  write.csv(deltaStar, file="BootstrapDeltaStar_Example.csv")
  #The percentile method
  CIBp<-quantile(deltaStar,probs=c(0.025,0.975))
  LCIBp<-round(as.numeric(CIBp[1]),4)
  UCIBp<-round(as.numeric(CIBp[2]),4)
  #The BCa method
  ValnX<-length(ValXobs)
  ValmY<-length(ValYobs)
  MatrizXYobs<-meshgrid(ValXobs,ValYobs)
  MatrizYmenorX<-MatrizXYobs$Y-MatrizXYobs$X
  IndicaYmenorX<-(MatrizYmenorX<0)
  SumRow<-rowSums(IndicaYmenorX)
  SumCol<-colSums(IndicaYmenorX)
  SumTotal<-sum(SumRow)
  deltahat<-SumTotal/(ValnX*ValmY)
  deltaSumRow<-(SumTotal-SumRow)/((ValmY-1)*ValnX)
  deltaSumCol<-(SumTotal-SumCol)/((ValnX-1)*ValmY)
  deltahatJack<-c(deltaSumRow,deltaSumCol)
  centereddeltahatJack<-mean(deltahatJack)-deltahatJack
  ahat<-(1/6)*sum(centereddeltahatJack^3)/((sum((centereddeltahatJack^2)))^(3/2))
  z0hat<-qnorm(sum(deltaStar<deltahat)/MBootstrap,0,1)
  z1<-qnorm(0.025,0,1)
  z2<-qnorm(0.975,0,1)
  V1Jack<-z0hat+(z0hat+z1)/(1-ahat*(z0hat+z1))
  V2Jack<-z0hat+(z0hat+z2)/(1-ahat*(z0hat+z2))
  alpha1Jack<-pnorm(V1Jack,0,1)
  alpha2Jack<-pnorm(V2Jack,0,1)
  CIBCa<-quantile(deltaStar,probs=c(alpha1Jack,alpha2Jack))
  LCIBCa<-round(as.numeric(CIBCa[1]),4)
  UCIBCa<-round(as.numeric(CIBCa[2]),4)
  return(c(LCIBp,UCIBp,z0hat,ahat,alpha1Jack,alpha2Jack,LCIBCa,UCIBCa))
}