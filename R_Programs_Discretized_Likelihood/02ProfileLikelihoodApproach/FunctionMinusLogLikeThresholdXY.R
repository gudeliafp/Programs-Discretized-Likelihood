FunctionMinusLogLikeThresholdXY<-function(VecPar,datosX,datosY,precision,Threshold)
{
  ScaleX<-VecPar[1]
  ScaleY<-VecPar[2]
  Shape<-VecPar[3]
  MenosLogVeroX<-FunctionMinusLogLike(c(Threshold,ScaleX,Shape),datosX,precision)
  MenosLogVeroY<-FunctionMinusLogLike(c(Threshold,ScaleY,Shape),datosY,precision)
  MenosLogVero<-MenosLogVeroX+MenosLogVeroY
  return(MenosLogVero)
}
