FunctionMinusLogLikeGEV<-function(VecPar,datos,precision)
{
  aPar<-VecPar[1]
  bPar<-VecPar[2]
  cPar<-VecPar[3]
  Xmin<-min(datos)
  m<-sum(datos==Xmin)
  Xm<-sort(datos)[-c(1:m)]
  Xmmenosh<-Xm-precision
  Xmmash<-Xm+precision
  MenosLogVeroXmin<--m*log(FunctionProbabilityGEVXmin(VecPar,Xmin,precision))
  Am<-mpfr(-(1+cPar*(Xmmenosh-aPar)/bPar)^(1/cPar),500)
  Bm<-mpfr(-(1+cPar*(Xmmash-aPar)/bPar)^(1/cPar),53)
  MenosLogVeroXm<--sum(log(exp(Am)-exp(Bm)))
  MenosLogVero<-MenosLogVeroXmin+MenosLogVeroXm
  return(as.numeric(MenosLogVero))
}

