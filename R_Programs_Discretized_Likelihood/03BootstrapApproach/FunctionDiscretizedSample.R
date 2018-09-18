FunctionDiscretizedSample<-function(a0,b0,c0,SampleSize,MiddlePointR,ExtremesR)
{
  p<-runif(SampleSize,0,1)
  Xsim<-a0+(b0/c0)*(((-log(1-p))^c0)-1)
  Xobs<-as.numeric(as.vector(cut(Xsim,ExtremesR,labels=MiddlePointR)))
  return(Xobs)
}