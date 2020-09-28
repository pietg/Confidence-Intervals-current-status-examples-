	library(Rcpp)
	A<-read.table("hepatitis-A-data.txt") 
	sourceCpp("curstat_bootstrap_HepatitisA.cpp")
  
	output <- ComputeIntervals(A)
   
	B <- output$MLE
	C <- output$SMLE
	D <- output$CI_SMLE
  
   x1<-B[,1]
   z1<-B[,2]
   x2<-C[,1]
   z2<-C[,2]
   x3<-D[,1]
   t1<-D[,2]
   u1<-D[,3]
   
   plot(c(-100,-100),xlim=c(0,max(x1)), ylim=c(0,max(u1)), main= "",ylab="",xlab="",bty="n",las=1)
   lines(x1,z1,lwd=2,type ="s",col="red")
   lines(x2, z2,lwd=2,col="blue")
   segments(x3,t1,x3,u1)

   
 