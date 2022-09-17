	library(Rcpp)
	sourceCpp("Banerjee_Wellner_simulation.cpp")
  	output <- ComputeIntervals(N=500,NumB=1000)
   
	B1 <- output$MLE
	C <- output$CI_Banerjee_Wellner
	D <- output$percentages
	
	B<-read.table("MLE.txt")

	f <- function(x) {(1-exp(-x))/(1-exp(-2))}
	x0<-seq(0,2,by=0.01)
	y0<-f(x0)
	
  
   x1<-B[,1]
   z1<-B[,2]
   x2<-C[,1]
   t1<-C[,2]
   u1<-C[,3]
   
   plot(c(-100,-100),xlim=c(0,2), ylim=c(0,max(u1)), main= "",ylab="",xlab="",bty="n",las=1)
   lines(x1,z1,lwd=2,type ="s",col="red")
   lines(x0,y0,lwd=2,lty=2)
   segments(x2,t1,x2,u1)


   x1<-D[,1]
   y1<-D[,2]
   	  	
   plot(c(-10000,-10000),xlim=c(0.0,2), ylim=c(0.0,0.15), main= "", ylab="",xlab="",bty="n",las=1)
   lines(x1,y1,lty=1)
   lines(c(0,2),c(0.05,0.05),col="red")
