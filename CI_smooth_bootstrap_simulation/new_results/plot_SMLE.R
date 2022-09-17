   B<-read.table("SMLE.txt")
   x<-B[,1]
   y<-B[,2]
   f <- function(x) {1-exp(-x)}/{1-exp(-2)}
   x0 <-seq(0,max(x),by=0.01)
   y0<-f(x0)
   plot(c(-1000,-1000),xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), main= "",ylab="",xlab="",bty="n",las=1)
   lines(x,y,lwd=2,col="blue")
   #lines(x0,y0,lwd=2,lty=2,col="red")  
	
	
	
	

    