### WARNING: Run PsychometricFunction1.R first and do not close the R window. 

#### PLOT FOR EXERCISE 12.1.2 

dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
	{
	scale <- seq(x[i,1],x[i,nstim[i]], by=.1)
	plot(x[i,],rprop[i,],main=paste("Subject",as.character(i)),xlab="",ylab="",pch=15,col="dark grey",ylim=c(0,1),yaxt="n",xaxt="n")
	for (g in 1:20)
	  {
	  lines(scale,F3(scale,i,g),type="l",col="light grey")
	  }
	lines(scale,F1(scale,i),type="l")
	if (i==1 | i==5) {
		axis(2,las=1,yaxp=c(0,1,2))
		#axis(2,at=0.84,las=1)
		}
	if (i>4) axis(1)
	}
	mtext("Proportion 'Long' Response",side=2,line=2,outer=T,cex=1.4)
	mtext("Test Interval (ms)",side=1,outer=T,line=3,cex=1.4)


#### PLOT FOR EXERCISE 12.1.4 

dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
{
	plot(density(JND[,i]), col="dark grey", type="h", main=paste("Subject",as.character(i)), ylab="", xlab="", xlim=c(10,100), ylim=c(0,0.12), las=1,yaxt="n",xaxt="n")
	par(new=F)
	if (i==1 | i==5) axis(2,las=1, yaxp=c(0,0.12,3))
	if (i>4) axis(1)
}	
	mtext("Posterior Density",side=2,line=2,outer=T,cex=1.4)
	mtext("JND (ms)",side=1,line=3,outer=T,cex=1.4)

	