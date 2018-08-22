rm(list=ls())  # Careful! This clears all of R's memory!
# Simple Bayesian model to infer Coin Biases 
source('DBDA2E-utilities.R')
# Resolution of all parameters. How fine-grained is the grid?
res = 50 
# Set Unit interval based on resolution
UnitInt = seq(1/res^2,1-1/res^2,length.out=res) 

# Parameters for priot beta distribution
omega = 0.5
kappa = 1.9 

a = omega * (kappa - 2) + 1 
b = (1 - omega) * (kappa - 2) + 1 

# the Data.  
mixedD =matrix(c(20,17,20,5,20,11,1,1),ncol=2, byrow=TRUE)
allD = matrix(c(20,20,20,0,20,20,1,1),ncol=2, byrow=TRUE)

# Prior distribution over theta
nThetas = 4
Theta = UnitInt # Possible theta values
pThetas = matrix(nrow=nThetas, ncol=res)
for (i in 1:nThetas) {
    pThetas[i,] = dbeta(Theta,a,b)*abs(Theta[2]-Theta[1]) # convert density to a mass function
}

pDataGivenThetas = matrix(nrow=nThetas, ncol=res)
pThetasGivenData = matrix(nrow=nThetas, ncol=res)

## data = mixedD
data = allD
for (i in 1:nThetas) {
    N = data[i,1]
    z = data[i,2]
    pDataGivenThetas[i,] = Theta^z * (1 - Theta)^(N-z) 
    pData = sum( pDataGivenThetas[i,] * pThetas[i,]) 
    pThetasGivenData[i,] = pDataGivenThetas[i,] * pThetas[i,] / pData # Bayes Theorem!
}


# Calculate the pseudo combined posterior
pseudoN = apply(data[1:nThetas-1,],2,sum)[1]
pseudoZ = apply(data[1:nThetas-1,],2,sum)[2]

pTheta = dbeta(Theta,a,b)*abs(Theta[2]-Theta[1])
pDataGivenTheta = Theta^pseudoZ * (1 - Theta)^(pseudoN-pseudoZ) 
pData = sum( pDataGivenTheta * pTheta) 
pseudoPrior = pDataGivenTheta * pTheta / pData # Bayes Theorem!

# Now update the Prior on just the last case!
pDataGivenTheta2 = Theta 
pData2 = sum( pDataGivenTheta2 * pseudoPrior) 
pseudoPosterior = pDataGivenTheta2 * pseudoPrior / pData2

## openGraph(width=5,height=7)
## # Plot the results.
layout(matrix( 1:((nThetas*2)+2), nrow=nThetas+1, ncol=2, byrow=TRUE))
## layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par( mar=c(3,3,1,0) , mgp=c(2,0.7,0) , mai=c(0.5,0.5,0.3,0.1) ) # margins

plotType = "l" # Plot Bars
for (i in 1:nThetas) {
    plot( Theta , pThetas[i,] , type=plotType , lwd=3,
         yaxt='n', xlab=bquote(theta) , ylab=bquote(p(theta)),
         main="Prior" , col="plum1" )
    plot( Theta, pThetasGivenData[i,], type=plotType , lwd=3,
         yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(" * theta * "|D)" ) ,
         main=paste("Posterior, N =",data[i,1]," z =",data[i,2]), col="plum1" )
}
plot( Theta , pseudoPrior , type=plotType , lwd=3,
     yaxt='n', xlab=bquote(theta) , ylab=bquote(p(theta)),
     main=paste("Prior, N =",pseudoN," z =",pseudoZ), col="plum1" )
plot( Theta, pseudoPosterior, type=plotType , lwd=3,
     yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(" * theta * "|D)" ) ,
     main=paste("Posterior, N =",1," z =",1), col="plum1" )


## plot( Theta , pTheta , type=plotType , lwd=3,
##      yaxt='n', xlab=bquote(theta) , ylab=bquote(p(theta)),
##      main="Prior" , col="plum1" )
## plot( Theta , pDataGivenTheta , type=plotType , lwd=3,
##      yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(D|" * theta * ")" ) ,
##      main="Likelihood" , col="plum1" )
## plot( Theta, pThetaGivenData, type=plotType , lwd=3,
##      yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(" * theta * "|D)" ) ,
##      main="Posterior" , col="plum1" )
## abline(lwd=3, v=weighted.mean(Theta,pThetaGivenData), col='plum1')
## weighted.mean(Theta,pThetaGivenData) 
## saveGraph(file="SimpleBayes",type="pdf")
