# Simple Bayesian model to infer bias of a coin.

source('DBDA2E-utilities.R')
# Resolution of all parameters. How fine-grained is the grid?
res = 1000 
# Set Unit interval based on resolution
UnitInt = seq(0.01,0.99,length.out=res) 

# Parameters for a beta distribution
alpha = 1 
beta = c(0.5,0.5) 

a = alpha*beta[1]
b = alpha*beta[2]
## omega = 0.5
## kappa = 4

## a = omega * (kappa - 2) + 1 
## b = (1 - omega) * (kappa - 2) + 1 

# the 'Data'
N = 20 
z = 17
N2 = 1
z2 = 1

# Prior distribution over theta
Theta = UnitInt # Possible theta values
dThetaBeta = dbeta(Theta,a,b) # The prior on theta as a density
pTheta = dThetaBeta*abs(Theta[2]-Theta[1]) # convert density to a mass function

pDataGivenTheta = Theta^z * (1 - Theta)^(N-z) # Likelihood
pData = sum( pDataGivenTheta * pTheta) # In place of an integral 
pThetaGivenData = pDataGivenTheta * pTheta / pData # Bayes Theorem!

pDataGivenTheta = Theta^z2 * (1 - Theta)^(N2-z2) # Likelihood
pData = sum( pDataGivenTheta * pTheta) # In place of an integral 
pThetaGivenData2 = pDataGivenTheta * pTheta / pData # Bayes Theorem!

openGraph(width=5,height=7)
# Plot the results.
layout( matrix( c( 1,2,3 ) ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par( mar=c(3,3,1,0) , mgp=c(2,0.7,0) , mai=c(0.5,0.5,0.3,0.1) ) # margins

plotType = "l" # Plot Bars

plot( Theta , dThetaBeta, type=plotType , lwd=3,
     xlab=bquote(theta) , ylab=bquote(p(theta)),
     yaxt='n',
     xlim=c(0.02,0.98),
     ylim=c(0,8),
     main="Prior" , col="plum1" )
plot( Theta, pThetaGivenData2, type=plotType , lwd=3,
     yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(" * theta * ")" ) ,
     main="Posterior for N=1, z=1" , col="plum1" )
abline(lwd=3, v=weighted.mean(Theta,pThetaGivenData2), col='plum1')
plot( Theta, pThetaGivenData, type=plotType , lwd=3,
     yaxt='n', xlab=bquote(theta) , ylab=bquote( "p(" * theta * ")" ) ,
     main="Posterior for N=20, z=17" , col="plum1" )
abline(lwd=3, v=weighted.mean(Theta,pThetaGivenData), col='plum1')
weighted.mean(Theta,pThetaGivenData)

saveGraph(file="SimpleBayes",type="pdf")
