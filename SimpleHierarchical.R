# Simplest hierarchical Bayesian model to infer bias of a single coin.
res = 100 # Resolution of all parameters. How fine-grained is the grid?
UnitInt = seq(1/res^2,1-1/res^2,length.out=res) # Set Unit interval based on resolution

# Parameters for the prior distribution on omega:
A = 2 
B = 2 

# Parameters for the prior distribution on theta
# Kappa is fixed, while omega is a variable
kappa = 10
Omega = UnitInt 
# the 'Data'
N = 20
z = 17

Theta = UnitInt 

dOmega = dbeta(Omega,A,B) # Prior density on Omega
pOmega = dOmega * abs(Omega[1]-Omega[2]) # Convert to mass function 

# Generate the conditional probabilities of theta given each parameter value
# of omega as a matrix
pThetaGivenOmega = matrix (nrow=res, ncol=res)
for (i in 1:res) {
    omega = Omega[i]
    a = omega * (kappa - 2) + 1
    b = (1 - omega) * (kappa - 2) + 1 # For kappa > 2 only!
    dThetaGivenOmega = dbeta(Theta,a,b)
    pThetaGivenOmega[i,] = dThetaGivenOmega * abs(Omega[1]-Omega[2])
}

# Calculate the joint probability by POINTWISE multiplication.
# This is NOT matrix multiplication.
pOmegaTheta = pThetaGivenOmega * pOmega

# Calculate the marginal priors by summing out over the columns.
# i.e. do a weighted sum.
pTheta = colSums(pThetaGivenOmega * pOmega)

pDataGivenTheta = Theta^z * (1 - Theta)^(N-z) # Likelihood

# Calculate the Posterior joint distribution
pData = choose(N,z) * sum(pOmegaTheta * pDataGivenTheta)
pData
pOmegaThetaGivenData =  pOmegaTheta * pDataGivenTheta / pData

pThetaPost = colSums(pOmegaThetaGivenData)
pOmegaPost = rowSums(pOmegaThetaGivenData)

# Plot!
layout( matrix( c(1,2,3,0,4,5,6,0) ,nrow=4 ,ncol=2 ,byrow=TRUE ) ) # 3x1 panels
par( mar=c(3,3,1,0) , mgp=c(2,0.7,0) , mai=c(0.5,0.5,0.3,0.1) ) # margins

contour(Omega, Theta, pOmegaTheta, xlab=expression(omega),
        drawlabels = FALSE, main='Prior', 
        ylab=expression(theta), col="violet")
plot( Theta, pTheta, xlab=expression(theta),lwd=3,col="plum1",
     yaxt='n', ylab=expression(paste("p(",theta,")")), type='l')
weighted.mean(Theta,pTheta)
abline(v=weighted.mean(Theta,pTheta),  lty=2, col='plum1',lwd=3)
plot( Omega, pOmega, xlab=expression(omega),lwd=3, col="plum1",
     yaxt='n', ylab=expression(paste("p(",omega,")")), type='l')

contour(Omega, Theta, pOmegaThetaGivenData, xlab=expression(omega),
        drawlabels = FALSE, main='Posterior',
        ylab=expression(theta), col="violet")
plot(Theta, pThetaPost, type='l', lwd=3, col="plum1",
     yaxt='n', ylab=expression(paste("p(",theta,"|D)")) )
abline(v=weighted.mean(Theta,pThetaPost), col='plum1',lwd=3)
weighted.mean(Theta,pThetaPost)
plot( Omega, pOmegaPost, type="l", lwd=3, col="plum1", 
     yaxt='n', xlab=expression(omega), ylab=expression(paste("p(",omega,"|D)"))) 
