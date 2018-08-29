

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())

require(rjags)
source('DBDA2E-utilities.R')
Niter = 200000
Nburn = 500
# Zooms the plot label texts
cx =1.9 

#################################################################################

# Creates the JAGS sample as mcmc.list from the specified model and data
mcmcf = function(data, params = c("theta","alpha","beta"),model) {
    
    jags = jags.model( model,
                            data = data,
                            n.chains = 4,
                            n.adapt = 500)
    ## Burn-In
    update(jags, n.iter=Nburn)

    mcmclist = coda.samples (jags,
                  variable.names = params,
                  n.iter = Niter,
                  thin=1)
}
#################################################################################
# Creates a data set of 10 0-bags and 10 1-bags

emptyBags = function() {
    dataList = list(
        Nbags = 21 
    )
}

uniformBags = function() {
# All bags are completely uniform

    y = c(rep(1,200),rep(0,200),1)
    Nbags = 20
    bags = c()
    for(bIdx in 1:Nbags) {
        bags = c(bags, rep(bIdx,20))
    }
    Nbags = Nbags + 1
    bags = c(bags,Nbags)

    dataList = list(
        y = y,
        bags = bags,
        Ntotal = length(y),
        Nbags = Nbags 
    ) 
}

# Creates 20 mixed bags, skewed towards 0-side
mixedBags = function() {

    Nbags = 20

    y = c()
    for(i in 1:(Nbags/2)) {
        add = c(rep(1,i), rep(0,Nbags-i)) 
        y = c(y, add, add)
    }
    y = c(y,1)

    bags = c()
    for(bIdx in 1:Nbags) {
        bags = c(bags, rep(bIdx,20))
    }
    Nbags = Nbags + 1
    bags = c(bags,Nbags)

    dataList = list(
        y = y,
        bags = bags,
        Ntotal = length(y),
        Nbags = Nbags 
    )
}

#################################################################################

# Plots alpha, beta, theta1, theta11, theta21 
ploth = function(mcmcl, main = "Prior") {

    ylim = c(0,8)

    plot(density(x=log(as.matrix(mcmcl[,'alpha']))),
         xlim=c(-8,2),
         ylim=c(0,1),
         xlab = expression("log("*alpha*")"),
         ylab = expression("p("*alpha*")"),
         main = main,
         cex.lab=cx,
         col="plum1", lwd=3)

    plot(density(x=as.matrix(mcmcl[,'beta']),
                 from=0.02,
                 to=0.98),
         ylim=ylim,
         xlab = expression(beta),
         ylab = expression("p("*beta*")"),
         main = "",
         cex.lab=cx,
         col="plum1", lwd=3)


    ## plot(density(x=as.matrix(mcmcl[,'theta[1]']),
    ##              from=0.02,
    ##              to=0.98),
    ##      ylim=ylim,
    ##      xlab = expression(theta[1]),
    ##      ylab = expression("p("*theta[1]*")"),
    ##      main = "",
    ##      cex.lab=cx,
    ##      col="skyblue", lwd=3)

    plot(density(x=as.matrix(mcmcl[,'theta[11]']),
                 from=0.02,
                 to=0.98),
         ylim=ylim,
         xlab = expression(theta[11]),
         ylab = expression("p("*theta[11]*")"),
         main = "",
         cex.lab=cx,
         col="skyblue", lwd=3)
    plot(density(x=as.matrix(mcmcl[,'theta[21]']),
                 from=0.02,
                 to=0.98),
         ylim=ylim,
         xlab = expression(theta[21]),
         ylab = expression("p("*theta[21]*")"),
         main = "",
         cex.lab=cx,
         col="skyblue",  lwd=3)
}

plotf = function(mcmcl, main="Prior") {

    ylim = c(0,8)

    plot(density(x=as.matrix(mcmcl[,'theta[1]']),
                 from=0.02,
                 to=0.98),
         xlab = expression(theta[1]),
         ylab = expression("p("*theta[1]*")"),
         main = main,
         cex.lab=cx,
          ylim=ylim,
         col="skyblue", lwd=3)

    plot(density(x=as.matrix(mcmcl[,'theta[11]']),
                 from=0.02,
                 to=0.98),
         xlab = expression(theta[11]),
         ylab = expression("p("*theta[11]*")"),
         main = "",
         ylim=ylim,
         cex.lab=cx,
         col="skyblue", lwd=3)
    plot(density(x=as.matrix(mcmcl[,'theta[21]']),
                 from=0.02,
                 to=0.98),
         xlab = expression(theta[21]),
         ylab = expression("p("*theta[21]*")"),
         main = "",
         cex.lab=cx,
         ylim=ylim,
         col="skyblue", lwd=3)
}

#################################################################################
# Hierarchical Model 
# Trick Jags into generating priors by providing no data
hierPri = mcmcf(data = emptyBags(),
               model = "MultiHierJAGSPriors.bug")

# Hierarchical Posterior w/ Uniform Data
hierPosUnif = mcmcf(data = uniformBags(),
                    model = "MultiHierJAGS.bug")

# Hierarchical Posterior w/ Mixed Data
hierPosMixed = mcmcf(data = mixedBags(),
                    model = "MultiHierJAGS.bug")

# Flat Model
flatPri = mcmcf(data = emptyBags(),
                model = "MultiJAGSPriors.bug")

flatPosUnif = mcmcf(data = uniformBags(),
                    model = "MultiJAGS.bug")

flatPosMixed = mcmcf(data = mixedBags(),
                    model = "MultiJAGS.bug")
#################################################################################

openGraph(width=10,height=12)
## layout(matrix( c(1:20,
##                  rbind(c(0,0,21,22,23),
##                        c(0,0,24,25,26)),
##                  ncol=6)))
## m = rbind(c( 1, 6,11),
##           c( 2, 7,12),
##           c( 3, 8,13),
##           c( 4, 9,14),
##           c( 5,10,15))
m = rbind(c(1,5,9),
          c(2,6,10),
          c(3,7,11),
          c(4,8,12))
layout(m)

mainUnif = "Posterior, uniform bags."
mainMix = "Posterior, mixed Bags."
par(mar=c(4,5,4 ,5))
ploth(hierPri)
ploth(hierPosMixed, main=mainMix)
ploth(hierPosUnif, main=mainUnif)

saveGraph(file="Hierarchical20Bags",type="pdf")
## # Open a new plot window
X11()

openGraph(width=10,height=10)
m = rbind(c(1,4,7),
          c(2,5,8),
          c(3,6,9))
par(mar=c(4,5,4,5))
layout(m)
plotf(flatPri)
plotf(flatPosMixed, main=mainMix)
plotf(flatPosUnif, main=mainUnif)
saveGraph(file="Flat20Bags",type="pdf")
