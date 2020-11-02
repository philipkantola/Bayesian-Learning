data = read.table("rainfall.dat")

# Setup 
my0 <- 1
#Increased value of tao0 results in bigger value of w hence we myn will go more towards the data mean
tao0squared <- 10
v0 <- 1
sigmaSquared <- 1  #init value for sigma, gets updated in Gibbs sampling
sigma0Squared <- 1  
nDraws <- 500 # Number of draws
dataMean = mean(data[,1])
n = length(data[,1])
vn = v0 + n

# Gibbs sampling
gibbsDraws <- matrix(0,nDraws,2)
theta2 <- 0
for (i in 1:nDraws){
  
  ## from lecture 2
  taonSquared = 1/(n/sigmaSquared + 1/tao0squared)
  w = (n/sigmaSquared)/(n/sigmaSquared + 1/tao0squared)
  myn = w*dataMean + (1-w)*my0
  ##
  
  ## from lecture 7
  
  #Compute Full conditional posterior (Normal model with conditionally conjugate prior)
  my <- rnorm(1, mean = myn, sd = sqrt(taonSquared))
  gibbsDraws[i,1] <- my
  
  sigmaSquaredVariance = ((v0*sigma0Squared + sum((data[,1]-my)^2))/n+v0)
  #Make it inverse chai squared
  sigmaSquared <- vn*sigmaSquaredVariance/(rchisq(1, vn))
  gibbsDraws[i,2] <- sigmaSquared
  ## 
  }

# Plot the values from the conditional posteriors of my|sigma and sigma|my
plot(gibbsDraws[,2], type="lines", xlab = "draw", ylab = "sigma", main = "Sigma conditioned on my")

plot(gibbsDraws[,1], type = "lines", xlab = "draw", ylab = "my", main = "My conditioned on sigma")
lines(dataMean, col = "red")


# (1b)


##Mattias Code
##########    BEGIN USER INPUT #################
# Data options
data(data)
rawData <- faithful
x <- as.matrix(data[,1])

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(10,nComp) # Prior mean of mu
tau2Prior <- rep(5,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(2,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green")
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))


for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

#########################   Enf of mattias code   ##############################################

#Extract what we are interested in
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)





# (1c)

my_posterior_mean = mean(gibbsDraws[,1])
sigma_posterior_mean = mean(gibbsDraws[,2])

# plot density of data
hist(x, breaks=20, freq=FALSE, ylim=c(0,0.025) ,main="Graphical comparison of data vs densities", xlab="Draw")
# plot posterior density from a)
lines(xGrid, dnorm(xGrid, mean=my_posterior_mean, sd = sqrt(sigma_posterior_mean)), type="lines", col="purple")
# plot mixture model density from b)
lines(xGrid, mixDensMean, type="lines", col = "red")
legend("topright", box.lty = 1, legend = c("Data density","Mixture density","Normal density"), col=c("black","red","purple"), lwd = 2)


#Assignment 2

#2a)


bids_DATA = read.table('eBayNumberOfBidderData.dat', header = TRUE)
#data = bids_DATA[,-2]

#Create the glm model. Gives us valus on the parameters which maximizes the likelihood for the given data. 
#Column 2 is taken away since glm model adds its own intercept

glm_model = glm(nBids ~ ., data = bids_DATA[,-2], family = poisson(link="log"))
Beta_covariates = glm_model$coefficients
summary(glm_model)
Beta_covariates


#2b
library(mvtnorm)


#covariance_matrix = c(1:8)
Y = as.vector(bids_DATA[,1])
X = as.matrix(bids_DATA[,-1])
col = ncol(X)

# parameters for the prior distribution of beta
my = as.vector(rep(0, col))
sigma_prior = 100*solve(t(X)%*%X)

# function which returns an expression proportional to log beta posterior, this can be optimized for 
# values of beta mode and hessian J (observed hessian evaluated at posterior mode)
#Log likelikihood is given by the log produkt of the poisson densituy distribution
log.posterior = function(betas, x,Y,my,sigma_prior){
  
  # logproduct of poission pdf
  log_likelihood = sum(-log(factorial(Y))+(betas%*%t(X))*Y - exp(betas%*%t(X)))
  
  # log of multivariate pdf => dmvnorm
  log_prior = dmvnorm(x=betas,mean=my,sigma=sigma_prior, log=TRUE)
  
  if(abs(log_likelihood == Inf)){
    log_likelihood = -5000
  }
  
  return(log_likelihood + log_prior)
}

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,col)); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode

#Find beta and hessian values in mode 
OptParams<-optim(initVal,log.posterior,gr=NULL,X,Y,my,sigma_prior,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

#Beta values
betas_posterior = OptParams$par
betas_posterior

# Takes negative so that the posterior can be approx. as normal
hessian_posterior = -solve(OptParams$hessian)
hessian_posterior



#2c
# Metropolis algorithm Le 8


#Draw samples from Betas posterior distribution. 
set.seed(12345)

#metopolis algorithm, simulates from posterior of beta. 
# args:  c = arbitrary tuning parameter, ... is a wildcard argument and represents the arguments of the posterior_density
# higher value of c gives bigger variance / bigger steps
# n = number of iterations (return n values for betas)
# "posterior density, ..." => unkown paramters will come in
metropolis = function(n, c, hessian, betas,  posterior_density, ...){

  # this step depends on previous position. Previous position becomes this turns mean. 
  proposal_draws_previous = betas;
  acceptedDraws = matrix(0, ncol=9,nrow=n)
  
  set.seed(12345)
  for(i in 1:n){
      # draws (theta_p) from the proposal distribution ~ N(theta_p-1, c*hessian)
      proposal_draws = rmvnorm(1, proposal_draws_previous, c*hessian)
      # create a ratio depending on if this draw is better than previous, take exp to remove logarithm (logposterior)
      # posterior_density = log.posterior => exp of the division => logA -logB 
      acceptance_ratio = min(1,exp(posterior_density(proposal_draws, ...)-posterior_density(proposal_draws_previous, ...)))
      # draw a random uniformed variable to compare wiht acceptance ratio
      random_acceptance = runif(1,0,1)
      # if acceptance ratio is bigger than random variable than we move to the new position, otherwise we stay
      if(acceptance_ratio >= random_acceptance){
        proposal_draws_previous = proposal_draws
        betas = proposal_draws
      }
      acceptedDraws[i,] = betas
      
  }
   acceptedDraws;
}

testDraws = metropolis(10000, c=1, hessian= hessian_posterior, betas = my, posterior_density = log.posterior, X, Y, my, sigma_prior)

# plot the result for each dimension of beta. 
for(i in 1:9){
  plot(testDraws[,i], type='s')
  a = c(rep(betas_posterior[i],nrow(testDraws)))
  lines(a, col='red')
}

#2d

# our x-values that we want to predict number of bids for
x_values <- c(1,1,1,1,0,0,0,1,0.5)

# mean of our beta values from the posterior
mean_betas = c()

# dont count in the first 1000 iterations in mean calculation (burn-in phase)
for(i in 1:9){
  mean_betas = append(mean_betas, mean(testDraws[1000:10000,i]))
}

#when given a Poisson regression model beta and an input vector x, 
#the predicted mean of the associated Poisson distribution is given by e^beta*x
#theta can be estimated by MLE, our as in our case, by metropolis.
poisson_regression <- exp(mean_betas%*% x_values)[1]

# simulate value for y considering our x-values
simulation <- rpois(50000,poisson_regression)
barplot(table(simulation))

# percentage of draws where y = 0
probZero = sum(simulation==0)/length(simulation)








