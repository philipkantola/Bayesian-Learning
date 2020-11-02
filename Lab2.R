# for multivariate normal functions
library(mvtnorm)

V_0 = 4
sigma_0 = 1
my_0 = c(-10,100,-100)
set.seed(12345)
omega_0 = 0.01*diag(3)
omega_0_Inverse = solve(omega_0)

set.seed(12345)
#Draw 10 draws from out given chi sqared distribution
XDraw = rchisq(10,V_0)
#Transform our 10 draws to the scaled inverse chi square
deviationDraw = V_0*sigma_0/XDraw 

# For each of the draws above we create the joint prior obtaining the distribution for Beta given sigma
# each draw gives us a set of values for the betas
#rmvnorm because of that my_0 is a vector => multivariate
# sigma is covarance matrix
joint_prior = list()
for(i in 1:10){
  set.seed(12345)
  joint_prior[[i]] =  rmvnorm(1, mean = my_0, sigma = deviationDraw[i]*omega_0_Inverse)
}

# Function for computing the regression function for different values of x (time)
function_temperature = function(betas,x){
  betas[1] + betas[2]*x + betas[3]*x^2
}

#Plot the temperature regression curve for every set of values for betas
time = seq(0,1,1/365)
plot(function_temperature(joint_prior[[1]][1,],time), type = 'l', ylim = c(-25,35), xlab = "Days since start of 2018", ylab = "Temperature")
for(i in 2:10){
  lines(function_temperature(joint_prior[[i]][1,],time), type = 'l')
  
}


# (b)

data = read.table("TempLinkoping.txt", header = TRUE)
X = data$time
data$xsquare = X^2  
data$intercept = 1
Y = data$temp
X_matrix = matrix(0, nrow = 365, ncol = 3)
X_matrix[,1] = data$intercept
X_matrix[,2] = data$time
X_matrix[,3] = data$xsquare

beta_hat = solve((t(X_matrix)%*%X_matrix))%*%t(X_matrix)%*%Y

# get value for my_n
my_n = solve(((t(X_matrix)%*%X_matrix)) + omega_0)%*%(((t(X_matrix)%*%X_matrix%*%beta_hat))+omega_0%*%my_0)

omega_n = (t(X_matrix)%*%X_matrix) + omega_0
omega_n_inverse = solve(omega_n)

# degrees of freedom
Vn = V_0 + length(X)


sigma_n = ((V_0*sigma_0 + (t(Y)%*%Y + t(my_0)%*%omega_0%*%my_0)[1] - (t(my_n)%*%omega_n%*%my_n)[1]))/Vn

# number of draws
n = 1000

# draw values from the marginal posterior distribution of sigma
set.seed(12345)
#Draw n draws from our given chi sqared distribution
XPostDraw = rchisq(n,Vn)
#Transform our n draws to the scaled inverse chi square
deviationPostDraw = Vn*sigma_n/XPostDraw 

# For each of the draws above we create the joint marginal posterior obtaining the distribution for Beta given sigma
# each draw gives us a set of values for the betas
joint_posterior = matrix(0, nrow = n, ncol = 3)
beta1 = c()
beta2 = c()
beta3 = c()


set.seed(12345)
for(i in 1:n){
  joint_posterior[i,] =  rmvnorm(1, mean = my_n, sigma = omega_n_inverse*deviationPostDraw[i])
  beta1 = append(beta1,joint_posterior[i,1])
  beta2 = append(beta2,joint_posterior[i,2])
  beta3 = append(beta3,joint_posterior[i,3])
  temp = temperature_posterior(joint_posterior[i,], X_matrix)
}


# visualize the simulated parameters in histograms
hist(deviationPostDraw)
hist(beta1)
hist(beta2)
hist(beta3)

# Function f(time) for certain beta values on every values of x (time)
temperature_posterior = function(betas,x){
  x%*%betas
}

temperature_matrix = matrix(0, nrow = 365, ncol = 1000)
for(i in 1:n){
  temperature_matrix[,i] = temperature_posterior(joint_posterior[i,], X_matrix)
}

#Get median value of f(x) computed for each day
median_temperature_vec = apply(temperature_matrix,1, median) 

CI_temperature_vec = apply(X = temperature_matrix, MARGIN=1, FUN=function(x) quantile(x,c(0.05,0.95)))


#Plot posterior median, upper/lower credible interval for the Beta values and SMHI data
plot(median_temperature_vec, type = 'l', ylim = c(-15,25), main = "Posterior median of Betas", xlab = "Days since start of 2018", ylab = "Temperature", col = "blue")
lines(Y, col = "Black")
lines(CI_temperature_vec[1,], col = "green")
lines(CI_temperature_vec[2,], col = "red")
legend("topleft", c("SMHI data", "Posterior median", "Upper", "Lower"), col = c("black", "blue", "red", "green"), pch = 21:22, lty = 1:2)


# 1 c)
# function which takes in beta values and returns which x (which day) gives the maximum value
maximal_time = function(beta2,beta3){
  
   x =  -beta2/(2*beta3)
}

# collect all the max x-values for every set of betas
max_vector = c()
for(i in 1:n){
  # multiply with 366 to get nr of days instead of range 0,1
  max_vector = append(max_vector,maximal_time((joint_posterior[i,2]),(joint_posterior[i,3]))*365)
}

hist(max_vector, main="Distribution for mode of x", xlab="Number of days since start of year")



#Assignment 2 

women_data = read.table('WomenWork.dat', header = TRUE)
women_data
n = nrow(women_data)

tao = 10
#covariance_matrix = c(1:8)
Y = as.vector(women_data[,1])
X = as.matrix(women_data[,-1])
col = ncol(X)

# parameters for the prior distribution of beta
my = rep(0, col)
sigma_prior = diag(tao^2, nrow = col)

# function which returns an expression proportional to log beta posterior, this can be optimized for 
# values of beta mode and hessian J (observed hessian evaluated at posterior mode)
log.posterior = function(betas, x,Y,my,sigma_prior){
  
  # is simply log of the product of density function
  log_likelihood = sum((X%*%betas)*Y-log(1+exp(X%*%betas)))
  
  #log of deensity for nultivariate normal => dmvnorm(density multivariate)
  log_prior = dmvnorm(x=betas,mean=my,sigma=sigma_prior, log=TRUE)
  
  return(log_likelihood + log_prior)
  
}

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,col)); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode
OptParams<-optim(initVal,log.posterior,gr=NULL,X,Y,my,sigma_prior,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

betas_posterior = OptParams$par
# takes negative so that the posterior can be approx. as normal
# J = -second derivate evaluated in theta_hat
hessian_posterior = -OptParams$hessian

# take inverse for using it in the formula
hessian_posterior = solve(hessian_posterior)

#Draw samples from Betas posterior distribution. 
set.seed(12345)
posterior_distribution = rmvnorm(n = 1000, mean = betas_posterior, sigma = (hessian_posterior))

#Distribution gives for each draw a beta vector of length 8
#Take out vecotr of interest (no of small children)
n_child = posterior_distribution[,7]

#95% credible interval
quantiles = quantile(n_child, probs = seq(0,1,0.025))
quantiles
interval = c(quantiles[2], quantiles[40])
interval

#Comparision with glm model
glmModel <- glm(Work ~ ., data = women_data, family = binomial)
glmModel$coefficients

#2b
#Write a function that simulates from the predictive distribution of the response
#variable in a logistic regression. Use your normal approximation from 2(a).
#Use that function to simulate and plot the predictive distribution for the Work
#variable for a 40 year old woman, with two children (3 and 9 years old), 8 years
#of education, 10 years of experience. and a husband with an income of 10.

# expression for the response variable y (logistic regression)
work_distribution = function(x_data, betas){
  exp(t(x_data)%*%betas) / (1+exp(t(x_data)%*%betas))
}

# the properties of the type of woman we want the work distribution for 
x_data = c(1,10.000,8,10,1,40,1,1)

outcome = c()

# Simulating from teh predictive distribution (hence giving teh actual outcomes) for y with respect to our x values by simulating beta values from the beta posterior. 
set.seed(12345)
for (i in 1:1000){
  betas = rmvnorm(n = 1, mean = betas_posterior, sigma = (hessian_posterior))
  betas = as.vector(betas)
  prob = work_distribution(x_data, betas)
  # y is binary, therefore bernoulli distributed. Same thing as binomial with one draw. 
  outcome = append(outcome, rbinom(1,1, prob))
  
  
}
# histogram representing prob of woman work
hist(outcome, main="Predictive distribution for if woman works", xlab="Woman work or not", ylab="Acumulation of women who works out of 1000 draws")


# 2c

#Now, consider 10 women which all have the same features as the woman in 2(b).
#Rewrite your function and plot the predictive distribution for the number of
#women, out of these 10, that are working. [Hint: Which distribution can be
#described as a sum of Bernoulli random variables?]

# the distribution will be a binomial

# expression for the binomial dist of the quantity out of 10 woman of same type that works. Depends on our 
# posterior dist for the probability that one woman of this type works
work_binomial = function(x_data, betas){
  p = work_distribution(x_data, betas)
  rbinom(1,10,p)  
}

work_binomial_result = c()

# Simulating the distribution for quantity out of 10 woman works (binomial)
# using posterior dist for beta vectors and for work probability.(work distribution)
set.seed(12345)
for (i in 1:1000){
  betas = rmvnorm(n = 1, mean = betas_posterior, sigma = (hessian_posterior))
  betas = as.vector(betas)
  work_binomial_result = append(work_binomial_result,work_binomial(x_data,betas))
}

hist(work_binomial_result, main="Quantity out of 10 women who works", xlab="Number of women working out of 10", ylab="Acumulation of 1000 draws")

# one can check these quantiles to compare with prior histogram mode
#quantile(work_binomial_result, probs = seq(0,1,0.025))


