library(rstan)


t = seq(1,200, 1)
my = 10
sigma_square = 2


set.seed(12345)

ARprocess = function(t,my,sigma_square,phi){
  time_points = c()
  for (i in 1:length(t)){
    if(i == 1){
      time_points[i] = my
      next
    }
    
    error_term = rnorm(1, mean = 0, sd = sqrt(sigma_square))
    
    time_points[i] = my + phi*(time_points[i-1] - my) + error_term
  }
  return(time_points)
}

time_points0.5 = ARprocess(t,my,sigma_square,0.5)
time_points0.99 = ARprocess(t,my,sigma_square,0.99)
time_pointsneg0.99 = ARprocess(t,my,sigma_square,-0.99)


plot(t, time_points0.5, type="l",main = "phi = 0.5", ylab = "Xt")
plot(t, time_points0.99, type="l",main = "phi = 0.99", ylab = "Xt")
plot(t, time_pointsneg0.99, type="l",main = "phi = -0.99", ylab = "Xt")



#1b)


time_serie = function(phi){
array = c()
  
for (i in 1:length(t)){
  if(i == 1){
    array[i] = my
    next
  }
  
  error_term = rnorm(1, mean = 0, sd = sqrt(sigma_square))
  array[i] = my + phi*(array[i-1] - my) + error_term
}
  return(array)
}


time_serie_0.95 = time_serie(0.95)
time_serie_0.3 = time_serie(0.3)

# aplha = my
# beta = phi
# sigma = sigma

alpha = 1
beta = 0.5
sigma = 2

stanModel = ' data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  for (n in 2:N)
    y[n] ~ normal(alpha + beta * (y[n-1]-alpha), sigma);
}'

burnin = 1000
niter = 3000

# for phi = 0.95
data = list(N=length(t), y=time_serie_0.95)
fit_0.95 = stan(model_code=stanModel,data=data,
           warmup=burnin,iter=niter,chains=2)
# Print the fitted model
print(fit_0.95,digits_summary=3)

# for phi = 0.3
data = list(N=length(t), y=time_serie_0.3)
fit_0.3 = stan(model_code=stanModel,data=data,
           warmup=burnin,iter=niter,chains=2)
# Print the fitted model
print(fit_0.3,digits_summary=3)


# Extract posterior samples
postDraws_0.95 <- extract(fit_0.95)
postDraws_0.3 <- extract(fit_0.3)

par(mfrow = c(1,1))
plot(postDraws_0.95$beta,type="l",ylab="Phi", main="Sample convergence for phi when phi = 0.95")
plot(postDraws_0.3$beta,type="l",ylab="Phi", main="Sample convergence for phi when phi = 0.3")

plot(postDraws_0.95$alpha,type="l",ylab="My", main="Sample convergence for my when phi = 0.95")
plot(postDraws_0.3$alpha,type="l",ylab="My", main="Sample convergence for my when phi = 0.3")


plot(postDraws_0.95$alpha, postDraws_0.95$beta,type="l",ylab="Phi",xlab = "My", main="Joint posterior for phi = 0.95")
plot(postDraws_0.3$alpha, postDraws_0.3$beta,type="l",ylab="Phi",xlab = "My", main="Joint posterior for phi = 0.3")


# Assignment 3

c_data = read.table('campy.dat', header = TRUE)



C_DataRStan<-
  list(N = nrow(c_data),
       c = c_data$c) 

# The poissons model
C_Model = '
data {
  int<lower=0> N;
  int c[N];
}



parameters {
  real alpha;
  real <lower=-1, upper= 1> beta;
  real<lower=0> sigma;
  real Xt[N];
}

model {

  // Prior
  for (n in 2:N){
    Xt[n] ~ normal(alpha + beta * (Xt[n-1]-alpha), sigma);
  }
  // Model/likelihood
  for(n in 1:N){
  
  c[n] ~ poisson(exp(Xt[n]));

  }                 
}'

fit_C_Model<-stan(model_code=C_Model,
               data=C_DataRStan,
               warmup=1000,
               iter=2000,
               chains=2)

print(fit_C_Model,digits=4)

# Plot some results
res<-extract(fit_C_Model)

x <- res$x

poster_means_vector <- c()

for (i in 1:length(x[1,])) {
  poster_means_vector <- c(poster_means_vector, mean(exp(x[,i])))
}


cred_lower = fit_C_Model@.MISC[["summary"]][["c_quan"]][,1,1][-c(1,2,3,144)]
cred_upper = fit_C_Model@.MISC[["summary"]][["c_quan"]][,9,1][-c(1,2,3,144)]

N = nrow(c_data)
xGrid = seq(1,140,1)
plot(xGrid, c_data$c, col = "red", main = "Possion model & Data", ylab = "nr of Cases", xlab = "Draw")
lines(exp(poster_means_vector[1:N]), ylim = c(0,60), col = "black", type = 'l')
lines(exp(cred_lower), col = "blue")
lines(exp(cred_upper), col = "green")

legend('topleft',legend = c('Data', 'Mean','Lower','Upper'), col = c('red', 'black', 'blue', 'green'), lwd=2)





##2d)



# The poissons model
C_Model = '
data {
  int<lower=0> N;
  int c[N];
}



parameters {
  real alpha;
  real <lower=-1, upper= 1> beta;
  real<lower=0> sigma;
  real Xt[N];
}

model {

  // Prior for sigma
  sigma~inv_gamma(1,1)
  
  // Prior
  for (n in 2:N){
    Xt[n] ~ normal(alpha + beta * (Xt[n-1]-alpha), sigma);
  }
  // Model/likelihood
  for(n in 1:N){
  
  c[n] ~ poisson(exp(Xt[n]));

  }                 
}'

fit_C_Model<-stan(model_code=C_Model,
                  data=C_DataRStan,
                  warmup=1000,
                  iter=2000,
                  chains=2)

print(fit_C_Model,digits=4)

# Plot some results
res<-extract(fit_C_Model)


poster_means_vector = fit_C_Model@.MISC[["summary"]][["msd"]][-c(1,2,3)]

cred_lower = fit_C_Model@.MISC[["summary"]][["c_quan"]][,1,1][-c(1,2,3,144)]
cred_upper = fit_C_Model@.MISC[["summary"]][["c_quan"]][,9,1][-c(1,2,3,144)]

N = nrow(c_data)
xGrid = seq(1,140,1)
plot(xGrid, c_data$c, col = "red", main = "Possion model & Data", ylab = "nr of Cases", xlab = "Draw")
lines(exp(poster_means_vector[1:N]), ylim = c(0,60), col = "black", type = 'l')
lines(exp(cred_lower), col = "blue")
lines(exp(cred_upper), col = "green")

legend('topleft',legend = c('Data', 'Mean','Lower','Upper'), col = c('red', 'black', 'blue', 'green'), lwd=2)





