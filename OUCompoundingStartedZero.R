# exp-OU compounding process started at zero

n <- 10000 # no time steps
tau <- 0.01 # time step
gamma <- 0.1 # mean reversion
sigma <- 0.01 # volatility

rho <- 0.025
offset <- 0.025

dG <- sigma^2/(2*gamma)*(1-exp(-2*gamma*tau))
varZinfty <- sigma^2/(2*gamma) # stationary OU process variance
beta <- exp(-gamma*tau)
volinfty=sqrt(varZinfty)
volinfty

varZ <- function(t){
  
  x <- varZinfty*(1-exp(-2*gamma*tau*t))
  return(x)
}

# Generate random normals N(0,sqrt(dG))
e <- rnorm(n, mean=0, sd=sqrt(dG))

rate = numeric(n)
Z = numeric(n)

Z1 = 0  #rnorm(1,0,sqrt(varZ))

Z[1] = Z1
rate[1] = (rho*exp(Z1 - 0.5*varZ(1))-offset)/(1+offset)

for(i in 2:n) {
  Z[i]= beta*Z[i-1] + e[i]
  rate[i] = (rho*exp(Z[i] - 0.5*varZ(i))-offset)/(1+offset)
}

head(rate)
plot(rate, type="l", col="blue")
abline(h=0, col="red")

acf(rate,type="correlation", lag.max=100, main="ACF(r)")
