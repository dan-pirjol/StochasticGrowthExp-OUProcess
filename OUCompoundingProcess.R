# stationary exp-OU compounding process

n <- 1000 # no time steps
tau <- 0.01 # time step
gamma <- 0.1 # mean reversion
sigma <- 0.2 # volatility

rho <- 0.025
offset <- 0.025

dG <- sigma^2/(2*gamma)*(1-exp(-2*gamma*tau))
varZ <- sigma^2/(2*gamma) # stationary OU process variance
beta <- exp(-gamma*tau)

# Generate random normals N(0,sqrt(dG))
e <- rnorm(n, mean=0, sd=sqrt(dG))

rate = numeric(n)
Z = numeric(n)

Z1 = rnorm(1,0,sqrt(varZ))

Z[1] = Z1
rate[1] = (rho*exp(Z1 - 0.5*varZ)-offset)/(1+offset)

for(i in 2:n) {
  Z[i]= beta*Z[i-1] + e[i]
  rate[i] = (rho*exp(Z[i] - 0.5*varZ)-offset)/(1+offset)
}

head(rate)
plot(rate, type="l", col="blue",xlab="Time",main="Growth rate r(t)")
abline(h=0, col="red")

acf(rate,type="correlation", lag.max=1500, main="ACF(r)")
