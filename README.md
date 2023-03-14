# Stochastic growth process with exp-OU growth rates

Numerical simulations for the paper [Growth rate of a stochastic growth process driven by an exponential Ornstein-Uhlenbeck process (2022)](https://aip.scitation.org/doi/10.1063/5.0065342)

Consider a stochastic growth process in discrete time for a quantity $B_i$ with random growth rates $r_i$ $$B_{i+1}=(1+r_i) B_i$$
started at $B_0=1$.

Take the growth rate $r_i$ proportional to the exponential of a stationary Ornstein-Uhlenbeck (exp-OU) process
$$r_i = \rho e^{Z_{t_i} - \frac12 var(Z_{t_i})}\qquad (1)$$ where the OU process $Z_t$ is sampled on uniformly spaced
times $t_i = i\tau$ with time step $\tau$. The process $Z_t$ 
satisfies the SDE $$dZ_t = - \gamma Z_t dt + \sigma dW_t$$ and mean-reverts to zero.
Recall that for a stationary OU process $Z = N(0,\frac{\sigma^2}{2\gamma})$. Thus the expectation of the growth rate is constant 
$\mathbb{E}[r_i]=\rho$.

This process describes for example the growth of a bank account accruing interest each period at an interest rate which follows an exp-OU process. This is the interest rates process assumed in the [Black-Karasinski model](https://en.wikipedia.org/wiki/Black%E2%80%93Karasinski_model). 
A similar model can be used to describe the growth of a positive quantity with correlated growth rates, for example a population undergoing growth in a random environment. More precisely the growth rates $r_i$ have Markovian dependence, which is appropriate for example for an environmental variable such as temperature, oxygen supply or food resources.

The process can be generalized such that the growth rates $r_i$ can be both positive and negative, for example as 
$$B_{i+1} = \frac{1}{1+\rho} (1 + \rho e^{Z_i - \frac12 var(Z_i) }) B_i\qquad (2)$$ which corresponds to a growth rate
$$r_i = \frac{1}{1+\rho} (\rho e^{Z_i - \frac12 var(Z_i) } - \rho )$$ This has $\mathbb{E}[r_i]=0$. For simplicity we will consider in the following the original version (1).

Let us study the expectation $M_t = \mathbb{E}[B_t]$. If the growth rates $r_t$ were uncorrelated, the expectation $M_t$ would have an exponential growth $M_n = \Pi_{i=0}^{n-1}(1 + \mathbb{E}[r_i]) = (1+\rho)^n$. The exponential growth rate is $\lambda_n := \frac{1}{n} \log M_t=\log(1+\rho)$.

We will show that the serial correlation among growth rates $r_i$ changes the growth pattern in an unexpected way. 

## **Simulation**

```
# stationary exp-OU compounding process

n <- 1000 # no time steps
tau <- 0.01 # time step
gamma <- 0.1 # mean reversion
sigma <- 0.2 # volatility

rho <- 0.025

dG <- sigma^2/(2*gamma)*(1-exp(-2*gamma*tau))
varZ <- sigma^2/(2*gamma) # stationary OU process variance
beta <- exp(-gamma*tau)

# Generate random normals N(0,sqrt(dG))
e <- rnorm(n, mean=0, sd=sqrt(dG))

rate = numeric(n)
Z = numeric(n)

Z1 = rnorm(1,0,sqrt(varZ))

Z[1] = Z1
rate[1] = (rho*exp(Z1 - 0.5*varZ)-rho)/(1+rho)

for(i in 2:n) {
  Z[i]= beta*Z[i-1] + e[i]
  rate[i] = (rho*exp(Z[i] - 0.5*varZ)-rho)/(1+rho)
}

plot(rate, type="l", col="blue",xlab="Time",main="Growth rate r(t)")
abline(h=0, col="red")

acf(rate,type="correlation", lag.max=1500, main="ACF(r)")

```

The autocorrelation of the growth rates for $\gamma=0.1,\sigma=0.2$ shows long-range correlation

<img width="312" alt="test2" src="https://user-images.githubusercontent.com/60016102/225134807-94f3f63e-deb7-4c25-be32-d7296e619ee2.png">

The left plot below shows the expectation $M_n$ for $n=100$ and time step $\tau=0.01$ with $\gamma=0.1$ as $\sigma$ increases. The parameter $\rho=0.025$. The horizontal black line is at 1. The expectation increases with $\sigma$, first smoothly and then more erratically.

<img width="592" alt="gamma0p1" src="https://user-images.githubusercontent.com/60016102/225136677-364c4a73-ea9f-4199-8d6e-f429bb033aaa.png">

The growth rate $\lambda_n = \frac{1}{n} \log M_n$ is shown in the right plot. For small $\sigma$, it is close to $\log(1+\rho) \simeq \rho$,
but as $\sigma$ increases, it becomes larger. 
