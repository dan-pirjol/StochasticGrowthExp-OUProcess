# StochasticGrowthExp-OUProcess

Consider a stochastic growth process in discrete time for a quantity $B_i$ defined by $$B_{i+1}=(1+r_i) B_i$$
Assume for simplicity $B_0=1$.

Take the growth rate $r_i$ proportional to the exponential of a stationary Ornstein-Uhlenbeck (exp-OU) process
$$r_i = \rho e^{Z_{t_i} - \frac12 var(Z_{t_i})}\qquad (1)$$ where the OU process $Z_t$ is sampled on uniformly spaced
times $t_i = i\tau$. The process $Z_t$ 
satisfies the SDE $$dZ_t = - \gamma Z_t dt + \sigma dW_t$$ and mean-reverts to zero.
Recall that for a stationary OU process $Z = N(0,\frac{\sigma^2}{2\gamma})$. Thus the expectation of the growth rate is constant 
$\mathbb{E}[r_i]=\rho$.

This process describes for example the growth of a bank account accruing interest each period at an interest rate which follows an exp-OU process. This is the interest rates process assumed in the Black-Karasinski model. A similar model can be used to describe the growth of a positive quantity with correlated growth rates, for example a population undergoing growth in a random environment. More precisely the growth rates $r_i$ have Markovian dependence, which is appropriate for example for an environmental variable such as temperature, oxygen supply or food resources.

The process can be generalized such that the growth rates $r_i$ can be both positive and negative, for example as 
$$B_{i+1} = \frac{1}{1+\rho} (1 + \rho e^{Z_i - \frac12 var(Z_i) }) B_i$$ which corresponds to a growth rate
$$r_i = \frac{1}{1+\rho} (\rho e^{Z_i - \frac12 var(Z_i) } - \rho )$$ This has $\mathbb{E}[r_i]=0$. For simplicity we will consider in the following the original version (1).

Let us study the expectation $M_t = \mathbb{E}[B_t]$. If the growth rates $r_t$ were uncorrelated, the expectation $M_t$ would have an exponential growth $M_n = \Pi_{i=0}^{n-1}(1 + \mathbb{E}[r_i]) = (1+\rho)^n$. However the serial correlation among growth rates $r_i$ changes the growth pattern in an unexpected way. 
