# StochasticGrowthExp-OUProcess

Consider a stochastic growth process in discrete time for a quantity $B_t$ defined by $$B_{t+1}=(1+r_t) B_t$$. Assume for simplicity $B_0=1$.
The growth rate $r_t$ is assumed to be proportional to the exponential of a stationary Ornstein-Uhlenbeck process
$$r_t = \rho e^{Z_t - \frac12 var(Z_t)}$$ where the OU process is given by $$dZ_t = - \gamma Z_t dt + \sigma dW_t$$.
Recall that for a stationary process $Z = N(0,\frac{\sigma^2}{2\gamma})$.

This process describes for example the growth of a bank account accruing interest each period at an interest rate which follows an exp-OU process. This is the interest rates process assumed in the Black-Karasinski model. 

Let us study the expectation $M_t = \mathbb{E}[B_t]$.
