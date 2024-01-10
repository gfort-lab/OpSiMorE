Given
- T observations $Z_1, \cdots, Z_T$ taking values in $\mathbb{Z}_{\geq 0}$
- T mean values $\Phi_1, \cdots, \Phi_T$ taking values in $\mathbb{R}_{\geq 0}$
- Two initial values $R_{-1}, R_0$ for the reproduction number  taking values in $\mathbb{R}_{>0}$

two distributions are defined.

**Model "no mixture"** $\pi$ 

A target distribution $\pi$ on $(\mathbb{R}_{>0} \times \mathbb{R})^{T}$. The log-density $\ln \pi$ is given by

> $$ 
> {\small \begin{split}
> (R_1, O_1, \cdots, R_T, O_T) & \mapsto \sum_{t=1}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
> & - \frac{\lambda_R}{4} \sum_{t=1}^T | R_t - 2 R_{t-1} + R_{t-2} | + T \ln \lambda_R   \\
> & - \lambda_O  \sum_{t=1}^T |O_t| + T \ln \lambda_O
\end{split} }
> $$

up to an additive constant. The support of the distribution is  

> $$
 {\small \text{for} \ t=1, \cdots,T, \qquad R_t > 0; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi$ depends on two positive parameters $\lambda_R$ and $\lambda_O$.

**Model "mixture"**  $\pi_m$

A target distribution $\pi_m$ on $(\mathbb{R}_{>0} \times \mathbb{R} \times \\{0,1\\} )^{T}$. The log-density $\ln \pi_m$ is given by

>  $$ 
 {\small \begin{split} (R_1, O_1, B_1, \cdots, R_T, O_T, B_T) & \mapsto \sum_{t=1}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
& - \frac{\lambda_R}{4} \sum_{t=1}^T | R_t - 2 R_{t-1} + R_{t-2} | + T \ln \lambda_R   \\
& - \lambda_{O,0}  \sum_{t=1}^T (1-B_t) |O_t| + (T-\sum_{t=1}^T B_t) \ln \lambda_{O,0} \\
& - \lambda_{O,1}  \sum_{t=1}^T B_t |O_t| + (\sum_{t=1}^T B_t) \ln \lambda_{O,1} \\
& + (T-\sum_{t=1}^T B_t) \ln (1-\omega) +  (\sum_{t=1}^T B_t) \ln \omega 
\end{split}}
> $$ 

up to an additive constant. The support of the distribution is

> $$
> {\small \text{for} \ t=1, \cdots,T, \qquad R_t > 0; \qquad \qquad B_t \in \\{0,1\\}; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi_m$ depends on three positive parameters $\lambda_R$, $\lambda_{O,0}$ and $\lambda_{O,1}$; and a weight $\omega \in (0,1)$.


## ${\color{blue} \text{GibbsPGdual\\_nomixture}}$

This MATLAB code runs a Metropolis-within-Gibbs sampler with target distribution $\pi$ and returns a Monte Carlo approximation for each expectation. The proposal mechanism depends on design parameters: they are adapted during the burnin phase in order to target a given mean acceptance ratio. 

> $$ {\small
>   I_R := \frac{1}{4} \int \sum_{t=1}^T |r_t - 2 r_{t-1} + r_{t-2}| \ \  \mathrm{d} \pi(r_1,o_1, \cdots, r_T, o_T) \qquad \qquad   I_O := \int \sum_{t=1}^T |o_t| \ \  \mathrm{d} \pi(r_1,o_1, \cdots, r_T, o_T) 
> } $$

### ${\color{violet} \text{Input structures}}$
A structure _data_ with fields
- _Z_ : Tx1, the sequence $Z_1, \cdots, Z_T$
- _Phi_ : Tx1, the sequence $\Phi_1, \cdots, \Phi_T$
- Rinit : 2x1, the initial values $R_{-1}$ and $R_0$
- LambdaR : 1x1, the value of $\lambda_R$
- LambdaO : 1x1, the value of $\lambda_O$ 

A structure _MCMC_ with fields
- _chain\_length_ : 1x1, number of MCMC iterations; default value 1e7
-  _chain\_burnin_ : 1x1, length of the burnin period; default value 0.5*1e7
-  _initial\_pointR_ : Tx1, initial value of the R chain
-  _initial\_pointO_ : Tx1, initial value of the O chain
-  _GammaO_: 1x1, initial value of the step size when proposing a candidate for O; default value 1e3
-  _GammaTildeR_ : 1x1, initial value of the step size when proposing a candidate for the second derivative of R; default value 1e-12
-  _target\_ratioAR_ : 1x1, targeted acceptance ratio during the adaptation phase; default value 0.25
- _Qvec_ : vector of order of quantiles; default value (0.025 0.05 0.5 0.95 0.975)

  
### ${\color{violet} \text{Output structures}}$
A structure _output_ with fields
- _StatR_ : the Monte Carlo approximation of $I_R$
- _StatO_ : the Monte Carlo approximation of $I_O$
- _GammaTildeR_ : step size for the proposal mechanism when sampling the second derivative of the R_t variables
- _GammaO_ : step size for the proposal mechanism when sampling the O_t variables
- _empirical_meanR_ : Tx1, a Monte Carlo approximation of the expectation of $R_1, \cdots, R_T$ under the distribution $\pi$
- _empirical_meanO_ : Tx1, a Monte Carlo approximation of the expectation of $O_1, \cdots, O_T$ under the distribution $\pi$
- _quantilesR_ : length(Qvec) x T, the quantiles of $R_1, \cdots, R_T$ under the distribution $\pi$
- _quantilesO_ : length(Qvec) x T, the quantiles of $O_1, \cdots, O_T$ under the distribution $\pi$
- _lastsampleR_ : Tx1, the last MCMC sample R
- _lastsampleO_ : Tx1, the last MCMC sample O
- _LogPi_ : 1xchain\_length, the values of $\log \pi$ along the MCMC iterations 
