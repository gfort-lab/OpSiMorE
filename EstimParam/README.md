Given
- T observations $Z_1, \cdots, Z_T$ taking values in $\mathbb{Z}_{\geq 0}$
- T mean values $\Phi_1, \cdots, \Phi_T$ taking values in $\mathbb{R}_{\geq 0}$
- Two initial values $R_1, R_2$ for the reproduction number  taking values in $\mathbb{R}_{>0}$

two distributions are defined.

**Model "no mixture"** $\pi$ 

A target distribution $\pi$ on $(\mathbb{R}_{>0} \times \mathbb{R})^{T-2}$. The log-density $\ln \pi$ is given by

> $$ 
> {\tiny \begin{split}
> (R_3, O_3, \cdots, R_T, O_T) & \mapsto \sum_{t=3}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
> & - \frac{\lambda_R}{4} \sum_{t=3}^T | R_t - 2 R_{t-1} + R_{t-2} | + (T-2) \ln \lambda_R   \\
> & - \lambda_O  \sum_{t=3}^T |O_t| + (T-2) \ln \lambda_O
\end{split} }
> $$

up to an additive constant. The support of the distribution is  

> $$
 {\tiny \text{for} \ t=3, \cdots,T, \qquad R_t > 0; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi$ depends on two positive parameters $\lambda_R$ and $\lambda_O$.

**Model "mixture"**  $\pi_m$

A target distribution $\pi_m$ on $(\mathbb{R}_{>0} \times \mathbb{R} \times \\{0,1\\} )^{T-2}$. The log-density $\ln \pi_m$ is given by

>  $$ 
 {\tiny \begin{split} (R_3, O_3, B_3, \cdots, R_T, O_T, B_T) & \mapsto \sum_{t=3}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
& - \frac{\lambda_R}{4} \sum_{t=3}^T | R_t - 2 R_{t-1} + R_{t-2} | + (T-2) \ln \lambda_R   \\
& - \lambda_{O,0}  \sum_{t=3}^T (1-B_t) |O_t| + (T-2-\sum_{t=3}^T B_t) \ln \lambda_{O,0} \\
& - \lambda_{O,1}  \sum_{t=3}^T B_t |O_t| + (\sum_{t=3}^T B_t) \ln \lambda_{O,1} \\
& + (T-2-\sum_{t=3}^T B_t) \ln (1-\omega) +  (\sum_{t=3}^T B_t) \ln \omega 
\end{split}}
> $$ 

up to an additive constant. The support of the distribution is

> $$
> {\tiny \text{for} \ t=3, \cdots,T, \qquad R_t > 0; \qquad \qquad B_t \in \\{0,1\\}; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi_m$ depends on three positive parameters $\lambda_R$, $\lambda_{O,0}$ and $\lambda_{O,1}$; and a weight $\omega \in (0,1)$.


## ${\color{blue} \text{GibbsPGdual\\_nomixture}}$

This MATLAB code runs a Metropolis-within-Gibbs sampler with target distribution $\pi$ and returns a Monte Carlo approximation for each expectation 

> $$ {\tiny
>   I_R := \frac{1}{4} \int \sum_{t=3}^T |r_t - 2 r_{t-1} + r_{t-2}| \ \  \mathrm{d} \pi(r_3,o_3, \cdots, r_T, o_T) \qquad \qquad   I_O := \int \sum_{t=3}^T |o_t| \ \  \mathrm{d} \pi(r_3,o_3, \cdots, r_T, o_T) 
> } $$

### ${\color{violet} \text{Input structures}}$
A structure _data_ with fields
- _Z_ : (T-2)x1, the sequence $Z_3, \cdots, Z_T$
- _Phi_ : (T-2)x1, the sequence $\Phi_3, \cdots, \Phi_T$
- Rinit : 2x1, the initial values $R_1$ and $R_2$

A structure _MCMC_ with fields
- _NbrMC_ : 1x1, number of MCMC iterations
-  _burnin_ : 1x1, length of the burnin period
  
### ${\color{violet} \text{Output structures}}$
A structure _output_ with fields
- _StatR_ : the Monte Carlo approximation of $\mathcal{I}_R$
- _StatO_ : the Monte Carlo approximation of $\mathcal{I}_O$
- _gammaR_ : step size for the proposal mechanism when sampling the R_t variables
- _gammaO_ : step size for the proposal mechanism when sampling the O_t variables
