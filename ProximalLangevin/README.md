


# Description of the matlab files

The codes below describe two Hastings-Metropolis samplers and two Gibbs samplers for the target density
  > $$ \log \pi(R_1, \cdots R_T, O_1, \cdots, O_T) =  - \sum_{t=1}^T \left( Z_t  \log(R_t \Phi_t + O_t) - (R_t \Phi_t + O_t)\right)+ \lambda_R \| {\bf D_2 R}\|_1 + \lambda_0 \|{\bf O}\|_1 $$
  > 
  > where ${\bf R} = (R_1, \cdots, R_T)$, ${\bf O} = (O_1, \cdots, O_T)$ and $\|\cdot\|_1$ denotes the L1-norm.   $D_2$ is a (T-2) x T Laplacian matrix.    
  > This density is positive when (i) $R_t \geq 0$ and (ii) $R_t\Phi_t+ O_t > 0$ when $Z_t >0$ and $R_t \Phi_t + O_t \geq 0$ when $Z_t \geq 0$.  The density is zero otherwise.
  
  
# ${\color{blue} \text{PGdual algorithm}}$

### ${\color{violet} \text{Input structures}}$

**data.Z**: T x 1.  The counts from time t=1 to time t=T, $(Z_1, \cdots, Z_T)$.    
**data.Zphi**: T x 1. The convolution of the counts and the serial function, $(\Phi_1, \cdots, \Phi_T)$.  

**MAP.lambdaR**: 1 x 1. Coefficient for the penalty term on the Laplacian of the reproduction numbers  
**MAP.lambdaO**: 1 x 1. Coefficient of the penalty term on the outliers.  
**MAP.method_augmentation**: a string, either 'orthogonal', or 'invert'. Defines how to extend the Laplacian matrix.  *The default value is 'orthogonal'*

**MCMC.chain_length**: 1 x 1. The total length of the Markov chain, including the burn-in samples.  *The default value is 1e7*.  
**MCMC.chain_burnin**: 1 x 1. The length of the burn-in phase.  *The default value is 50% of the chain_length*.  
**MCMC.initial_point**: (2T) x 1. The initial value of the bivariate chain $({\bf R,O})$.  *The default value is the constant vector 1 for R, and the null vector for O*  
**MCMC.target_ratio**: 1 x 1. A number in (0,1). For the adaptation of the step size, provide a target value for the acceptance-rejection ratio.  
**MCMC.gamma_init**: 1 x 1. A positive number, provides the initial value of the step size.  *The default value is 1e-7 when T=70*

**param.vectQ**: a row vector, that collects the order of the quantiles   of the distribution of each component of the bivariate chain, to be estimated from the sampler (by discarding the samples of the burn in phase).  
**param.displayglobal**: a binary variable set to '1' or '0' for the display (or not) of the graphical controls during the run of the algorithm.  *The default value is '0'*  
**param.frequency**: an interger, which defines the number of iterations between two updates of the step size.  *The default value is 1e4*

### ${\color{violet} \text{Output structures}}$
**generic.empirical_mean**: 2 x T. The empirical expectation of the bivariate chain $({\bf R,O})$, computed by discarding the samples of the burn-in phase.  
**generic.R_quantiles**: q x T. The q quantiles for each of the T components of ${\bf R}$,  computed by discarding the samples of the burn-in phase.  
**generic.O_quantiles**: q x T. The q quantiles for each of the T components of ${\bf O}$,  computed by discarding the samples of the burn-in phase.  
**additional.gamma**: collects the successive values of the step size, adapted during the burn-in phase (and no more adapted, after the burn-in phase).  
**additional.logPi**: collects the successive values of the log-density along iterations.   

### ${\color{violet} \text{Example}}$



