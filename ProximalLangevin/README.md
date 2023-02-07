


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

### ${\color{violet} \text{Pseudo-code}}$
* Define $\bar D$ a $T \times T$ invertible matrix, extension of $D_2$. 
* Denote by ${\bf A}$ the $(2T) \times (2T)$ block diagonal matrix with blocks $\bar D$ and $\lambda_R/\lambda_0  \, I_T$.  
* Initialize the chain in the image space: ${\bf Z} = {\bf A} ({\bf R}, {\bf O})$.
* Store the value of the log-density at the initial point. Denote by $\tilde \pi$ the image of $\pi$ by the linear transformation ${\bf A}$.
* Repeat
  * Sample a candidate in the image space 
    * Compute *Grad*, the gradient of the smooth component of $\log \tilde \pi$, evaluated at the current point $Z$.
    * Compute a gradient step: *GradStep* $= Z +\gamma$ *Grad*.
    * Modify components $3$ to $2T$: compute the proximal associated to the function $\lambda_R \gamma |\cdot|$, and applied to components $3$ to $2T$ of *GradStep*.
    * Add i.i.d. random variables with variance $2 \gamma$ on each of the $2T$ components. 
    * Denote by *Proposal* this Gaussian proposal. 
  * Test if *Proposal* is in the support of the distribution $\tilde \pi$.
  * If it is not: reject *Proposal* and set the new value of chain equal to $Z$.
  * It it is: implement an accept-rejection step
    * starting from *Proposal*, compute the proximal-gradient step as above (replacing $Z$ with *Proposal* in the method above), and obtain *DriftRev*
    * evaluate a Gaussian density with covariance matrix $2 \gamma I_{2T}$ at the point *DriftRev*
    * evaluate a Gaussian density with covariance matrix $2 \gamma I_{2T}$ at the point *Proposal*
    * evaluate $\log \tilde \pi$ at the point *DriftRev*
    * evaluae $\log \tilde \pi$ at the point *Proposal*
    * compute the acceptance-rejection ratio $\alpha$.
    * define the new value of the chain: with probability $\alpha$, set $Z = Proposal$; and otherwise do not change $Z$.
  * Store: the current value of the chain and the current value of $\log \pi$.
  * Update $\gamma$: during the burn-in phase, regularly update $\gamma$ in order to reach a targeted empirical acceptance-rejection rate.
* Discard the sample of the burn-in phase. With the other ones: 
  * move the chain to the original space by applying ${\bf A}^{-1}$.
  * for each of the $(2T)$ components of the chain, compute the quantiles and the empirical expectation. 



