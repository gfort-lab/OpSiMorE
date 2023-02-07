


# Description of the matlab files and the data sets

The codes below describe two Hastings-Metropolis samplers and two Gibbs samplers for the target density
  > $$ \log \pi(R_1, \cdots R_T, O_1, \cdots, O_T) =  - \sum_{t=1}^T \left( Z_t  \log(R_t \Phi_t + O_t) - (R_t \Phi_t + O_t)\right)+ \lambda_R \| {\bf D_2 R}\|_1 + \lambda_0 \|{\bf O}\|_1 $$
  > 
  > where ${\bf R} = (R_1, \cdots, R_T)$, ${\bf O} = (O_1, \cdots, O_T)$ and $\|\cdot\|_1$ denotes the L1-norm.   $D_2$ is a (T-2) x T Laplacian matrix.    
  > This density is positive when (i) $R_t \geq 0$ and (ii) $R_t\Phi_t+ O_t > 0$ when $Z_t >0$ and $R_t \Phi_t + O_t \geq 0$ when $Z_t \geq 0$.  The density is zero otherwise.
  > 

The algorithm below designs Markov chain Monte Carlo samplers with target distribution $\pi$. They return, for each of the $2T$ components, (i) the empirical expectation, (ii) the empirical quantiles computed along the path of the Markov chain (after a burn in phase). 
  
The Covid19 data provided in the data sets are those made available at the Johns Hopkins University (https://coronavirus.jhu.edu)  
  
# ${\color{blue} \text{PGdual algorithm}}$

### ${\color{violet} \text{Input structures}}$

**data.Z**: T x 1.  The counts from time t=1 to time t=T, $(Z_1, \cdots, Z_T)$.    
**data.Zphi**: T x 1. The convolution of the counts and the serial function. The serial function $\phi_1, \cdots, \phi_T$ is equal to the probability density function of a Gamma distribution with shape parameter   and scale parameter , evaluated at $t=1, \cdots, t=T$. 

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
**output.empirical_mean**: 2 x T. The empirical expectation of the bivariate chain $({\bf R,O})$, computed by discarding the samples of the burn-in phase.  
**output.R_quantiles**: q x T. The q quantiles for each of the T components of ${\bf R}$,  computed by discarding the samples of the burn-in phase.  
**output.O_quantiles**: q x T. The q quantiles for each of the T components of ${\bf O}$,  computed by discarding the samples of the burn-in phase.  
**output.gamma**: collects the successive values of the step size, adapted during the burn-in phase (and no more adapted, after the burn-in phase).  
**output.logPi**: collects the successive values of the log-density along iterations.   

### ${\color{violet} \text{Example}}$
```
%% load data.Z and data.Zphi and MCMC.initial_point
% The initial value was obtained as follows: 
    % set:  MCMC.initial_point = [ones(T,1); zeros(T,1)]; 
    % run GibbsPGdual with MAP.method_augmentation = 'invert';
    % save the last value of the Markov chain (output.lastsample)
    % set: MCMC.initial_point = output.lastsample;
load FranceDataSet1.mat 

MAP.lambdaR = 3.5 * std(data.Z);
MAP.lambdaO = 0.05;
MAP.method_augmentation = 'orthogonal'; % could be 'invert'

param.T = length(data.Z);
param.displayglobal = 0; % could be '1' for the display of  graphical controls during the run
param.vecQ = [0.025 0.5 0.075];

MCMC.chain_length = 1.2e7;
MCMC.chain_burnin = ceil(0.5*MCMC.chain_length);
MCMC.gamma_init = 1e-7;
MCMC.target_ratio = 0.25;

outputPGdual = PGdual_github(data,MAP,MCMC,param);
```


# ${\color{blue} \text{PGdec algorithm}}$

### ${\color{violet} \text{Input structures}}$ 
The same as **PGdual**, except that   
**MAP.method_augmentation** is not required.  
**MCMC.covariance** is required. It describes the covariance matrix of the Gaussian proposal. Its value is either 'orthogonal', 'invert', or 'identity'.

### ${\color{violet} \text{Output structures}}$   
The same as **PGdual**  


### ${\color{violet} \text{Example}}$
```
%% load data.Z and data.Zphi and MCMC.initial_point
% The initial value was obtained as follows: 
    % set:  MCMC.initial_point = [auxinitR ;(data.Z-auxinitR.*data.Zphi)/2]
    %       where auxinitR = min(2,max(coeff(2)*(1:T)+coeff(1),0.5)) and coeff are the regression of vectors of data.Z on [data.Zphi (1:T)'.*data.Zphi]
    % run GibbsPGdual with MAP.method_augmentation = 'orthogonal' and chain_length = 2.5e7
    % save the last value of the Markov chain (output.lastsample)
    % set: MCMC.initial_point = output.lastsample;
load FranceDataSet1.mat 

MAP.lambdaR = 3.5 * std(data.Z);
MAP.lambdaO = 0.05;

param.T = length(data.Z);
param.displayglobal = 0; % could be '1' for the display of  graphical controls during the run
param.vecQ = [0.025 0.5 0.075];

MCMC.chain_length = 1.2e7;
MCMC.chain_burnin = ceil(0.5*MCMC.chain_length);
MCMC.gamma_init = 1e-7;
MCMC.target_ratio = 0.25;
MCMC.covariance = 'identity'; % or 'invert' or 'orthogonal'

outputPGdec = PGdec_github(data,MAP,MCMC,param);
```



# ${\color{blue} \text{Descrition of the data sets}}$

All the data sets contain: data.Z, data.Zphi and MCMC.initial_point

**FranceDataSet1.mat** the data from France, on $T=70$ days, available on Feb 7 2023.  
**FranceDataSet2.mat** the data from France, on $T=70$ days, available on May 1st 2022.  
**JapanDataSet1.mat** the data from Japan, on $T=70$ days, available on   
**IndiaDataSet1.mat** the data from India, on $T=70$ days, available on 


