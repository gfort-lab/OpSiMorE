 Given
- T observations $Z_1, \cdots, Z_T$ taking values in $\mathbb{Z}_{\geq 0}$
- T mean values $\Phi_1, \cdots, \Phi_T$ taking values in $\mathbb{R}_{\geq 0}$
- Two initial values $R_{-1}, R_0$ for the reproduction number  taking values in $\mathbb{R}_{>0}$

a distribution $\pi$ is defined on $(\mathbb{R}\_{>0} \times \mathbb{R})^{T} \times \mathbb{R}\_{>0} \times \mathbb{R}_{>0}$

The log-density $\log \pi$ is given by 
> $$ 
> {\small \begin{split}
> (R_1, O_1, \cdots, R_T, O_T, \lambda_R, \lambda_O) & \mapsto \sum_{t=1}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
> & - \frac{\lambda_R}{4} \sum_{t=1}^T | R_t - 2 R_{t-1} + R_{t-2} | + T \ln \lambda_R   \\
> & - \lambda_O  \sum_{t=1}^T |O_t| + T \ln \lambda_O \\
> & -\beta_R \lambda_R + (\alpha_R -1) \ln \lambda_R \\
> & - \beta_O \lambda_O + (\alpha_O -1) \ln \lambda_O
\end{split} } 
> $$

up to an additive constant on the support of the distribution, and $-\infty$ otherwise. The support of the distribution is the set $\mathcal{D} \times \mathbb{R}\_{>0} \times \mathbb{R}\_{>0}$ where $\mathcal{D}$ is defined by

> $$
> {\small \begin{split}
> \text{for} \ t=1, \cdots, T :  \qquad & R_t > 0;  \\
> & \text{and} \quad   R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.
> \end{split}}
> $$

$\pi$ depends on four positive parameters $\alpha_R, \beta_R, \alpha_O,\beta_O$ which are (hyper)parameters of the Gamma priors on $\lambda_R$ and $\lambda_O$ (see [ICASSP 2025](<https://hal.science/hal-04695138>)). 




## ${\color{blue} \text{FullBayesian\\_PriorGamma}}$

This MATLAB code runs a Gibbs sampler with target distribution $\pi$: more precisely, the variables $(R_1,O_1, R_2, O_2, \cdots, R_T, O_T)$ are sampled via a Metropolis-within-Gibbs, and the variables $(\lambda_R, \lambda_O)$ are sampled via two independent Gamma distributions with shape and rate parameters given by Eq(13) in [ICASSP 2025](<https://hal.science/hal-04695138>).

It returns   XXXX 

The proposal mechanism of the Metropolis-within-Gibbs step, depends on design parameters: they are adapted during the burnin phase in order to target a given mean acceptance ratio. 

### ${\color{violet} \text{Input structures}}$
A structure _data_ with fields
- _Z_ : Tx1, the sequence $Z_1, \cdots, Z_T$
- _Phi_ : Tx1, the sequence $\Phi_1, \cdots, \Phi_T$
- _Rinit_ : 2x1, the initial values $R_{-1}$ and $R_0$
- _shapeR_ : 1x1, the value of $\alpha_R$, the shape parameter of the Gamma prior on $\lambda_R$
- _shapeO_ : 1x1, the value of $\alpha_O$, the shape parameter of the Gamma prior on $\lambda_O$
- _inversescaleR_ : 1x1, the value of $\beta_R$, the inverse scale parameter of the Gamma prior on $\lambda_R$
- _inversescaleO_ : 1x1, the value of $\beta_O$, the inverse scale parameter of the Gamma prior on $\lambda_O$

A structure _MCMC_ with fields
- _chain\_length_ : 1x1, number of MCMC iterations; default value 1e7
-  _chain\_burnin_ : 1x1, length of the burnin period; default value 0.5*_chain\_length_
-  _initial\_pointR_ : Tx1, initial value $(R_1, \cdots, R_T)$ of the R chain 
-  _initial\_pointO_ : Tx1, initial value $(O_1, \cdots, O_T)$ of the O chain
-  _initial\_pointLR_ : 1x1, initial value of the $\lambda_R$ chain; default value is $3.5  \mathrm{std}(Z)$
- _initial\_pointLO_ : 1x1, initial value of the $\lambda_O$ chain; default value is $0.05$
-  _GammaO_: 1x1, initial value of the step size when proposing a candidate for the $O_t$ variables; default value 1e3
-  _GammaTildeR_ : 1x1, initial value of the step size when proposing a candidate for the second derivative of the $R_t$ variables; default value 1e-12
-  _adapt_frequency_ : 1x1, how frequent the adaptation mechanism of the parameters $\gamma_{\tilde R}$ and $\gamma_O$ is; default value is 1e4 iterations
-  _target\_ratioAR_ : 1x1, targeted acceptance ratio during the adaptation mechanism; default value 0.25
- _Qvec_ : vector of order of quantiles; may be empty by setting _MCMC.Qvec = []_. Default value (0.025 0.05 0.1 0.5 0.9 0.95 0.975)

  
### ${\color{violet} \text{Output structures}}$
A structure _output_ with fields
- _GammaTildeR_ : 1x1, step size when proposing a candidate for the second derivative of the R_t variables
- _GammaO_ : 1x1, step size when proposing a candidate for the O_t variables
- _empirical_meanR_ : Tx1, a Monte Carlo approximation of the expectation of $(R_1, \cdots, R_T)$ under the distribution $\pi$ -- computed after discarding burn-in samples.
- _empirical_meanO_ : Tx1, a Monte Carlo approximation of the expectation of $(O_1, \cdots, O_T)$ under the  distribution $\pi$ -- computed after discarding burn-in samples.
- _empirical_meanLR_ : 1x1, a Monte Carlo approximation of the expectation of $\lambda_R$ under the distribution $\pi$ -- computed after discarding burn-in samples.
- _empirical_meanLO_ : 1x1, a Monte Carlo approximation of the expectation of $\lambda_O$ under the  distribution $\pi$ -- computed after discarding burn-in samples.
- _lastR_ : Tx1, the last MCMC sample $(R_1, \cdots, R_T)$
- _lastO_ : Tx1, the last MCMC sample $(O_1, \cdots, O_T)$
- _lastLR_ : 1x1, the last MCMC sample $\lambda_R$
- _lastLO_ : 1x1, the last MCMC sample $\lambda_O$
- _logPi_ : 1xchain\_length, the values of $\log \pi$ along the MCMC iterations
- _logMarginal_ : 1xchain\_length, the values of $\log \tilde \pi$ along the MCMC iterations
- 
and, if _MCMC.Qvec_ is not empty,
- _quantilesR_ : length(Qvec) x T, the quantiles of $R_1, \cdots, R_T$ under the marginal distributions of $\pi$
- _quantilesO_ : length(Qvec) x T, the quantiles of $O_1, \cdots, O_T$ under the marginal distributions of $\pi$
- _quantilesLR_ : length(Qvec) x 1, the quantiles of $\lambda_R$ under the marginal distributions of $\pi$
- _quantilesLO_ : length(Qvec) x 1, the quantiles of $\lambda_O$ under the marginal distributions of $\pi$



### ${\color{violet} \text{Example}}$
(see [camsap23 paper](https://hal.science/hal-04174245v2)) for details on data.Z, data.Phi, data.Rinit

```
%% load data.Z, data.Phi, data.Rinit and MCMC.initial_pointR, MCMC.initial_pointO
% data.Z is part of a time series downloaded from JHU repository
% data.Phi is built from this time series 
% data.Rinit is built from this time series 
% The initial value MCMC.initial_pointR of the vector R was obtained from: 
    % a code by [B. Pascal](https://bpascal-fr.github.io/), which computes the MAP of \pi  given a set of values for (\lambda_R, \lambda_O)
    % Here, $\lambda_R$ and $\lambda_O$ are fixed to 3.5 std(data.Z) and 0.05 respectively.
% The initial valueMCMC.initial_pointO  of O_t is chosen as a linear convex combination of Z_t and R_t Phi_t.

load FranceDataSet1.mat

data.LambdaR = 3.5*std(data.Z);
data.LambdaO = 0.05;

MCMC.Qvec= [0.025 0.05 0.1 0.5 0.9 0.95 0.975];

MCMC.chain_length = 1e7;
MCMC.chain_burnin = ceil(0.5*MCMC.chain_length);

MCMC.GammaTildeR = 1e-12;
MCMC.GammaO = 1e3;
MCMC.target_ratioAR = 0.25;
MCMC.adapt_frequency = 1e4;

[output] = GibbsPGdual_nomixture(data,MCMC);

% Estimates of the reproduction number:
% - point estimate via the expectation
% - credibility intervals at level 80%, 90% and 95%
figure(1)
clf 
T = size(data.Z,1);
plot(output.empirical_meanR,'k','LineWidth',2);
hold on
grid on
Qlower = output.quantilesR(1,:);
Qupper = output.quantilesR(7,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.2);
Qlower = output.quantilesR(2,:);
Qupper = output.quantilesR(6,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.4);
Qlower = output.quantilesR(3,:);
Qupper = output.quantilesR(5,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.6);
plot(-1,data.Rinit(1),'kd','MarkerSize',4,'MarkerFaceColor','k');
plot(0,data.Rinit(2),'kd','MarkerSize',4,'MarkerFaceColor','k');
yyaxis right
plot(1:T,data.Z,'-o','Color',[1 0 0 0.2],'MarkerSize',2);
caption = sprintf('France');
title(caption,'FontSize',10);
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';

% Estimates of the denoised data Z_t-O_t, from
% - a point estimate of O, via the expectation
% - credibility intervals at level 80%, 90% and 95%
figure(2)
clf 
plot(1:T,data.Z-output.empirical_meanO,'r','LineWidth',2);
hold on
grid on
plot(1:T,data.Z,'k--o','MarkerSize',2);
Qlower = data.Z'-output.quantilesO(1,:);
Qupper = data.Z'-output.quantilesO(7,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'m','EdgeColor','m','LineStyle','--','FaceAlpha',0.2);
Qlower = data.Z'-output.quantilesO(2,:);
Qupper = data.Z'-output.quantilesO(6,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'r','EdgeColor','r','LineStyle','--','FaceAlpha',0.4);
Qlower = data.Z'-output.quantilesO(3,:);
Qupper = data.Z'-output.quantilesO(5,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'r','EdgeColor','r','LineStyle','--','FaceAlpha',0.6);
caption = sprintf('France');
title(caption,'FontSize',10);
```  


## ${\color{blue} \text{FullBayesian\\_nomixture}}$

This MATLAB code runs a Metropolis-within-Gibbs sampler with target distribution proportional to 

> $$
> {\small (R_1,O_1, \cdots, R_T,O_T, \lambda_R, \lambda_0) \mapsto \pi(R_1,O_1, \cdots, R_T,O_T; \lambda_R, \lambda_O) \qquad \text{on} \quad \mathcal{D} \times (0,\infty) \times (0, \infty).}
> $$

 It returns a Monte Carlo approximation of quantiles and expectation of the distributions 

> $$
 {\small \begin{align} \pi^{(1)}: \quad (R_1,O_1, \cdots, R_T,O_T) & \mapsto \int_0^{\infty} \int_0^\infty   \  \pi(R_1,O_1, \cdots, R_T,O_T; \lambda_R, \lambda_O) \ \ \mathrm{d} \lambda_R \ \mathrm{d} \lambda_O \qquad \qquad \text{on $\mathcal{D}$}   \\
\pi^{(2)}: \quad  (\lambda_R, \lambda_O) & \mapsto \int_{\mathcal{D}} \ \pi(R_1,O_1, \cdots, R_T,O_T; \lambda_R, \lambda_O) \ \ \mathrm{d}r_1 \mathrm{d} o_1 \cdots \mathrm{d} r_T \mathrm{d} o_T \qquad \qquad  \text{on}  \quad (0, \infty) \times (0,\infty);   \end{align}}
 > $$ 

and a Monte Carlo approximation of the distribution $\pi^{(2)}$.

### ${\color{violet} \text{Input structures}}$
A structure _data_ with fields
- _Z_ : Tx1, the sequence $Z_1, \cdots, Z_T$
- _Phi_ : Tx1, the sequence $\Phi_1, \cdots, \Phi_T$
- Rinit : 2x1, the initial values $R_{-1}$ and $R_0$

A structure _MCMC_ with the same fields as in **GibbsPGdual\_nomixture**, and in addition
- _initial_pointLR_ : 1x1, initial value of the $\lambda_R$ chain
- _initial_pointLO_ : 1x1, initial value of the $\lambda_O$ chain

### ${\color{violet} \text{Output structures}}$
A structure _output_ with the same fields _GammaTildeR_ and _GammaO_ as in **GibbsPGdual\_nomixture** and
- _empirical_meanR_ : Tx1, a Monte Carlo approximation of the expectation of $(R_1, \cdots, R_T)$ under the  distribution $\pi^{(1)}$ 
- _empirical_meanO_ : Tx1, a Monte Carlo approximation of the expectation of $(O_1, \cdots, O_T)$ under the  distribution $\pi^{(1)}$
- _empirical_meanLR_ : 1x1, a Monte Carlo approximation of the expectation of $\lambda_R$ under the  distribution $\pi^{(2)}$ 
- _empirical_meanLO_ : 1x1, a Monte Carlo approximation of the expectation of $\lambda_O$ under the  distribution $\pi^{(2)}$
- _quantilesR_ : length(Qvec) x T, the quantiles of $R_1, \cdots, R_T$ under the marginal distributions of $\pi^{(1)}$
- _quantilesO_ : length(Qvec) x T, the quantiles of $O_1, \cdots, O_T$ under the marginal distributions of $\pi^{(2)}$
- _Lambdachain_ : 2xL, the bivariate Markov chain approximating $\pi^{(2)}$; the samples from the burnin period are discarded so that L = MCMC.chain\_length-MCMC.chain\_burnin.
- _lastR_ : T x 1, the last sample of the R-chain
- _lastO_ : T x 1, the last sample of the O-chain
- _lastLR_ : 1 x 1, the last sample of the $\lambda_R$ chain
- _lastLO_ : 1 x 1, the last sample of the $\lambda_O$ chain

### ${\color{violet} \text{Example}}$
(see [camsap23 paper](https://hal.science/hal-04174245v2)) for details on data.Z, data.Phi, data.Rinit
```
%% load data.Z, data.Phi, data.Rinit and MCMC.initial_pointR, MCMC.initial_pointO
% data.Z is part of a time series downloaded from JHU repository
% data.Phi is built from this time series 
% data.Rinit is built from this time series 
% The initial value MCMC.initial_pointR of the vector R was obtained from: 
    % a code by [B. Pascal](https://bpascal-fr.github.io/), which computes the MAP of \pi  given a set of values for (\lambda_R, \lambda_O)
    % Here, $\lambda_R$ and $\lambda_O$ are fixed to 3.5 std(data.Z) and 0.05 respectively.
% The initial value MCMC.initial_pointO  of O_t is chosen as a linear convex combination of Z_t and R_t Phi_t.

load FranceDataSet1.mat

MCMC.Qvec= [0.025 0.05 0.1 0.5 0.9 0.95 0.975];

MCMC.chain_length = 1e7;
MCMC.chain_burnin = ceil(0.5*MCMC.chain_length);

MCMC.GammaTildeR = 1e-12;
MCMC.GammaO = 1e3;
MCMC.target_ratioAR = 0.25;
MCMC.adapt_frequency = 1e4;

MCMC.initial_pointLR = 3.5*std(data.Z);
MCMC.initial_pointLO = 0.05; 

[outputFB] = FullBayesian_nomixture(data,MCMC);

% Estimates of the reproduction number:
% - point estimate via the expectation
% - credibility intervals at level 80%, 90% and 95%
figure(1)
clf 
T = size(data.Z,1);
plot(outputFB.empirical_meanR,'k','LineWidth',2);
hold on
grid on
Qlower = outputFB.quantilesR(1,:);
Qupper = outputFB.quantilesR(7,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.2);
Qlower = outputFB.quantilesR(2,:);
Qupper = outputFB.quantilesR(6,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.4);
Qlower = outputFB.quantilesR(3,:);
Qupper = outputFB.quantilesR(5,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'b','EdgeColor','b','LineStyle','--','FaceAlpha',0.6);
plot(-1,data.Rinit(1),'kd','MarkerSize',4,'MarkerFaceColor','k');
plot(0,data.Rinit(2),'kd','MarkerSize',4,'MarkerFaceColor','k');
yyaxis right
plot(1:T,data.Z,'-o','Color',[1 0 0 0.2],'MarkerSize',2);
caption = sprintf('France');
title(caption,'FontSize',10);
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';

% Estimates of the denoised data Z_t-O_t, from
% - a point estimate of O, via the expectation
% - credibility intervals at level 80%, 90% and 95%
figure(2)
clf 
plot(1:T,data.Z-outputFB.empirical_meanO,'r','LineWidth',2);
hold on
grid on
plot(1:T,data.Z,'k--o','MarkerSize',2);
Qlower = data.Z'-outputFB.quantilesO(1,:);
Qupper = data.Z'-outputFB.quantilesO(7,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'m','EdgeColor','m','LineStyle','--','FaceAlpha',0.2);
Qlower = data.Z'-outputFB.quantilesO(2,:);
Qupper = data.Z'-outputFB.quantilesO(6,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'r','EdgeColor','r','LineStyle','--','FaceAlpha',0.4);
Qlower = data.Z'-outputFB.quantilesO(3,:);
Qupper = data.Z'-outputFB.quantilesO(5,:);
fill([1:T fliplr(1:T)],[Qlower fliplr(Qupper)],'r','EdgeColor','r','LineStyle','--','FaceAlpha',0.6);
caption = sprintf('France');
title(caption,'FontSize',10);


% Histogram of the LambdaR, LambdaO
figure(3)
clf
subplot(2,1,1);
histogram(outputFB.Lambdachain(1,:),'Normalization','pdf');
title('distribution of \lambda_R under \pi^{(2)}')
subplot(2,1,2);
histogram(outputFB.Lambdachain(2,:),'Normalization','pdf');
title('distribution of  \lambda_O under \pi^{(2)}')
```


## ${\color{blue} \text{SAEM\\_nomixture}}$
This MATLAB code runs a _Stochastic Approximation Expectation Maximization_ (SAEM) algorithm (see Section 3 in [camsap23 paper](https://hal.science/hal-04174245v2), and references therein) in order to solve the optimization problem

>$$
{\small  \mathrm{argmax}  \int_{\mathcal{D}} \ \pi(R_1,O_1, \cdots, R_T,O_T; \lambda_R, \lambda_O) \ \ \mathrm{d}r_1 \mathrm{d} o_1 \cdots \mathrm{d} r_T \mathrm{d} o_T  \qquad (\lambda_R, \lambda_O) \in (0,\infty) \times (0, \infty).}
>$$

SAEM is an iterative algorithm: **SAEM\_nomixture** returns a sequence $\\{ (\lambda_R^k, \lambda_O^k), k \geq 0 \\}$ of parameters converging to a solution.


### ${\color{violet} \text{Input structures}}$

The structure _data_ with fields
- _Z_: Tx1, the sequence $Z_1, \cdots, Z_T$
- _Phi_ : Tx1, the sequence $\Phi_1, \cdots, \Phi_T$
- _Rinit_ : 2x1, the initial values $R_{-1}$ and $R_0$

The _SAEM_structure with fields
- _NbrIter_: 1x1, number of iterations of SAEM; default value is 5e5
- _LambdaRinit_: 1x1, initial value of the $\lambda_R$ sequence; default value is 3.5*std(data.Z)
- _LambdaOinit_: 1x1, initial value of the $\lambda_O$ sequence; default value is 0.05
- _pas\_vect\_R_: 1xSAEM.NbrIter, sequence of SAEM learning rate; default value is 0.05*((ones(1,10) 2*ones(1,100) 4./sqrt(100:100+(SAEM.NbrIter-110)));
- _pas\_vect\_O_: 1xSAEM.NbrIter, sequence of SAEM learning rate; default value is 0.5*(0.1*ones(1,10) 0.05*ones(1,100) 0.1./sqrt(100:100+(SAEM.NbrIter-110)))
- _controldisplay_: 1x1, a binary variable set to '1' for the display of the SAEM sequence every 5e2 iterations, and '0' otherwise.

The _MCMC_structure with fields
- _chain_length_: 1xSAEM.NbrIter, the length of the MCMC chain at each iteration of SAEM; default value is (1e7 5e6  3e3*ones(1,SAEM.NbrIter-2))
- _chain_burnin_: 1xSAEM.NbrIter, the length of the burnin phase when running the MCMC chain  at each iteration of SAEM; default value is 0.5*_MCMC.chain_length_
- _initial\_pointR_: Tx1, the initial value of the R-chain for the first iteration of SAEM
- _initial\_pointO_: Tx1, the initial value of the O-chain for the first iteration of SAEM
- _GammaTildeR_: 1x1, the initial value of the $\gamma_{\tilde R}$ parameter when running the MCMC chain for the first iteration of SAEM; default value is 1e-12
- _GammaO_: 1x1, the initial value of the $\gamma_O$ parameter when running the MCMC chain for the first iteration of SAEM; default value is 1e3
- _adapt\_frequency_:  1xSAEM.NbrIter, at eah iteration of SAEM, how often is the adaptation mechanism during the burnin phase; the devault value is max(MCMC.chain_burnin/500,500)
- _target\_ratioAR_: 1x1, the targeted mean acceptance ratio when adapting the parameters $\gamma_{\tilde R}$ and $\gamma_O$.
  


### ${\color{violet} \text{Output structures}}$
A structure _output_ with fields
- _LambdaRpath_: 1xSAEM.NbrIter, the $\lambda_R$ sequence produced by SAEM
- _LambdaOpath_: 1xSAEM.NbrIter, the $\lambda_O$ sequence produced by SAEM
- _GammaTildeRpath_: 1xSAEM.NbrIter, the sequence of  parameters $\gamma_{\tilde R}$ used in the MCMC chain at each iteration of SAEM
- _GammaOpath_: 1xSAEM.NbrIter, the sequence of parameters $\gamma_O$ used in the MCMC chain at each iteration of SAEM


### ${\color{violet} \text{Example}}$
(see [camsap23 paper](https://hal.science/hal-04174245v2)) for details on data.Z, data.Phi, data.Rinit
```
%% load data.Z, data.Phi, data.Rinit and MCMC.initial_pointR, MCMC.initial_pointO
% data.Z is part of a time series downloaded from JHU repository
% data.Phi is built from this time series 
% data.Rinit is built from this time series 
% The initial value MCMC.initial_pointR of the vector R was obtained from: 
    % a code by [B. Pascal](https://bpascal-fr.github.io/), which computes the MAP of \pi  given a set of values for (\lambda_R, \lambda_O)
    % Here, $\lambda_R$ and $\lambda_O$ are fixed to 3.5 std(data.Z) and 0.05 respectively.
% The initial value MCMC.initial_pointO  of O_t is chosen as a linear convex combination of Z_t and R_t Phi_t.

load FranceDataSet1.mat


SAEM.LambdaRinit = 3.5*std(data.Z); 
SAEM.LambdaOinit = 0.05;

SAEM.NbrIter = 5e5;

SAEM.pas_vect_R = 0.05*[ones(1,10) 2*ones(1,100) 4./sqrt(100:100+(SAEM.NbrIter-110))];
SAEM.pas_vect_O = 0.5*[0.1*ones(1,10) 0.05*ones(1,100) 0.1./sqrt(100:100+(SAEM.NbrIter-110))];

SAEM.controldisplay = 1;

MCMC.chain_length = [1e7 5e6  3e3*ones(1,SAEM.NbrIter-2)];
MCMC.chain_burnin =ceil(0.5*MCMC.chain_length);

MCMC.GammaTildeR = 1e-12;
MCMC.GammaO = 1e3; 

MCMC.target_ratioAR = 0.25;    
MCMC.adapt_frequency = max(MCMC.chain_burnin/500,500); 

outputSAEM = SAEM_nomixture(data,SAEM,MCMC);

%% Display the LambdaR sequence produced by SAEM
figure(1)
clf
plot(2:SAEM.NbrIter+1,outputSAEM.LambdaRpath(1,2:SAEM.NbrIter+1),);
grid on
caption = sprintf('SAEM sequence: \lambda_R');
title(caption,'FontSize',10);

%% Display the LambdaO sequence produced by SAEM
figure(2)
clf
plot(2:SAEM.NbrIter+1,outputSAEM.LambdaOpath(1,2:SAEM.NbrIter+1),);
grid on
caption = sprintf('SAEM sequence: \lambda_O');
title(caption,'FontSize',10);
```
