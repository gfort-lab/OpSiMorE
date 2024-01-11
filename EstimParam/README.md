 Given
- T observations $Z_1, \cdots, Z_T$ taking values in $\mathbb{Z}_{\geq 0}$
- T mean values $\Phi_1, \cdots, \Phi_T$ taking values in $\mathbb{R}_{\geq 0}$
- Two initial values $R_{-1}, R_0$ for the reproduction number  taking values in $\mathbb{R}_{>0}$

two distributions are defined.

**Model "no mixture"** $\pi$ 

A target distribution $\pi(\cdot; \lambda_R,\lambda_O)$ on $(\mathbb{R}_{>0} \times \mathbb{R})^{T}$. The log-density $\ln \pi(\cdot; \lambda_R,\lambda_O)$ is given by

> $$ 
> {\small \begin{split}
> (R_1, O_1, \cdots, R_T, O_T) & \mapsto \sum_{t=1}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
> & - \frac{\lambda_R}{4} \sum_{t=1}^T | R_t - 2 R_{t-1} + R_{t-2} | + T \ln \lambda_R   \\
> & - \lambda_O  \sum_{t=1}^T |O_t| + T \ln \lambda_O
\end{split} }
> $$

up to an additive constant. The support of the distribution is the set $\mathcal{D}$ defined by

> $$
 {\small \text{for} \ t=1, \cdots,T, \qquad R_t > 0; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi(\cdot; \lambda_R,\lambda_O)$ depends on two positive parameters $\lambda_R$ and $\lambda_O$.

**Model "mixture"**  $\pi_m$

A target distribution $\pi_m(\cdot; \lambda_R,\lambda_{O,O},\lambda_{O,1},\omega)$ on $(\mathbb{R}_{>0} \times \mathbb{R} \times \\{0,1\\} )^{T}$. The log-density $\ln \pi_m$ is given by

>  $$ 
 {\small \begin{split} (R_1, O_1, B_1, \cdots, R_T, O_T, B_T) & \mapsto \sum_{t=1}^T \Bigl( Z_t \ln( R_t \Phi_t+O_t) - (R_t \Phi_t + O_t) \Bigr)  \\
& - \frac{\lambda_R}{4} \sum_{t=1}^T | R_t - 2 R_{t-1} + R_{t-2} | + T \ln \lambda_R   \\
& - \lambda_{O,0}  \sum_{t=1}^T (1-B_t) |O_t| + (T-\sum_{t=1}^T B_t) \ln \lambda_{O,0} \\
& - \lambda_{O,1}  \sum_{t=1}^T B_t |O_t| + (\sum_{t=1}^T B_t) \ln \lambda_{O,1} \\
& + (T-\sum_{t=1}^T B_t) \ln (1-\omega) +  (\sum_{t=1}^T B_t) \ln \omega 
\end{split}}
> $$ 

up to an additive constant. The support of the distribution is the set $\mathcal{D}_m$ defined by 

> $$
> {\small \text{for} \ t=1, \cdots,T, \qquad R_t > 0; \qquad \qquad B_t \in \\{0,1\\}; \qquad \qquad  R_t \Phi_t + O_t > 0  \quad \text{if} \quad Z_t >0, \quad \text{or} \quad R_t \Phi_t + O_t \geq 0  \quad \text{if} \quad Z_t  \geq 0.}
> $$

$\pi_m(\cdot; \lambda_R,\lambda_{O,O},\lambda_{O,1},\omega)$ depends on three positive parameters $\lambda_R$, $\lambda_{O,0}$ and $\lambda_{O,1}$; and a weight $\omega \in (0,1)$.


## ${\color{blue} \text{GibbsPGdual\\_nomixture}}$

This MATLAB code runs a Metropolis-within-Gibbs sampler with target distribution $\pi(\cdot; \lambda_R,\lambda_{O})$ and returns a Monte Carlo approximation for each expectation. The proposal mechanism depends on design parameters: they are adapted during the burnin phase in order to target a given mean acceptance ratio. 

 >$$ 
 {\small \begin{split}  I_R(\lambda_R,\lambda_O) & := \frac{1}{4} \int_{\mathcal{D}} \ \sum_{t=1}^T |r_t - 2 r_{t-1} + r_{t-2}| \ \   \pi(r_1,o_1, \cdots, r_T, o_T; \lambda_R,\lambda_O) \ \mathrm{d}r_1 \mathrm{d} o_1 \cdots \mathrm{d} r_T \mathrm{d} o_T \\   
 I_O(\lambda_R,\lambda_O) &:= \int_{\mathcal{D}} \  \sum_{t=1}^T |o_t| \ \  \mathrm{d} \pi(r_1,o_1, \cdots, r_T, o_T; \lambda_R,\lambda_O)  \  \mathrm{d}r_1 \mathrm{d} o_1 \cdots \mathrm{d} r_T \mathrm{d} o_T.
 \end{split}}
 >$$

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
-  _GammaO_: 1x1, initial value of the step size when proposing a candidate for the $O_t$ variables; default value 1e3
-  _GammaTildeR_ : 1x1, initial value of the step size when proposing a candidate for the second derivative of the $R_t$ variables; default value 1e-12
-  _target\_ratioAR_ : 1x1, targeted acceptance ratio during the adaptation phase; default value 0.25
- _Qvec_ : vector of order of quantiles; default value (0.025 0.05 0.5 0.95 0.975)

  
### ${\color{violet} \text{Output structures}}$
A structure _output_ with fields
- _StatR_ : the Monte Carlo approximation of $I_R(\lambda_R,\lambda_O)$
- _StatO_ : the Monte Carlo approximation of $I_O(\lambda_R,\lambda_O)$
- _GammaTildeR_ : step size when proposing a candidate for the second derivative of the R_t variables
- _GammaO_ : step size when proposing a candidae for the O_t variables
- _empirical_meanR_ : Tx1, a Monte Carlo approximation of the expectation of $(R_1, \cdots, R_T)$ under the distribution of $\pi(\cdot; \lambda_R,\lambda_0)$
- _empirical_meanO_ : Tx1, a Monte Carlo approximation of the expectation of $(O_1, \cdots, O_T)$ under the  distribution of $\pi(\cdot; \lambda_R,\lambda_0)$
- _quantilesR_ : length(Qvec) x T, the quantiles of $R_1, \cdots, R_T$ under the marginal distributions of $\pi(\cdot; \lambda_R,\lambda_0)$
- _quantilesO_ : length(Qvec) x T, the quantiles of $O_1, \cdots, O_T$ under the marginal distributions of $\pi(\cdot; \lambda_R,\lambda_0)$
- _lastsampleR_ : Tx1, the last MCMC sample R
- _lastsampleO_ : Tx1, the last MCMC sample O
- _LogPi_ : 1xchain\_length, the values of $\log \pi$ along the MCMC iterations 


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

and a Monte Carlo approximation of the second distribution 

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

MCMC.Qvec= [0.025 0.05 0.1 0.5 0.9 0.95 0.975];

MCMC.chain_length = 1e7;
MCMC.chain_burnin = ceil(0.5*MCMC.chain_length);

MCMC.GammaTildeR = 1e-12;
MCMC.GammaO = 1e3;
MCMC.target_ratioAR = 0.25;

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
title('distribution of  \lambda_0 under \pi^{(2)}')
```

