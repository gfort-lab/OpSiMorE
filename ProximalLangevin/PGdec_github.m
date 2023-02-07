function [output] = PGdec_github(data,MAP,MCMC,param)
%% 
%% Developed by G. Fort, February 2023. 
%% 
%% PGdec algorithm
%%
%% 
%% INPUT
%% structure array data
%%
% data.Z:  T x 1 integers, the counts
% data.Zphi:  T x 1, the convolution for the Poisson intensity
%
%% structure array MAP
%%
% MAP.lambdaR: 1 x 1,  a positive regularization parameter for the R-block
% MAP.lambdaO: 1 x 1, a positive regularization parameter for the O-block
%
%% structure array MCMC
%%
% MCMC.chain_length: 1 x 1, the total length of the Markov chain (including
% burn-in)
%
% MCMC.chain_burnin: 1 x 1, the length of the burn-in
%
% MCMC.initial_point: (2T) x 1, the initial value of the chain  
%
% MCMC.target_ratio: 1x1, target value of the accept-reject ratio
%
% MCMC.gamma_init: 1x1, positive step size  
%
% MCMC.covariance: 1x1, either 'invert', or 'orthogonal' or 'identity'.
%%
%% structure array param 
% param.vectQ: a row vector, collects the order of the quantiles
%
% param.displayglobal:  1 x 1, a binary variable set to '1' when graphical
% controls are displayed.
%
% param.frequency: 1x1, frequency of the updates, during burn-in, for the step size gamma.
%%
%%%%%%%%%%
%% OUTPUT
%%%%%%%%%$
%%
%% structure output
%%
% output.empirical_mean: 2 x T. The empirical expectation of the bivariate chain , computed by discarding the samples of the burn-in phase.
%
% output.R_quantiles: q x T. The q quantiles for each of the T components of R, computed by discarding the samples of the burn-in phase.
% 
% output.O_quantiles: q x T. The q quantiles for each of the T components of O, computed by discarding the samples of the burn-in phase.
% 
% output.gamma: collects the successive values of the step size, adapted during the burn-in phase (and no more adapted, after the burn-in phase).
%
% output.logPi: collects the successive values of the log-density along iterations.
%
% output.lastsample: 2T x 1. the last sample of the Markov chain.



%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Read the inputs
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Read the inputs'); 
Zdata = data.Z;
Zphi = data.Zphi;
T = length(Zdata);
lambdaR = MAP.lambdaR;
lambdaO = MAP.lambdaO;
ratioAR = MCMC.target_ratio;
vectQ      = param.vectQ;


if isfield(MCMC,'chain_length')
    NbrMC  = MCMC.chain_length; 
else 
    NbrMC  = 1e7;
end

if isfield(MCMC,'chain_burnin')
    forget = MCMC.chain_burnin; 
else 
    forget  = ceil(0.5*NbrMC);
end

if isfield(MCMC,'initial_point')
     InitPointR = MCMC.initial_point(1:T);
     InitPointO = MCMC.initial_point(T+1:2*T);
else 
    InitPointR = ones(T,1);
    InitPointO = zeros(T,1);
end

if isfield(MCMC,'gamma_init')
    Gamma = MCMC.gamma_init; 
else 
    Gamma = 1e-7;
end
 

if isfield(param,'displayglobal')
    displayglobal = param.displayglobal;
else 
    displayglobal = 0;
end

if isfield(param,'frequency')
    frequency = frequency;
else 
    frequency = 1e4;
end


format long e

%%%%%%%%%%%%%%%
%% 1 - DECIDE 
%%%%%%%%%%%%%%%
if displayglobal == 1,
    % Initialize the counter of figures
    CntFig = 1;
end;
% How often controls on the behavior of the Markov Chain are displayed
displayMCfrequency = 1e5;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 - Definitions of variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Definition of variables '); 

% Define the (T-2) x T matrix  D2  
D2 = zeros(T-2,T); % (T-2) x T
v = [1/4 -2/4 1/4 zeros(1,T-3)];
vcirculant = toeplitz([v(1) fliplr(v(2:end))],v);
D2 = vcirculant(1:T-2,:);
D2 = D2*4/sqrt(6); % (T-2) x T

lambdaR = lambdaR*sqrt(6)/4;


% Define the three blocks from D2
D2block{1} = D2(1:3:T-2,:);   
D2block{2} = D2(2:3:T-2,:);   
D2block{3} = D2(3:3:T-2,:);  

% Define a projector operator
Projblock{1} = eye(T)-D2block{1}'*D2block{1};    % T x T
Projblock{2} = eye(T)-D2block{2}'*D2block{2};   % T x T
Projblock{3} = eye(T)-D2block{3}'*D2block{3};   % T x T


% Define the (2T-2) x (2T) matrix A
A = zeros(2*T-2,2*T);
A(1:T-2,1:T) = D2;  % (T-2) x T
A(T-1:2*T-2,T+1:2*T) = (lambdaO/lambdaR)*diag(Zphi);  % T x T

% Define the covariance matrices 
    % Define the barA matrix (2T) x (2T)
    barA = zeros(2*T,2*T);
    barA(3:2*T,:) = A;    % (2T-2) x 2T
    
    barA(1,1) =1;
    barA(2,1) = -2;
    barA(2,2) =1;
    barA(2,:) = barA(2,:)/sqrt(5);

switch MCMC.covariance
    case 'identity'
        Cov = eye(2*T);   % (2T) x (2T)
        SqrtCov = eye(2*T); % (2T) x (2T)
        InvCov = eye(2*T); % (2T) x (2T)

    case 'orthogonal'
        for ii=2:-1:1
            barA(ii,1:T) = (barA(ii,1:T)'-barA(ii+1:T,1:T)'*inv(barA(ii+1:T,1:T)*barA(ii+1:T,1:T)')*barA(ii+1:T,1:T)*barA(ii,1:T)')';
            barA(ii,1:T) = barA(ii,1:T)/norm(barA(ii,1:T));
        end;
        Cov = inv(barA)*(inv(barA))';  % (2T) x (2T)
        SqrtCov = sqrtm(Cov);  % (2T) x (2T)
        InvCov = barA'*barA; % (2T) x (2T)

    case 'invert'
        Cov = inv(barA)*(inv(barA))';  % (2T) x (2T)
        SqrtCov = sqrtm(Cov);  % (2T) x (2T)
        InvCov = barA'*barA; % (2T) x (2T)
end


% Where the data are positive
ZdataPos = Zdata>0; % T x 1




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3- Prepare the STORAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Allocate memory for storage'); 

% Store the Markov chain 
StoreMarkovChainR = zeros(T,1+NbrMC);  % T x (1+NbrMC)  
StoreMarkovChainO = zeros(T,1+NbrMC);  % T x (1+NbrMC) -- NOT normalized

% A binary vector, '1' when a proposed point is accepted
VectorAccept =  zeros(1,NbrMC);

% The successive values of Gamma during the burn-in phase
Gamma_store = zeros(1,1+NbrMC);

% The successive values of logpi
logpi_store = zeros(1,1+NbrMC);

% Sample all the Gaussian r.v. for the proposal step
GaussRnd = SqrtCov*randn(2*T,NbrMC); 

% Sample Uniform r.v. for the acceptance-rejection step
UnifRnd = rand(1,NbrMC);

% Sample a r.v. in (1,2,3) for selecting the block of D2
BlockRnd = ceil(3*rand(1,NbrMC)); 


%% Display the data
if displayglobal==1,
    %% Display the counts Z_t 
    figure(CntFig);
    clf;
    plot(1:T,Zdata,'ro-');
    title('Daily counts Zt');
    
    CntFig = CntFig+1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4- Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Initialization of the MCMC sampler');


% Initialization of the current value of the chain 
InitPointOn =  InitPointO./Zphi; % T x 1
currentMC = [InitPointR; InitPointOn];   % (2T) x 1  

% Store the initial value of the chain, in the original space
StoreMarkovChainR(:,1) = InitPointR;   % T x 1 
StoreMarkovChainO(:,1) = InitPointO;   % T x 1  -- NOT normalized

% The log-Like, the penalty term and the log-density at the current point
    % the penalty term 
NegPenaltyCurrent = -lambdaR*sum(abs(A*currentMC));   % 1 x 1
    % the log-likelihood term
auxintensity = Zphi.*(currentMC(1:T)+currentMC(T+1:2*T));   % T x 1    
LogLikeCurrent = Zdata(ZdataPos)'*log(auxintensity(ZdataPos))-sum(auxintensity); % 1 x 1
    % log pi
logpiCurrent = LogLikeCurrent+NegPenaltyCurrent;    % 1 x 1
logpi_store(1,1) = logpiCurrent;

% Store the initial learning rate
Gamma_store(1,1) = Gamma;   % 1 x 1

if displayglobal == 1
     fprintf('\n \t \t Initial value of the log target density term: %f,', logpiCurrent);  
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5- Run the Markov chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Loop of the MCMC sampler \n');

tic

for nn = 1:NbrMC
    % display or not some controls every 'displayMCfrequency' iterations
    display = displayglobal*mod(nn,displayMCfrequency);

    %% Select one of the block
    ell = BlockRnd(nn);
    D2local = D2block{ell};
    D2localt = D2local';
    Projlocal = Projblock{ell};


    %% propose a candidate for R and On 
    Grad(1:T,1) = -Zphi;
    aux = currentMC(T+1:2*T);
    Grad(ZdataPos,1) = Zdata(ZdataPos)./(currentMC(ZdataPos)+aux(ZdataPos))+Grad(ZdataPos,1);
    Grad(T+1:2*T,1) = Grad(1:T,1);  % T x 1
    
    GradStep = currentMC+Gamma*Grad;   % (2T) x 1

    auxprox = D2local*GradStep(1:T); % x x 1;
    Drift(1:T,1) = Projlocal*GradStep(1:T)+D2localt*(max(abs(auxprox)-Gamma*lambdaR,0).*sign(auxprox));
    Drift(T+1:2*T,1) = max(abs(GradStep(T+1:2*T))-Gamma*lambdaO*Zphi,0).*sign(GradStep(T+1:2*T));
    
    Proposal = Drift + sqrt(2*Gamma)*GaussRnd(:,nn);    % (2T) x 1
    

    %% Test if the proposal is in the domain D
    intensity = Zphi.*(Proposal(1:T)+Proposal(T+1:2*T));    % T X 1
    test1 = (intensity)>0;    % T x 1
    test2 = ((intensity)==0).*(Zdata==0); % T x 1
    auxsum = sum((Proposal(1:T)>=0).*(test1+test2));
    testsign = (auxsum == T); 

    % Display
    if display==1      
        %% Construction of one proposition, O-part
        figure(CntFig+1);
        clf;
        plot(1:T,InitPointOn,'c-','LineWidth',2);
        hold on;
        grid on;
        plot(1:T,currentMC(T+1:2*T),'r--','LineWidth',3);
        plot(1:T,GradStep(T+1:2*T),'gs','LineWidth',2);
        plot(1:T,Drift(T+1:2*T),'bo','LineWidth',2);
        plot(1:T,Proposal(T+1:2*T),'k--','LineWidth',2);
        title('Onorm : one iteration mechanism');
        legend('Init','Current','Gdt step','Gdt-Prox step','Proposal','location','best');        
        
        %% Construction of one proposition, R-part
        figure(CntFig+2);
        clf;
        plot(1:T,InitPointR(1:T),'c-','LineWidth',2);
        hold on;
        grid on;
        plot(1:T,currentMC(1:T),'r--','LineWidth',3);
        plot(1:T,GradStep(1:T),'gs','LineWidth',2);
        plot(1:T,Drift(1:T),'bo','LineWidth',2);
        plot(1:T,Proposal(1:T),'k--','LineWidth',2);
        title('R: one iteration mechanism');
        legend('Init','Current','Gdt step','Gdt-Prox step','Proposal','location','best');        
        
        
        %% Display examples of points, O-part
        figure(CntFig+3);
        plot(1:T,Zdata,'co-');
        hold on; grid on ; 
        plot(1:T,Zdata-(Zphi.*currentMC(T+1:2*T)),'k-');
        legend('Data','Data - Ocurr');
        title('Denoised data');
       
        %% Display examples of points, R-part
        figure(CntFig+4);
        plot(1:T,currentMC(1:T),'k-');
        hold on; grid on ; 
        title('R: many paths');
       
             
        %% Display the current intensity 
        figure(CntFig+5);
        plot(1:T,Zphi.*(currentMC(1:T)+currentMC(T+1:2*T)),'k-','Linewidth',2);
        grid on ; 
        hold on;
        plot(1:T,Zdata,'ro--','Linewidth',2);
        legend('Intensity','Data Z');
        title('Current intensity of the Poisson distribution');
        
        %% Print quantitative informations
        fprintf('\n At iteration %f',nn);
        fprintf('\n \t \t proposal in the set D: %f',testsign);
        fprintf('\n \t \t nbr constraints OK: %f', auxsum);
        fprintf('\n \t \t percent of accepted move: %1.4e', sum(VectorAccept(1,1:(nn-1)))/(nn-1)); 
       
    end;% of the display
      
      
    
    
    % if the proposed point has a positive probability
    if testsign == 1  
        %% compute the acceptance-rejection log-probability (log-alpha)
        % computation of the means of the Gaussian distribution
        Gradrev(1:T,1) = -Zphi;    % T x 1
        aux = Proposal(T+1:2*T);    % T x 1
        Gradrev(ZdataPos,1) = Zdata(ZdataPos)./(Proposal(ZdataPos)+aux(ZdataPos))+Gradrev(ZdataPos,1);    % T x 1
        Gradrev(T+1:2*T,1) = Gradrev(1:T,1);  % T x 1
      
        GradSteprev = Proposal+Gamma*Gradrev;   % (2T) x 1
        
        auxproxrev = D2local*GradSteprev(1:T); % x x 1;
        Driftrev(1:T,1) = Projlocal*GradSteprev(1:T)+D2localt*(max(abs(auxproxrev)-Gamma*lambdaR,0).*sign(auxproxrev));
        Driftrev(T+1:2*T,1) = max(abs(GradSteprev(T+1:2*T))-Gamma*lambdaO*Zphi,0).*sign(GradSteprev(T+1:2*T));
       
        % log-density, numerator
            % the penalty
        NegPenaltyProposal = -lambdaR*sum(abs(A*Proposal));   % 1 x 1
            % the log-likelihood term
        auxintensity = Zphi.*(Proposal(1:T)+Proposal(T+1:2*T));   % T x 1    
        LogLikeProposal = Zdata(ZdataPos)'*log(auxintensity(ZdataPos))-sum(auxintensity); % 1 x 1
            % logpi
        logpiProposal = NegPenaltyProposal+LogLikeProposal;
        
        % log-proposal, numerator
        QuadNum = (currentMC-Driftrev)'*InvCov*(currentMC-Driftrev);
        logGaussNum =  -QuadNum/(4*Gamma);
        
        % log-proposal, denominator
        QuadDenom = (Proposal - Drift)'*InvCov*(Proposal - Drift);
        logGaussDenom =  - QuadDenom/(4*Gamma);
     
        % log-ratio
        logalpha = min(0,logpiProposal-logpiCurrent +logGaussNum-logGaussDenom);
   
               
        %% test the value of the ratio
        if log(UnifRnd(1,nn))<=logalpha
            logpiCurrent = logpiProposal;
            NegPenaltyCurrent = NegPenaltyProposal;
            LogLikeCurrent = LogLikeProposal;
            currentMC = Proposal;   % (2T) x 1
            VectorAccept(1,nn) = 1;    
       end;
    end; % of testsign==1
    
    
    %% Adapt the step size Gamma during the burnin phase
    if ((mod(nn,frequency) ==1 )&& (nn>=frequency)&&(nn<forget)) 
        Gamma_aux =  Gamma+(sum(VectorAccept(1,nn-frequency:(nn-1)))/(frequency-1)-ratioAR)*Gamma;
        if display == 1
            fprintf('\n \t \t Gamma : current %1.4e \t  next %1.4e',Gamma, Gamma_aux);
            fprintf('\n \t \t Local acceptance rate is %1.4e',sum(VectorAccept(1,(nn-frequency):(nn)))/(frequency)); 
        end;
        Gamma = Gamma_aux;   
    end;
    Gamma_store(1,nn+1) = Gamma;    
        
  %%----------------------------------
  %% Store the chain and the energies
  %%-----------------------------------
   
    StoreMarkovChainR(:,nn+1) = currentMC(1:T);  % T x 1  -- in the original space
    StoreMarkovChainO(:,nn+1) = Zphi.*currentMC(T+1:2*T);  % T x 1 -- in the original space, not normalized
    logpi_store(nn+1) = logpiCurrent;   % 1 x 1
     
end;    % loop over the iterations of the Markov Chain
toc

clear GaussRnd UnifRnd BlockRnd

auxR = StoreMarkovChainR(:,forget+1:NbrMC);
auxO = StoreMarkovChainO(:,forget+1:NbrMC);
clear StoreMarkovChainO
clear StoreMarkovChainR

StoreMarkovChainR = auxR;
StoreMarkovChainO = auxO;
clear auxR auxO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of the quantiles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Computation of the quantiles');
MatrixQuantileR = zeros(length(vectQ),T);
MatrixQuantileO = zeros(length(vectQ),T);
for tt=1:T,
    MatrixQuantileR(:,tt) = quantile(StoreMarkovChainR(tt,:),vectQ');  % length(vectQ) x 1
    MatrixQuantileO(:,tt) = quantile(StoreMarkovChainO(tt,:),vectQ');  % length(vectQ) x 1
end;


%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE THE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Prepare the output');

output.empirical_mean = [mean(StoreMarkovChainR(:,:),2); mean(StoreMarkovChainO(:,:),2)];

output.R_quantiles = MatrixQuantileR;
output.O_quantiles = MatrixQuantileO;

output.gamma = Gamma_store(1:forget);

output.logPi = logpi_store;

output.lastsample = [StoreMarkovChainR(:,end); StoreMarkovChainO(:,end)];
    
    
  
  
    
   