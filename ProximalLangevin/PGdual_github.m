function [generic, additional] = PGdual_github(data,MAP,MCMC,param)
%% 
%% Developed by G. Fort, January 2023. 
%% 
%% PGdual algorithm
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
% MAP.method_augmentation: 1 x 1, either 'orthogonal', or 'invert'
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
%%%%%%%%%%
%%
%% structure generic
%%
% generic.empirical_mean: 2 x T. The empirical expectation of the bivariate chain , computed by discarding the samples of the burn-in phase.
%
% generic.R_quantiles: q x T. The q quantiles for each of the T components of R, computed by discarding the samples of the burn-in phase.
% 
% generic.O_quantiles: q x T. The q quantiles for each of the T components of O, computed by discarding the samples of the burn-in phase.
% 
% additional_gamma: collects the successive values of the step size, adapted during the burn-in phase (and no more adapted, after the burn-in phase).
%
% additional.logPi: collects the successive values of the log-density along iterations.



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

if isfield(MAP,'method_augmentation')
    method  = MAP.method_augmentation; 
else 
    method  = 'orthogonal';
end


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

% Define the barA matrix (2T) x (2T)
barA = zeros(2*T,2*T);
barA(3:T,1:T) = D2;    % (T-2) x T
barA(T+1:2*T,T+1:2*T) = (lambdaO/lambdaR)*diag(Zphi);  % T x T

barA(1,1) =1;
barA(2,1) = -2;
barA(2,2) =1;
barA(2,:) = barA(2,:)/sqrt(5);


switch method
    case 'orthogonal'
        for ii=2:-1:1
            barA(ii,1:T) = (barA(ii,1:T)'-barA(ii+1:T,1:T)'*inv(barA(ii+1:T,1:T)*barA(ii+1:T,1:T)')*barA(ii+1:T,1:T)*barA(ii,1:T)')';
        end;
    
end


% Compute the inverse of barA
invbarA = zeros(2*T,2*T);   % (2T) x (2T)
invbarA(1:T,1:T) = inv(barA(1:T,1:T));
invbarA(T+1:2*T,T+1:2*T) = (lambdaR/lambdaO)*diag(1./Zphi);  

% Compute the transpose of the inverse of barA
invbarAt = invbarA';    % (2T) x (2T) 

% Where the data are positive
ZdataPos = Zdata>0; % T x 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3- Prepare the STORAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Allocate memory for storage'); 

% Store the Markov chain 
StoreMarkovChainR = zeros(T,1+NbrMC);  % T x (1+NbrMC)  -- in the original space
StoreMarkovChainO = zeros(T,1+NbrMC);  % T x (1+NbrMC) -- in the original space, NOT normalized

% A binary vector, '1' when a proposed point is accepted
VectorAccept =  zeros(1,NbrMC);

% The successive values of Gamma during the burn-in phase
Gamma_store = zeros(1,1+NbrMC);

% The successive values of logpi
logpi_store = zeros(1,1+NbrMC);

% Sample all the Gaussian r.v. for the proposal step
GaussRnd = randn(2*T,NbrMC); 

% Sample Uniform r.v. for the acceptance-rejection step
UnifRnd = rand(1,NbrMC);


%% Display the data
if displayglobal==1,
    %% Display the counts Z_t 
    figure(CntFig);
    clf;
    plot(datetime('today')-T+1:datetime('today'),Zdata,'ro-');
    title('Daily counts Zt');
    
    CntFig = CntFig+1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4- Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Initialization of the MCMC sampler');

% Initialization of the Markov chain in the image space
InitPointOn = InitPointO./Zphi; % T x 1
AuxInitChainTilde = barA * [InitPointR; InitPointOn];   % (2T) x 1


% Initialization of the current value of the chain 
currentMCtilde = AuxInitChainTilde; % (2T) x 1 
currentMCorig = [InitPointR; InitPointOn];   % (2T) x 1  

% Store the initial value of the chain, in the original space
StoreMarkovChainR(:,1) = InitPointR;   % T x 1 
StoreMarkovChainO(:,1) = InitPointO;   % T x 1  -- NOT normalized


% The log-Like, the penalty term and the log-density at the current point
    % the penalty term 
NegPenaltyCurrent = -lambdaR*sum(abs(currentMCtilde(3:2*T)));   % 1 x 1
    % the log-likelihood term
auxintensity = Zphi.*(currentMCorig(1:T)+currentMCorig(T+1:2*T));   % T x 1    
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
    
    %% propose a candidate for Rtilde and Ontilde
    currentMCorig = invbarA*currentMCtilde;
    Grad(1:T,1) = Zdata./(currentMCorig(1:T)+currentMCorig(T+1:2*T))-Zphi;    % T x 1
    Grad(T+1:2*T,1) = Grad(1:T,1);  % T x 1
    Grad = invbarAt*Grad;   % (2T) x 1
    
    GradStep = currentMCtilde+Gamma*Grad;   % (2T) x 1
    Drift(1:2,1) = GradStep(1:2);
    Drift(3:2*T,1) = max(abs(GradStep(3:2*T))-Gamma*lambdaR,0).*sign(GradStep(3:2*T));
    
    Proposal = Drift + sqrt(2*Gamma)*GaussRnd(:,nn);    % (2T) x 1
    

    %% Test if the proposal is in the domain D_T
    ProposalOrig = invbarA*Proposal;    % (2T) x 1
    intensity = Zphi.*(ProposalOrig(1:T)+ProposalOrig(T+1:2*T));    % T X 1
    test1 = (intensity)>0;    % T x 1
    test2 = ((intensity)==0).*(Zdata==0); % T x 1
    auxsum = sum((ProposalOrig(1:T)>=0).*(test1+test2));
    testsign = (auxsum == T); 

    % Display
    if display==1      
        %% Construction of one proposition, O-part
        figure(CntFig+1);
        clf;
        plot(1:T,AuxInitChainTilde(T+1:2*T),'c-','LineWidth',2);
        hold on;
        grid on;
        plot(1:T,currentMCtilde(T+1:2*T),'r--','LineWidth',3);
        plot(1:T,GradStep(T+1:2*T),'gs','LineWidth',2);
        plot(1:T,Drift(T+1:2*T),'bo','LineWidth',2);
        plot(1:T,Proposal(T+1:2*T),'k--','LineWidth',2);
        title('Onorm Tilde: one iteration mechanism');
        legend('Init','Current','Gdt step','Gdt-Prox step','Proposal','location','best');        
        
        %% Construction of one proposition, R-part
        figure(CntFig+2);
        clf;
        plot(1:T,AuxInitChainTilde(1:T),'c-','LineWidth',2);
        hold on;
        grid on;
        plot(1:T,currentMCtilde(1:T),'r--','LineWidth',3);
        plot(1:T,GradStep(1:T),'gs','LineWidth',2);
        plot(1:T,Drift(1:T),'bo','LineWidth',2);
        plot(1:T,Proposal(1:T),'k--','LineWidth',2);
        title('R Tilde: one iteration mechanism');
        legend('Init','Current','Gdt step','Gdt-Prox step','Proposal','location','best');        
        
        
        %% Display examples of points, O-part
        figure(CntFig+3);
        plot(datetime('today')-T+1:datetime('today'),Zdata,'co-');
        hold on; grid on ; 
        plot(datetime('today')-T+1:datetime('today'),Zdata-(Zphi.*currentMCorig(T+1:2*T)),'k-');
        legend('Data','Data - Ocurr');
        title('Denoised data');
       
        %% Display examples of points, R-part
        figure(CntFig+4);
        plot(datetime('today')-T+1:datetime('today'),currentMCorig(1:T),'k-');
        hold on; grid on ; 
        title('R-orig: many paths');
       
             
        %% Display the current intensity 
        figure(CntFig+5);
        plot(datetime('today')-T+1:datetime('today'),Zphi.*(currentMCorig(1:T)+currentMCorig(T+1:2*T)),'k-','Linewidth',2);
        grid on ; 
        hold on;
        plot(datetime('today')-T+1:datetime('today'),Zdata,'ro--','Linewidth',2);
        legend('Intensity','Data Z');
        title('Current intensity of the Poisson distribution');
        
        %% Print quantitative informations
        fprintf('\n At iteration %f',nn);
        fprintf('\n \t \t proposal in the set DT: %f',testsign);
        fprintf('\n \t \t nbr constraints OK: %f', auxsum);
        fprintf('\n \t \t percent of accepted move: %1.4e', sum(VectorAccept(1,1:(nn-1)))/(nn-1)); 
       
    end;% of the display
      
      
    
    
    % if the proposed point has a positive probability
    if testsign == 1  
        %% compute the acceptance-rejection log-probability (log-alpha)
        % computation of the means of the Gaussian distribution
        Gradrev(1:T,1) = Zdata./(ProposalOrig(1:T)+ProposalOrig(T+1:2*T))-Zphi;    % T x 1
        Gradrev(T+1:2*T,1) = Gradrev(1:T,1);  % T x 1
        Gradrev = invbarAt*Gradrev;   % (2T) x 1
    
        GradSteprev = Proposal+Gamma*Gradrev;   % (2T) x 1
        Driftrev(1:2,1) = GradSteprev(1:2,1);
        Driftrev(3:2*T,1) = max(abs(GradSteprev(3:2*T,1))-Gamma*lambdaR,0).*sign(GradSteprev(3:2*T,1));
           
        % log-density, numerator
            % the penalty
        NegPenaltyProposal = -lambdaR*sum(abs(Proposal(3:2*T)));   % 1 x 1
            % the log-likelihood term
        auxintensity = Zphi.*(ProposalOrig(1:T)+ProposalOrig(T+1:2*T));   % T x 1    
        LogLikeProposal = Zdata(ZdataPos)'*log(auxintensity(ZdataPos))-sum(auxintensity); % 1 x 1
            % logpi
        logpiProposal = NegPenaltyProposal+LogLikeProposal;
        
        % log-proposal, numerator
        QuadNum = (currentMCtilde-Driftrev)'*(currentMCtilde-Driftrev);
        logGaussNum =  -QuadNum/(4*Gamma);
        
        % log-proposal, denominator
        QuadDenom = (Proposal - Drift)'*(Proposal - Drift);
        logGaussDenom =  - QuadDenom/(4*Gamma);
     
        % log-ratio
        logalpha = min(0,logpiProposal-logpiCurrent +logGaussNum-logGaussDenom);
   
               
        %% test the value of the ratio
        if log(UnifRnd(1,nn))<=logalpha
            logpiCurrent = logpiProposal;
            NegPenaltyCurrent = NegPenaltyProposal;
            LogLikeCurrent = LogLikeProposal;
            currentMCtilde = Proposal;   % (2T) x 1
            currentMCorig = ProposalOrig; % (2T) x 1
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
   
    StoreMarkovChainR(:,nn+1) = currentMCorig(1:T);  % T x 1  -- in the original space
    StoreMarkovChainO(:,nn+1) = Zphi.*currentMCorig(T+1:2*T);  % T x 1 -- in the original space, not normalized
    logpi_store(nn+1) = logpiCurrent;   % 1 x 1
     
end;    % loop over the iterations of the Markov Chain
toc

clear GaussRnd UnifRnd 

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

generic.empirical_mean = [mean(StoreMarkovChainR(:,:),2); mean(StoreMarkovChainO(:,:),2)];

generic.R_quantiles = MatrixQuantileR;
generic.O_quantiles = MatrixQuantileO;

additional.gamma = Gamma_store(1:forget);

additional.logPi = logpi_store;


    
    
 
    
  
  
    
   