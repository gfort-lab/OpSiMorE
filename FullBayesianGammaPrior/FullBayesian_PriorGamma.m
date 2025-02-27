function [output] = FullBayesian_PriorGamma(data,MCMC)
%% An algorithm to approximate the joint distribution
%% of (R_1,O_1, ..., R_T, O_T, lambda_R, lambda_O)
%%
%% when independent Gamma priors on (lambda_R, lambda_O)
%%
%%  June 2024 - developed by G. Fort
%%  Revised in January 2025 - by G. Fort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ the input variables
%% or
%% SET to the default value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the data Z and Phi 
Z = data.Z; % T x 1
Phi = data.Phi; % T x 1
T = size(Z,1);


% Read the parameters of the Gamma priors
alphaR = data.shapeR;
alphaO = data.shapeO;
betaR= data.inversescaleR;
betaO = data.inversescaleO;

% Length of the Markov chain
if isfield(MCMC,'chain_length')
    NbrMC  = MCMC.chain_length; 
else 
    NbrMC  = 1e7;
end

% Length of the burnin
if isfield(MCMC,'chain_burnin')
    burnin = MCMC.chain_burnin; 
else 
    burnin  = ceil(0.5*NbrMC);
end

% Initial point of the chain, for R
if isfield(MCMC,'initial_pointR')
    Rcurrent = MCMC.initial_pointR;    % T x 1
else 
    Rcurrent = ones(T,1);
end

% Initial point of the chain, for O
if isfield(MCMC,'initial_pointO')
    Ocurrent = MCMC.initial_pointO;    % T x 1
else 
    Ocurrent = zeros(T,1);
end


% Initial point of the chain, for LambdaR
if isfield(MCMC,'initial_pointLR')
    LambdaRcurrent = MCMC.initial_pointLR;    % 1 x 1
else 
    LambdaRcurrent = 3.5*std(Z);
end

% Initial point of the chain, for LambdaR
if isfield(MCMC,'initial_pointLO')
    LambdaOcurrent = MCMC.initial_pointLO;    % 1 x 1
else 
    LambdaOcurrent = 0.05;
end

% The Gamma parameters
if isfield(MCMC,'GammaO')
    GammaO = MCMC.GammaO;    % 1 x 1
else 
    GammaO = 1e3;
end
if isfield(MCMC,'GammaTildeR')
    GammaTildeR = MCMC.GammaTildeR;    % 1 x 1
else 
    GammaTildeR = 1e-12;
end

% Adaptation of GammaTildeR and GammaO
if isfield(MCMC,'adapt_frequency')
    frequency = MCMC.adapt_frequency;    % 1 x 1
else 
    frequency = 1e4; 
end
if isfield(MCMC,'target_ratioAR')
    ratioAR = MCMC.target_ratioAR;    % 1 x 1
else 
    ratioAR = 0.25; 
end


% Quantiles
if isfield(MCMC,'Qvec')
    vectQ = MCMC.Qvec; 
else 
    vectQ = [0.025 0.05 0.1 0.5 0.9 0.95 0.975];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE quantities
%% constant over iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Translated shape parameters
TshiftR = T+alphaR; % 1 x 1
TshiftO = T+alphaO; % 1 x 1

%% Define the shift 
shift = zeros(T,1);
if isfield(data,'Rinit')
     aux  = data.Rinit;
     shift(1) = aux(1)/4-aux(2)/2;
     shift(2) = aux(2)/4;
else  % as if aux(1)=aux(2)=1;
     shift(1) = -1/4;
     shift(2) = 1/4;    
end

%% Define the D matrix and it inverse
D = zeros(T,T);
v = [1/4 -2/4 1/4 zeros(1,T-3)];
vcirculant = toeplitz([v(1) fliplr(v(2:end))],v);
D(3:T,:) = vcirculant(1:T-2,:);
D(2,:) = [-1/2 1/4 zeros(1,T-2)];
D(1,:) = [1/4 zeros(1,T-1)];
clear v vcirculant

invD = inv(D);  % T x T
invDt = invD';  % T x T


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and memory allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the TildeR chain
TildeRcurrent = D*Rcurrent; % T x 1
% Initialize the current Poisson intensity
Intensitycurrent = Rcurrent.*Phi+Ocurrent; % T x 1
% Initialize the log-target density
StatRcurrent = sum(abs(TildeRcurrent+shift));   % 1 x 1
StatOcurrent = sum(abs(Ocurrent));  % 1 x 1
StatLRcurrent = (TshiftR-1)*log(LambdaRcurrent); % 1 x 1
StatLOcurrent = (TshiftO-1)*log(LambdaOcurrent); % 1 x 1
LogPicurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent; % 1 x 1
LogMarginalCurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-TshiftR*log(betaR+StatRcurrent)-TshiftO*log(betaO+StatOcurrent);
          

% Compute global and local acceptance rate
localARrateR = 0;
globalARrateR = 0;
localARrateO = 0;
globalARrateO = 0;


% Allocate memory for storage
StoreRchain = zeros(T,NbrMC+1);
StoreRchain(:,1) = Rcurrent;

StoreOchain = zeros(T,NbrMC+1);
StoreOchain(:,1) = Ocurrent;

StoreLRchain = zeros(1,NbrMC+1);
StoreLRchain(:,1) = LambdaRcurrent;

StoreLOchain = zeros(1,NbrMC+1);
StoreLOchain(:,1) = LambdaOcurrent;

StorelogPi = zeros(1,NbrMC+1);
StorelogPi(1) = LogPicurrent ; 

StorelogMarginal = zeros(1,NbrMC+1);
StorelogMarginal(1) = LogMarginalCurrent; 

% Sample i.i.d. gamma variables with parameters ((T+alpha),1).
RndGammaR = gamrnd((TshiftR)*ones(NbrMC,1), ones(NbrMC,1));
RndGammaO = gamrnd((TshiftO)*ones(NbrMC,1), ones(NbrMC,1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iterate a Metropolis-within-Gibbs algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=1:NbrMC

    %% Update the TildeR component

    % Propose a candidate
    GdtStep = TildeRcurrent+GammaTildeR*invDt*((Z./Intensitycurrent-1).*Phi);   % T x 1
    Drift = (max(abs(GdtStep+shift)-GammaTildeR*LambdaRcurrent,0).*sign(GdtStep+shift)-shift); % T x 1
    TildeRproposal = Drift+sqrt(2*GammaTildeR)*randn(T,1); % T x 1
    Rproposal = invD*TildeRproposal;    % T x 1
    Intensityproposal = Rproposal.*Phi+Ocurrent;    % T x 1
   
    % Accept or Reject the candidate
    test = sum((Rproposal>0).*((Intensityproposal>0)+(Z==0).*(Intensityproposal==0)));

    if test==T
        % compute the acceptance-rejection ratio
        GdtSteprev = TildeRproposal+GammaTildeR*invDt*((Z./Intensityproposal-1).*Phi);   % T x 1
        Driftrev = (max(abs(GdtSteprev+shift)-GammaTildeR*LambdaRcurrent,0).*sign(GdtSteprev+shift)-shift); % T x 1
        NumLogGauss = -sum((Driftrev-TildeRcurrent).^2)/(4*GammaTildeR); % 1 x 1
        DenomLogGauss = -sum((Drift-TildeRproposal).^2)/(4*GammaTildeR); % 1 x 1
        StatRproposal = sum(abs(TildeRproposal+shift)); % 1 x 1 
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*(betaR+StatRproposal)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent; % 1 x 1
        ARlogratio = LogPiproposal-LogPicurrent+NumLogGauss-DenomLogGauss; % 1 x 1
   
        % accept or reject
        if log(rand(1,1))<=ARlogratio
            % accept
            Intensitycurrent = Intensityproposal;   % T x 1
            TildeRcurrent = TildeRproposal; % T x 1
            Rcurrent = Rproposal;   % T x 1
            StatRcurrent = StatRproposal;
            LogPicurrent = LogPiproposal;   % 1 x 1
            localARrateR = localARrateR+1;    % 1 x 1
            globalARrateR = globalARrateR+1;  % 1 x 1
        end
     end
  

    %% Update the O component

    % Propose a candidate
    GdtStep = Ocurrent+GammaO*(Z./Intensitycurrent-1);   % T x 1
    Drift = max(abs(GdtStep)-GammaO*LambdaOcurrent,0).*sign(GdtStep);  % T x 1
    Oproposal = Drift+sqrt(2*GammaO)*randn(T,1); % T x 1
    Intensityproposal = Rcurrent.*Phi+Oproposal;    % T x 1
 
    % Accept or Reject the candidate
    test = sum((Intensityproposal>0)+((Z==0).*(Intensityproposal==0)));

    if test==T
        % compute the acceptance-rejection ratio
        GdtSteprev = Oproposal+GammaO*(Z./Intensityproposal-1);   % T x 1
        Driftrev = max(abs(GdtSteprev)-GammaO*LambdaOcurrent,0).*sign(GdtSteprev);  % T x 1
        NumLogGauss = -sum((Driftrev-Ocurrent).^2)/(4*GammaO); % 1 x 1
        DenomLogGauss = -sum((Drift-Oproposal).^2)/(4*GammaO); % 1 x 1
        StatOproposal = sum(abs(Oproposal));
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOproposal)+StatLRcurrent+StatLOcurrent; % 1 x 1
         
        ARlogratio = LogPiproposal-LogPicurrent+NumLogGauss-DenomLogGauss; % 1 x 1
 
        % accept or reject
        if log(rand(1,1))<=ARlogratio
            % accept
            Intensitycurrent = Intensityproposal;   % T x 1
            Ocurrent = Oproposal; % T x 1
            StatOcurrent = StatOproposal;   % T x 1
            LogPicurrent = LogPiproposal;   % 1 x 1
            localARrateO = localARrateO+1;    % 1 x 1
            globalARrateO = globalARrateO+1;  % 1 x 1
        end
    end


    %% Adapt GammaO and GammaTildeR

    if ((nn<=burnin) && (mod(nn,frequency) == 1 ) && (nn>=frequency))
            GammaTildeR =  GammaTildeR+(localARrateR/(frequency+1)-ratioAR)*GammaTildeR;
            localARrateR = 0;
            GammaO =  GammaO+(localARrateO/(frequency+1)-ratioAR)*GammaO;
            localARrateO = 0;
    end;

    %% Sample the LamdbaR and LambdaO chains

    LambdaRcurrent = RndGammaR(nn,1)/(betaR+StatRcurrent); % 1 x 1
    LambdaOcurrent = RndGammaO(nn,1)/(betaO+StatOcurrent); % 1 x 1
    StatLRcurrent = (TshiftR-1)*log(LambdaRcurrent); % 1 x 1
    StatLOcurrent = (TshiftO-1)*log(LambdaOcurrent); % 1 x 1
    auxLog =  Z'*log(Intensitycurrent)-sum(Intensitycurrent);   % 1 x 1
    LogPicurrent = auxLog-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent; % 1 x 1
    LogMarginalCurrent = auxLog-TshiftR*log(betaR+StatRcurrent)-TshiftO*log(betaO+StatOcurrent);
    

    %% Store the chains
    StoreRchain(:,nn+1) = Rcurrent; % T x 1
    StoreOchain(:,nn+1) = Ocurrent; % T X 1
    StoreLRchain(1,nn+1) = LambdaRcurrent; % 1 x 1
    StoreLOchain(1,nn+1) = LambdaOcurrent; % 1 x 1
    
   %% Store logPi
    StorelogPi(:,nn+1) = LogPicurrent; % 1 x 1
    StorelogMarginal(:,nn+1) = LogMarginalCurrent; % 1 x 1
end
clear RndGammaO RndGammaR


%%%%%%%%%%
%% Output 
%%%%%%%%%%
% The Gamma parameters 
output.GammaTildeR = GammaTildeR;
output.GammaO = GammaO;

% The chains LambdaR and LambdaO, with NO discarding
output.Lambdachain(1,:) = StoreLRchain(1,:); % 1 x NbrMC
output.Lambdachain(2,:) = StoreLOchain(1,:); % 1 x NbrMC

% The empirical mean of the R, O, LamdbaR, LambdaO chains (burnin phase discarded)
auxRchain = StoreRchain(:,burnin+1:NbrMC+1); % T x ...
auxOchain = StoreOchain(:,burnin+1:NbrMC+1); % T x ...
auxLRchain = StoreLRchain(1,burnin+1:NbrMC+1); % 1 x ...
auxLOchain = StoreLOchain(1,burnin+1:NbrMC+1); % 1 x ...
clear StoreLRchain StoreLOchain StoreRchain StoreOchain 

output.empirical_meanR = mean(auxRchain,2);
output.empirical_meanO = mean(auxOchain,2);
output.empirical_meanLR = mean(auxLRchain,2);
output.empirical_meanLO = mean(auxLOchain,2);

% The last value of the chains
output.lastR = auxRchain(:,end);  % T x 1
output.lastO = auxOchain(:,end); % T x 1
output.lastLR = auxLRchain(1,end); % 1 x 1
output.lastLO = auxLOchain(1,end); % 1 x 1

%  The Log distributions, with NO discarding
output.logPi = StorelogPi; % 1 x NbrMC+1
output.logMarginal = StorelogMarginal; %1 x NbrMC+1


% Some quantiles, for lambdaR and lambdaO chains (burnin phase, discarded)
QuantileLR = quantile(auxLRchain,vectQ');  % length(vectQ) x 1
QuantileLO = quantile(auxLOchain,vectQ');  % length(vectQ) x 1
output.quantilesLR = QuantileLR;
output.quantilesLO = QuantileLO;
clear auxLRchain auxLOchain 

% Some quantiles, for each component of the R and O chains (burnin phase, discarded)
MatrixQuantileR = zeros(length(vectQ),T);   
MatrixQuantileO = zeros(length(vectQ),T);
for tt=1:T
    MatrixQuantileR(:,tt) = quantile(auxRchain(tt,:),vectQ');  % length(vectQ) x 1
    MatrixQuantileO(:,tt) = quantile(auxOchain(tt,:),vectQ');  % length(vectQ) x 1
end
output.quantilesR = MatrixQuantileR;
output.quantilesO = MatrixQuantileO;
clear auxRchain auxOchain 

