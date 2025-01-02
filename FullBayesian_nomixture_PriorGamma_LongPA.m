function [output] = FullBayesian_nomixture_PriorGamma_Long(data,MCMC)
%% An algorithm to approximate the joint distribution
%% of (R_1,O_1, ..., R_T, O_T, lambda_R, lambda_O)
%%
%% when there is independent Gamma priors on (lambda_R, lambda_O)
%% case : no mixture model
%%
%%  June 2024 - developed by G. Fort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the data Z and Phi 
Z = data.Z; % T x 1
Phi = data.Phi; % T x 1
T = size(Z,1);

% Read the parameters of the Gamma priors
alphaR = data.shapeR;
alphaO = data.shapeO;
betaR= data.inversescaleR;
betaO = data.inversescaleO;

TshiftR = T+alphaR;
TshiftO = T+alphaO;

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

% Quantiles
if isfield(MCMC,'Qvec')
    vectQ = MCMC.Qvec; 
else 
    vectQ = [0.025 0.05 0.1 0.5 0.9 0.95 0.975];
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
    LambdaRcurrent = MCMC.initial_pointLR;    % T x 1
else 
   LambdaRcurrent = 3.5*std(Z);
end

% Initial point of the chain, for LambdaR
if isfield(MCMC,'initial_pointLO')
    LambdaOcurrent = MCMC.initial_pointLO;    % T x 1
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
frequency = 1e4;
if isfield(MCMC,'target_ratioAR')
    ratioAR = MCMC.target_ratioAR;    % 1 x 1
else 
    ratioAR = 0.25; 
end


% Define the shift 
shift = zeros(T,1);
if isfield(data,'Rinit')
     aux  = data.Rinit;
     shift(1) = aux(1)/4-aux(2)/2;
     shift(2) = aux(2)/4;
else  % as if aux(1)=aux(2)=1;
     shift(1) = -1/4;
     shift(2) = 1/4;    
end


% Define the D matrix
D = zeros(T,T);
v = [1/4 -2/4 1/4 zeros(1,T-3)];
vcirculant = toeplitz([v(1) fliplr(v(2:end))],v);
D(3:T,:) = vcirculant(1:T-2,:);
D(2,:) = [-1/2 1/4 zeros(1,T-2)];
D(1,:) = [1/4 zeros(1,T-1)];
clear v vcirculant

invD = inv(D);  % T x T
invDt = invD';  % T x T



% Initialize the TildeR chain
TildeRcurrent = D*Rcurrent; % T x T
% Initialize the current Poisson intensity
Intensitycurrent = Rcurrent.*Phi+Ocurrent; % T x 1
% Initialize the log-target density
StatRcurrent = sum(abs(TildeRcurrent+shift));
StatOcurrent = sum(abs(Ocurrent));
StatLRcurrent = (TshiftR-1)*log(LambdaRcurrent);
StatLOcurrent = (TshiftO-1)*log(LambdaOcurrent);
LogPicurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent;
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
    TildeRproposal = Drift+sqrt(2*GammaTildeR)*randn(T,1);
    Rproposal = invD*TildeRproposal;    % T x 1
    Intensityproposal = Rproposal.*Phi+Ocurrent;    % T x 1
   
    % Accept or Reject the candidate
    test = sum((Rproposal>0).*((Intensityproposal>0)+(Z==0).*(Intensityproposal==0)));

    if test==T
        % compute the acceptance-rejection ratio
        GdtSteprev = TildeRproposal+GammaTildeR*invDt*((Z./Intensityproposal-1).*Phi);   % T x 1
        Driftrev = (max(abs(GdtSteprev+shift)-GammaTildeR*LambdaRcurrent,0).*sign(GdtSteprev+shift)-shift); % T x 1
        NumLogGauss = -sum((Driftrev-TildeRcurrent).^2)/(4*GammaTildeR);
        DenomLogGauss = -sum((Drift-TildeRproposal).^2)/(4*GammaTildeR);
        StatRproposal = sum(abs(TildeRproposal+shift));
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*(betaR+StatRproposal)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent;
        LogMarginalProposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-TshiftR*log(betaR+StatRproposal)-TshiftO*log(betaO+StatOcurrent);
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
    Oproposal = Drift+sqrt(2*GammaO)*randn(T,1);
    Intensityproposal = Rcurrent.*Phi+Oproposal;    % T x 1
 
    % Accept or Reject the candidate
    test = sum((Intensityproposal>0)+((Z==0).*(Intensityproposal==0)));

    if test==T
        % compute the acceptance-rejection ratio
        GdtSteprev = Oproposal+GammaO*(Z./Intensityproposal-1);   % T x 1
        Driftrev = max(abs(GdtSteprev)-GammaO*LambdaOcurrent,0).*sign(GdtSteprev);  % T x 1
        NumLogGauss = -sum((Driftrev-Ocurrent).^2)/(4*GammaO);
        DenomLogGauss = -sum((Drift-Oproposal).^2)/(4*GammaO);
        StatOproposal = sum(abs(Oproposal));
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOproposal)+StatLRcurrent+StatLOcurrent;
        ARlogratio = LogPiproposal-LogPicurrent+NumLogGauss-DenomLogGauss; % 1 x 1
 
        % accept or reject
        if log(rand(1,1))<=ARlogratio
            % accept
            Intensitycurrent = Intensityproposal;   % T x 1
            Ocurrent = Oproposal; % T x 1
            StatOcurrent = StatOproposal;   % T x 1
            LogPicurrent = LogPiproposal;   % 1 x 1
            LogMarginalCurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-TshiftR*log(betaR+StatRcurrent)-TshiftO*log(betaO+StatOcurrent);
            localARrateO = localARrateO+1;    % 1 x 1
            globalARrateO = globalARrateO+1;  % 1 x 1
        end
    end


    %% Adapt GammaX 
    if ((nn<=burnin) && (mod(nn,frequency) == 1 ) && (nn>=frequency))
        % Adapt GammaTildeR and GammaO
            GammaTildeR =  GammaTildeR+(localARrateR/(frequency+1)-ratioAR)*GammaTildeR;
            localARrateR = 0;
            GammaO =  GammaO+(localARrateO/(frequency+1)-ratioAR)*GammaO;
            localARrateO = 0;
    end;
        
    %% Sample the LamdbaR and LambdaO chains
    LambdaRcurrent = RndGammaR(nn,1)/(betaR+StatRcurrent);
    LambdaOcurrent = RndGammaO(nn,1)/(betaO+StatOcurrent);
    StatLRcurrent = (TshiftR-1)*log(LambdaRcurrent);
    StatLOcurrent = (TshiftO-1)*log(LambdaOcurrent);
    LogPicurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-LambdaRcurrent*(betaR+StatRcurrent)-LambdaOcurrent*(betaO+StatOcurrent)+StatLRcurrent+StatLOcurrent;

    

    %% store the chain
    StoreRchain(:,nn+1) = Rcurrent;
    StoreOchain(:,nn+1) = Ocurrent;
    StoreLRchain(1,nn+1) = LambdaRcurrent;
    StoreLOchain(1,nn+1) = LambdaOcurrent;
    
   %% store logPi
    StorelogPi(:,nn+1) = LogPicurrent;
    StorelogMarginal(:,nn+1) = LogMarginalCurrent;
end

OK = 1 ; 
save('/Users/pabry/Documents/MATLAB/OKtmp2.mat','OK')

%%%%%%%%%%
%% Output
%%%%%%%%%%
% The Gamma parameters 
output.GammaTildeR = GammaTildeR;
clear GammaTildeR
output.GammaO = GammaO;
clear GammaO 

% The empirical mean of the R, O, LamdbaR, LambdaO chains (burnin phase, discarded)
auxRchain = StoreRchain(:,burnin+1:NbrMC+1);
clear StoreRchain
auxOchain = StoreOchain(:,burnin+1:NbrMC+1);
clear StoreOchain

OK = 1 ; 
save('/Users/pabry/Documents/MATLAB/OKtmp3.mat','OK')

auxLRchain = StoreLRchain(1,burnin+1:NbrMC+1);
output.Lambdachain(1,:) = StoreLRchain(1,:); % auxLRchain ;
clear StoreLRchain 

auxLOchain = StoreLOchain(1,burnin+1:NbrMC+1);
output.Lambdachain(2,:) = StoreLOchain(1,:); % auxLOchain ;
clear  StoreLOchain

output.empirical_meanR = mean(auxRchain,2);
output.empirical_meanO = mean(auxOchain,2);
output.lastR = auxRchain(:,NbrMC-burnin);
output.lastO = auxOchain(:,NbrMC-burnin);
% Some quantiles, for lambdaR and lambdaO chains (burnin phase, discarded)
MatrixQuantileR = zeros(length(vectQ),T);   
MatrixQuantileO = zeros(length(vectQ),T);
for tt=1:T
    MatrixQuantileR(:,tt) = quantile(auxRchain(tt,:),vectQ');  % length(vectQ) x 1
    MatrixQuantileO(:,tt) = quantile(auxOchain(tt,:),vectQ');  % length(vectQ) x 1
end
output.quantilesR = MatrixQuantileR;
output.quantilesO = MatrixQuantileO;
clear MatrixQuantileR MatrixQuantileO
clear auxRchain auxOchain

OK = 1 ; 
save('/Users/pabry/Documents/MATLAB/OKtmp4.mat','OK')

output.empirical_meanLR = mean(auxLRchain,2);
output.empirical_meanLO = mean(auxLOchain,2);
output.lastLR = auxLRchain(NbrMC-burnin);
output.lastLO = auxLOchain(NbrMC-burnin);
% Some quantiles, for each component of the R and O chains (burnin phase, discarded)
QuantileLR = quantile(auxLRchain,vectQ');  % length(vectQ) x 1
QuantileLO = quantile(auxLOchain,vectQ');  % length(vectQ) x 1
output.quantilesLR = QuantileLR;
output.quantilesLO = QuantileLO;
clear MatrixQuantileO QuantileLR
clear auxLRchain auxLOchain

output.logPi = StorelogPi; % (burnin+1:NbrMC+1);
clear StorelogPi
output.logMarginal = StorelogMarginal(burnin+1:NbrMC+1);
clear StorelogMarginal

OK = 1 ; 
save('/Users/pabry/Documents/MATLAB/OKtmp5.mat','OK')


