function [output] = FullBayesian_nomixture(data,MCMC)
%% An algorithm to approximate the joint distribution
%% of (R_1,O_1, ..., R_T, O_T, lamnda_R, lambda_O)
%%
%% proportional to pi
%% case : no mixture model
%%
%%  January 2024 - developed by G. Fort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the data Z and Phi 
Z = data.Z; % T x 1
Phi = data.Phi; % T x 1
T = size(Z,1);


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
    GammaO = 1e-7;
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
StatLRcurrent = T*log(LambdaRcurrent);
StatLOcurrent = T*log(LambdaOcurrent);
LogPicurrent = Z'*log(Intensitycurrent)-sum(Intensitycurrent)-LambdaRcurrent*StatRcurrent-LambdaOcurrent*StatOcurrent+StatLRcurrent+StatLOcurrent;


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

% Sample i.i.d. gamma variables with parameters ((T+1),1).
RndGamma = gamrnd((T+1)*ones(NbrMC,2), ones(NbrMC,2));

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
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*StatRproposal-LambdaOcurrent*StatOcurrent+StatLRcurrent+StatLOcurrent;
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
        LogPiproposal = Z'*log(Intensityproposal)-sum(Intensityproposal)-LambdaRcurrent*StatRcurrent-LambdaOcurrent*StatOproposal+StatLRcurrent+StatLOcurrent;
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


    %% Adapt GammaX 
    if ((nn<=burnin) && (mod(nn,frequency) == 1 ) && (nn>=frequency))
        % Adapt GammaTildeR and GammaO
            GammaTildeR =  GammaTildeR+(localARrateR/(frequency+1)-ratioAR)*GammaTildeR;
            localARrateR = 0;
            GammaO =  GammaO+(localARrateO/(frequency+1)-ratioAR)*GammaO;
            localARrateO = 0;
    end;
        

    %% Sample the LamdbaR and LambdaO chains
    aux = RndGamma(nn,:)./[StatRcurrent StatOcurrent];
    LambdaRcurrent = aux(1);
    LambdaOcurrent = aux(2);
    StatLRcurrent = T*log(LambdaRcurrent);
    StatLOcurrent = T*log(LambdaOcurrent);

    %% store the chain
    StoreRchain(:,nn+1) = Rcurrent;
    StoreOchain(:,nn+1) = Ocurrent;
    StoreLRchain(1,nn+1) = LambdaRcurrent;
    StoreLOchain(1,nn+1) = LambdaOcurrent;
    

end




%%%%%%%%%%
%% Output
%%%%%%%%%%
% The Gamma parameters 
output.GammaTildeR = GammaTildeR;
output.GammaO = GammaO;


% The empirical mean of the R, O, LamdbaR, LambdaO chains (burnin phase, discarded)
auxRchain = StoreRchain(:,burnin+1:NbrMC+1);
auxOchain = StoreOchain(:,burnin+1:NbrMC+1);
auxLRchain = StoreLRchain(1,burnin+1:NbrMC+1);
auxLOchain = StoreLOchain(1,burnin+1:NbrMC+1);
output.empirical_meanR = mean(auxRchain,2);
output.empirical_meanO = mean(auxOchain,2);
output.empirical_meanLR = mean(auxLRchain,2);
output.empirical_meanLO = mean(auxLOchain,2);

% Some quantiles, for each component of the R and O chains (burnin phase, discarded)
MatrixQuantileR = zeros(length(vectQ),T);   
MatrixQuantileO = zeros(length(vectQ),T);
for tt=1:T
    MatrixQuantileR(:,tt) = quantile(auxRchain(tt,:),vectQ');  % length(vectQ) x 1
    MatrixQuantileO(:,tt) = quantile(auxOchain(tt,:),vectQ');  % length(vectQ) x 1
end
output.quantilesR = MatrixQuantileR;
output.quantilesO = MatrixQuantileO;


% The chains LambdaR and LambdaO, after discarding the burnin period
output.Lambdachain(1,:) = auxLRchain;
output.Lambdachain(2,:) = auxLOchain;
