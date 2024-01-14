function [output] = SAEM_nomixture(data,SAEM,MCMC)
%% SAEM_nomixture
%%
%% Runs a SAEM algorithm (see Delyon, Lavielle, Moulines, 1999)
%% for solving a maximum likelihood in a latent variable model. 
%%
%% Developed by Gersende Fort, January 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read the input variables 

% Read the data Z and Phi 
Z = data.Z; % T x 1
Phi = data.Phi; % T x 1
T = size(Z,1);

% Number of SAEM iterations
if isfield(SAEM,'NbrIter')
    NbrIter  = SAEM.NbrIter; 
else 
    NbrIter  = 5e5; 
end


% Initial value of the SAEM sequence
% for LambdaR
if isfield(SAEM,'LambdaRinit')
    LambdaRcurrent  = SAEM.LambdaRinit; 
else 
    LambdaRcurrent  = 3.5*std(Z); 
end
% for LambdaO
if isfield(SAEM,'LambdaOinit')
    LambdaOcurrent  = SAEM.LambdaOinit; 
else 
    LambdaOcurrent  = 0.05; 
end


% Step size of SAEM
if isfield(SAEM,'pas_vect_R')
    pas_vect_R  = SAEM.pas_vect_R; 
else 
    pas_vect_R = 0.05*[ones(1,10) 2*ones(1,100) 4./sqrt(100:100+(NbrIter-110))];
end
if isfield(SAEM,'pas_vect_O')
    pas_vect_O  = SAEM.pas_vect_O; 
else 
    pas_vect_O = 0.5*[0.1*ones(1,10) 0.05*ones(1,100) 0.1./sqrt(100:100+(NbrIter-110))];
end

   
% Display some control or not
controldisplay = SAEM.controldisplay;


% Length of the Markov chain
if isfield(MCMC,'chain_length')
    NbrMC  = MCMC.chain_length; 
else 
    NbrMC  = [1e7 5e6  3e3*ones(1,NbrIter-2)];
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

% The Gamma parameters
if isfield(MCMC,'GammaTildeR')
    GammaTildeR = MCMC.GammaTildeR;    % 1 x 1
else 
    GammaTildeR = 1e-12;
end
if isfield(MCMC,'GammaO')
    GammaO = MCMC.GammaO;    % 1 x 1
else 
    GammaO = 1e3;
end


% Adaptation of GammaTildeR and GammaO
if isfield(MCMC,'adapt_frequency')
    frequency = MCMC.adapt_frequency;    % 1 x NbrIter
else 
    frequency = max(burnin/500,500);  % 1 x NbrIter
end
if isfield(MCMC,'target_ratioAR')
    ratioAR = MCMC.target_ratioAR;    % 1 x 1
else 
    ratioAR = 0.25; 
end


% Store some variables
StoreLambdaR = zeros(1,NbrIter+1);
StoreLambdaR(1,1) = LambdaRcurrent;
StoreLambdaO = zeros(1,NbrIter+1);
StoreLambdaO(1,1) = LambdaOcurrent;
StoreGammaTildeR = zeros(1,NbrIter+1);
StoreGammaTildeR(1,1) = GammaTildeR;
StoreGammaO = zeros(1,NbrIter+1);
StoreGammaO(1,1) = GammaO;



%% Initialize the SAEM algorithm
data.LambdaR = LambdaRcurrent;  % 1 x 1
data.LambdaO = LambdaOcurrent;  % 1 x 1

MCMClocal.chain_length = NbrMC(1);  % 1 x 1
MCMClocal.burnin = burnin(1);  % 1 x 1
MCMClocal.initial_pointR =Rcurrent;  % T x 1
MCMClocal.initial_pointO = Ocurrent;  % T x 1
MCMClocal.GammaTildeR = GammaTildeR;
MCMClocal.GammaO = GammaO; 
MCMClocal.adapt_frequency = frequency(1);  % 1 x 1
MCMClocal.target_ratioAR = MCMC.target_ratioAR; % 1 x 1
MCMClocal.Qvec = [];

output = GibbsPGdual_nomixture(data,MCMClocal);


% Initialize the SAEM statistics
SAEMStatR = output.StatR;
SAEMStatO = output.StatO; 


% update the parameters 
data.LambdaR = T/SAEMStatR;   
data.LambdaO = T/SAEMStatO; 
% Store the parameters
StoreLambdaR(1,2) = data.LambdaR;
StoreLambdaO(1,2) = data.LambdaO;


%% Iterate the SAEM algorithm
for nn=2:NbrIter
    % Define the input structures for the next sampling step
    MCMClocal.chain_length = NbrMC(nn);  % 1 x 1
    MCMClocal.burnin = burnin(nn);  % 1 x 1
    MCMClocal.adapt_frequency = frequency(nn);  % 1 x 1
    MCMClocal.initial_pointR = output.lastsampleR;  % T x 1
    MCMClocal.initial_pointO = output.lastsampleO;  % T x 1
    MCMClocal.GammaTildeR = output.GammaTildeR;  % 1 x 1
    MCMClocal.GammaO = output.GammaO;    % 1 x 1
    
    StoreGammaTildeR(1,nn) = output.GammaTildeR;    % 1 x 1
    StoreGammaO(1,nn) = output.GammaO;      % 1 x 1

    output = GibbsPGdual_nomixture(data,MCMClocal); 

    % Update the SAEM statistics
    SAEMStatR = SAEMStatR + pas_vect_R(nn)*(output.StatR-SAEMStatR);    % 1 x 1
    SAEMStatO = SAEMStatO + pas_vect_O(nn)*(output.StatO-SAEMStatO);    % 1 x 1
   
    % update the parameters 
    data.LambdaR = T/SAEMStatR;     % 1 x 1 
    data.LambdaO = T/SAEMStatO;     % 1 x 1
    % Store the parameters
    StoreLambdaR(1,nn+1) = data.LambdaR; % 1 x 1
    StoreLambdaO(1,nn+1) = data.LambdaO; % 1 x 1

    if ((controldisplay==1) && (mod(nn,5e2)==0))
        figure(50);
        clf
        plot(2:nn+1,StoreLambdaR(1,2:nn+1));
        grid on
        title('SAEM sequence: \lambda_R');

        figure(51);
        clf
        plot(2:nn+1,StoreLambdaO(1,2:nn+1));
        grid on
        title('SAEM sequence: \lambda_O');

    end

end


%% Define the output structure
output.LambdaRpath =  StoreLambdaR;
output.LambdaOpath =    StoreLambdaO;
output.GammaTildeRpath = StoreGammaTildeR;
output.GammaOpath = StoreGammaO;
end