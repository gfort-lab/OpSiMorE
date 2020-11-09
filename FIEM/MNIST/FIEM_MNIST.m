%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FAST INCREMENTAL EM  
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

format long e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define binary variables to include some procedures or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunInit = 1;    % construct (0) or read in a file (1) an initialization vector of parameter.
DisplayPlot = 1; % display(1) or not(0) some plots at each iteration to control the behavior of the algo.



%%%%%%%%%%%%%%%%%%%%%
%% Load the data set
%%%%%%%%%%%%%%%%%%%%%
load Data.mat  % load Xred
[dinit n] = size(Xred);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix some quantities for the run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Number of selected features for each example
d = 20;
X = Xred(1:d,:);    % d x n
%% Number of classes
L = 12;
%% Nbr of epochs 
NbrEpoch = 53;
%% Size of the mini-batch
minibatch = 100;
if mod(n/minibatch,10)~=0,
    % check the minibatch size to be compatible with the implementation of
    % EvalDenGauss
    fprintf('Choose another size of the minibatch');
    return;
end;
%% Nbr of iterations (total nbr)
NbrIter = NbrEpoch*n/minibatch;
%% Nbr independent runs of the algorithm
NbrRun = 10;
%% Start like an 'online-EM' until  kswitch and then switch to 'FIEM'
kswitch = 6;
vectok = [zeros(1,kswitch) ones(1,NbrEpoch-kswitch)];
%% Constant step size sequence
vectpas = [5e-3*ones(1,kswitch) 5e-3*ones(1,NbrEpoch-kswitch)];

%% Compute the empirical cov matrix
COV = zeros(d,d,n);
for nn = 1:n,
    COV(:,:,nn) = X(:,nn)*X(:,nn)';
end;
EmpiricalCov = sum(COV,3)/n;
clear COV;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('global INITIALIZATION \n');
fprintf(' \t create or load the initial parameters \n');
if RunInit == 1,
   [VectMean InvCovMatrix] = Initialisation(X,L); 
   Weight = 1/L*ones(1,L); % 1 x L
   save InitParameter.mat Weight VectMean InvCovMatrix
else
    load InitParameter.mat
    % load : "Weight" (1 x L); "VectMean" (d x L); "InvCovMatrix" (d x
    % (dL)), the L cov matrices are concatened.
end;
Weight = Weight';   % L x 1
auxIC = InvCovMatrix(:,1:d);
clear InvCovMatrix
InvCovMatrix = auxIC;
clear auxIC



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save space for storing quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\t Save space for some storage \n');
% The Memory matrix for the a posteriori distribution of each class
AposterioriMemory = zeros(L,n); % L x n
% The mean Memory vector for the weight
MeanWeightMemory = zeros(L,1);  % L x 1;
% The mean Memory matrix for the expectations
MeanMuMemory = zeros(d,L);  % d x L

% Store the parameters along the path
Sweight = zeros(L,NbrIter);  % (weight)
Smu = zeros(d,NbrIter);  % (one expectation among L)
Seigen = zeros(d,NbrIter);    % (eigenvalues of the InvCov matrix)

% Store the successive parameters, at each epoch (+ the initialization)
StoreWeight = zeros(L,NbrEpoch+1);  % (weight)
StoreMean1 = zeros(d,NbrEpoch+1);   % (one expectation among L)
StoreEigen = zeros(d,NbrEpoch+1);  % (eigenvalues of one InvCov matrix)
LogLikelihood = zeros(1,NbrEpoch+1);

% Store the parameters at the end of each run
RunStoreWeight = zeros(L,NbrRun);
RunStoreMean = zeros(d,L,NbrRun);
RunStoreInvCov = zeros(d,d,NbrRun);
RunLogLikelihood = zeros(NbrRun,NbrEpoch+1); % store also at initialization




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of the statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the Stat-weight 
    Aux = zeros(L,n);   % L x n
    for ll=1:L,
        eee = EvalDensGauss(X,VectMean(:,ll),InvCovMatrix,n); % 1 x n
        Aux(ll,:) = Weight(ll)*eee;
    end;
    ss = sum(Aux);  % 1 x n
    
fprintf('\t Initialize the likelihood \n');
    LogLikelihood(1) = mean(log(ss));
    RunLogLikelihood(:,1) = ones(NbrRun,1)*LogLikelihood(1);
    fprintf('\t \t the likelihood at init: %f \n',LogLikelihood(1));
   
fprintf('\t Initialize the statistics \n');
    fprintf('\t \t the weights \n');
    Proba = Aux./ss;  % L x n
    StatWeight = mean(Proba,2);  % L x 1
    
    % The a posteriori Memory matrix
    fprintf('\t \t The A posteriori Memory matrix \n');
    AposterioriMemory = Proba;  % L x n
    
    % The mean Memory for the weights
    fprintf('\t \t The mean Memory (weight) \n');
    MeanWeightMemory = StatWeight;  % L x n
    
    % Initialize the Stat-Mu
    fprintf('\t \t the expectation \n');
    StatMu = X*Proba'/n;    % d x L
    
    % The mean Memory for the expectations
    fprintf('\t \t The mean Memory (expectation) \n');
    MeanMuMemory = StatMu;  % d x L
  
fprintf('\t Store the initialization, common to all the indep run \n');
    fprintf('\t \t the parameters \n');
    InitWeight = Weight;    % L x 1
    InitMu = VectMean;  % d x L
    InitInvCov = InvCovMatrix;  % d x d
   
    fprintf('\t \t the Memory quantities \n');
    InitAposterioriMemory = AposterioriMemory;
    InitMeanWeightMemory = MeanWeightMemory;   
    InitMeanMuMemory = MeanMuMemory;
     
    fprintf('\t \t the statistics \n');
    InitStatWeight = StatWeight;
    InitStatMu = StatMu;
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Independent runs of FiEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for rr = 1:NbrRun,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Ech run is defined as two nested loops
    %% the outer one ---> which is a block of iterations, with a total cost of n conditional expectation evaluations
    %% the inner ones ---> which is a set of n/minibatch iterations for FIEM and OnlineEM.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('RUN nbr %f \n',rr);  
    
    %% Sample the examples 
    fprintf('\t Sample the examples, with replacement \n');
        %% For the auxiliary variables Memory
        SelectExamples = ceil(n*rand(1,NbrIter*minibatch));   % 1 x (NbrIter x minibatch)
        %% For the Stochastic Approximation step
        SelectExamplesBis = ceil(n*rand(1,NbrIter*minibatch));  % 1 x (NbrIter x minibatch)
      
    fprintf('\t Initialize the Memory quantities \n');
    AposterioriMemory = InitAposterioriMemory;
    MeanWeightMemory = InitMeanWeightMemory;
    MeanMuMemory = InitMeanMuMemory;
    
    fprintf('\t Initialize the statistics \n');
    StatWeight = InitStatWeight;
    StatMu = InitStatMu;
    
    fprintf('\t Initialize the parameters \n');
    Weight = InitWeight;    % L x 1
    VectMean = InitMu;  % d x L
    InvCovMatrix = InitInvCov;  % d x L
    
     %% Store the value of some quantities at the beginning of the run
     StoreWeight(:,1) = Weight;   % L x 1
     StoreMean1(:,1) = VectMean(:,1); % d x 1
     StoreEigen(:,1) = sort(eig(InvCovMatrix)); % d x 1

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% LOOP over the epochs 
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for kout = 1:NbrEpoch,
        %% start as an 'onlineEM' and then switch to 'FIEM'
        ok = vectok(kout);
        %% step size
        pas = vectpas(kout);
        
        fprintf('\t Epoch: %f and OnlineEM %f or FIEM %f \n ',kout,1-ok,ok);
        
        if kout == kswitch,
            fprintf('\t \t Start FIEM (initialize the Memory) \n');
                Aux = zeros(L,n);   % L x n
                for ll=1:L,
                    eee = EvalDensGauss(X,VectMean(:,ll),InvCovMatrix,n); % 1 x n
                    Aux(ll,:) = Weight(ll)*eee;
                end;
                ss = sum(Aux);  % 1 x n
                AposterioriMemory = Aux./ss;  % L x n
                MeanWeightMemory = mean(AposterioriMemory,2);  % L x n
                MeanMuMemory = (X*AposterioriMemory')/n;
        end;    % end of initialization of "Memory" when FIEM starts
        
        
        %% Select the examples for the current epoch
        BlockIndexEx = SelectExamples((kout-1)*n+1:(kout*n)); % 1 x n
        BlockIndexSA = SelectExamplesBis((kout-1)*n+1:(kout*n)); % 1 x n

        % Loop over the inner ones
        for kin = 1:(n/minibatch),  
            %% Select the examples
            bloc = (kin-1)*minibatch+1:kin*minibatch;
            J = BlockIndexEx(bloc);
            I = BlockIndexSA(bloc); 
            clear bloc
            
            %% Compute the a posteriori distribution
            %% and also for the SA-step
            %% for the examples I and J
            auxW = zeros(L,minibatch);
            auxWSA = zeros(L,minibatch);
            for ll=1:L,
                eee = EvalDensGauss(X(:,I),VectMean(:,ll),InvCovMatrix,minibatch); % 1 x minibatch
                auxW(ll,:) = eee*Weight(ll);  % auxW : L x minibatch
                % prepare the SA-step
                eeeSA = EvalDensGauss(X(:,J),VectMean(:,ll),InvCovMatrix,minibatch); % 1 x minibatch
                auxWSA(ll,:) = eeeSA*Weight(ll);  % auxWSA : L x minibatch
             end;
            Aposteriori = auxW./sum(auxW);  % L x minibatch
            AposterioriSA = auxWSA./sum(auxWSA);  % L x minibatch
            clear auxW auxWSA eee eeeSA
            
            % compute the "old" and "new" memory statistics
            auxOld = AposterioriMemory(:,I);   % L x minibatch
            auxNew = Aposteriori;   % L x minibatch
           
            OldWeightMemory = auxOld;  % L x minibatch
            
            OldMuMemory = zeros(d,L);
            NewMuMemory = zeros(d,L);
             for ll=1:L,
                OldMuMemory(:,ll) = sum(X(:,I)*diag(auxOld(ll,:)),2);    % d x 1
                NewMuMemory(:,ll) = sum(X(:,I)*diag(auxNew(ll,:)),2);    % d x 1
             end;
            clear bloc auxOld auxNew
           
            % Update the Memory quantities
            AposterioriMemory(:,I) = Aposteriori;   % L x minibatch
            MeanWeightMemory = MeanWeightMemory+(sum(Aposteriori,2)-sum(OldWeightMemory,2))/n;  % L x 1
            MeanMuMemory = MeanMuMemory+(NewMuMemory-OldMuMemory)/n;   % d x L
            clear OldWeightMemory NewMuMemory OldMuMemory  bloc auxOld auxNew
           
            % Update the statistics : weight
            barsWeight = mean(AposterioriSA,2);  % L x 1
            TempStatWeight = StatWeight + pas*(barsWeight-StatWeight);  % L x 1 
            StatWeight = TempStatWeight + pas*ok*(MeanWeightMemory-mean(AposterioriMemory(:,J),2));   % L x 1
            if min(StatWeight)<0,
                fprintf('WEIGHT : epoch %f, iter %f\n',kout,kin);
                keyboard
            end;
            Weight = StatWeight/sum(StatWeight);    % L x 1
            
            % Update the statistics : expectations
            barsMu = zeros(d,L);
            barsMuAux = zeros(d,L);
              for ll=1:L,
                barsMu(:,ll) = mean(X(:,J)*diag(AposterioriSA(ll,:)),2); % d x 1
                barsMuAux(:,ll) = mean(X(:,J)*diag(AposterioriMemory(ll,J)),2); % d x 1
              end;
            TempStatMu = StatMu + pas*(barsMu-StatMu);  % d x L
            StatMu = TempStatMu + pas*ok*(MeanMuMemory-barsMuAux);  % d x L
            clear barsMu barsMuAux TempStatMu TempStatWeight
            
            % Update the parameters 
            VectMean = StatMu*diag(1./StatWeight);  %   d x L
            recursive = zeros(d,d);
            for ll=1:L
                recursive = recursive + StatWeight(ll)*VectMean(:,ll)*VectMean(:,ll)';
             end;
             CovMatrix = EmpiricalCov-recursive;    % d x d
             [vect value] = eig(CovMatrix);
             if (min(diag(value))<0),
                   fprintf('epoch %f, iter %f,  pas %f \n',kout,kin,ok);
                   keyboard
             end;
             InvCovMatrix = inv(CovMatrix);    % d x d
            
            %% STORE the weights and the first expectation, and some eigenvalues
            Sweight(:,(kout-1)*n/minibatch+kin) = Weight;  % L x 1
            Smu(:,(kout-1)*n/minibatch+kin) = VectMean(:,1);    % d x 1
            Seigen(:,(kout-1)*n/minibatch+kin) = sort(eig(InvCovMatrix));    % d x 1
        end;    % END of "kin"

        
        
         %% Store the value of some parameters at the end of every epoch
         StoreWeight(:,kout+1) = Weight;   % L x 1
         StoreMean1(:,kout+1) = VectMean(:,1); % d x 1
         StoreEigen(:,kout+1) = sort(eig(InvCovMatrix)); % d x 1

          %% Compute the LogLilelihood at the end of every epoch
          Aux = zeros(L,n);  
          for ll=1:L,
            eee = EvalDensGauss(X,VectMean(:,ll),InvCovMatrix,n); % 1 x n
            Aux(ll,:) = Weight(ll)*eee; % 1 x n
          end;
          LogLikelihood(kout+1) = mean(log(sum(Aux)));
          fprintf('\t \t  LogLikelihood %f \n',LogLikelihood(kout+1));
          RunLogLikelihood(rr,kout+1) = LogLikelihood(kout+1);
          
          
         %% CONTROLS : display the evolution of some quantities
             if DisplayPlot==1,
                    % the weight
                    figure(1);
                    clf
                    semilogy(1:n*kout/minibatch,Sweight(:,1:kout*n/minibatch)');

                    % a vector of mean of the 1st component
                    figure(2)
                    clf;
                    plot(1:n*kout/minibatch,Smu(:,1:kout*n/minibatch)');

                    % the eignevalues of the 1st inv-cov matrix
                    figure(3);
                    clf;
                    semilogy(1:n*kout/minibatch,Seigen(:,1:kout*n/minibatch)'); 
                    
                    % the likelihood 
                    figure(4);
                    clf;
                    semilogy(1:kout+1,LogLikelihood(1:kout+1)); 
             end
    end;    % END of the loops : kout
    
    
    fprintf('\t End of the run : store some estimated parameters \n');
    RunStoreWeight(:,rr) = Weight;  % L x 1
    RunStoreMean(:,:,rr) = VectMean;    % d x L
    RunStoreInvCov(:,:,rr) = InvCovMatrix;  % d x d
    
    save RunEstimationFIEM.mat L RunStoreWeight RunStoreMean RunStoreInvCov RunLogLikelihood Smu Sweight Seigen StoreWeight StoreMean1 StoreEigen
 
end; %%END  of the loops : NbrRun


function out = EvalDensGauss1(Y,M,IC)
    dd = det(IC);
    out = exp(-0.5*(Y-M)'*IC*(Y-M)+0.5*log(dd));
end

function out = EvalDensGauss(Y,M,IC,n)
    bloc = n/10;
    dd = det(IC);
    for jj=1:10
        Yaux = Y(:,(jj-1)*bloc+1:jj*bloc);
        out((jj-1)*bloc+1:jj*bloc) = exp(diag(-0.5*(Yaux-M)'*IC*(Yaux-M))+0.5*log(dd));
    end;
end