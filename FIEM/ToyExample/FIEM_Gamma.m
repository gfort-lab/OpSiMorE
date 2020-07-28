%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fast Incremental EM 
%% the optimal one (see section 2.3.4) in the paper HAL-03617725
%%
%% Codes by G. Fort, May 2020
%% paper "Fast Incremental Expectation Maximization for non-convex finite-sum optimization: 
%% non asymptotic convergence bounds", HAL-03617725
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long e

tic 
fprintf('\t \t *** opt-FIEM *** \n');
fprintf('\t when nothing happens: check the windows "MENU" to answer the questions and fix the values of the design parameters\n');


DisplayBool = menu('During the run, do you want to display some graphical controls ?','No','Yes');
DisplayBool = DisplayBool-1;    % 1 for display some plots; 0 otherwise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Definition of the model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the examples
fprintf('Load the data set (obtained by running SampleData.m) \n');
load Data.mat   % Ymatrix, A, X

%% Dimensions of the model
[dim_Y,n] = size(Ymatrix);
[dim_Z,dim_theta] = size(X);
fprintf('The sample size is n = %f \n', n);


%% Penalty term   
upsilon = input('Enter the regularization parameter upsilon (the default value is 0.1):\n');
if isempty(upsilon)==1
   upsilon = 0.1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary quantities for the Gaussian model
%% often called in the algorithm ---> computed here
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some inverse matrices
auxinv1 = inv(eye(dim_Z)+A'*A);
auxinv2 = inv(upsilon*eye(dim_theta)+X'*X);

% matrices Pi1, Pi2
Pi1 = X'*auxinv1*A';
Pi2 = X'*auxinv1*X*auxinv2;

% vmin, vmax, L, Lvdot
auxeig = eig(X'*X);
vmin = 1/(upsilon+max(auxeig));
vmax = 1/(upsilon+min(auxeig));
L = sqrt(max(eig(Pi2'*Pi2))); 
Lvdot = max(abs(eig(auxinv2*(Pi2-eye(dim_theta)))));

% Mean value of the observations
barY = mean(Ymatrix,2);

% Fixed parts of the update i-EM scheme 
auxEM = Pi1*barY;   % dim_theta x 1
auxPdt = Pi1*Ymatrix; % dim_theta x n

clear auxeig 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% True optimum (unique) : theta_star
%%
%% computed for comparison, to assess the results
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxinv3 = inv(eye(dim_Y)+A*A');
theta_star = inv(upsilon*eye(dim_theta)+X'*A'*auxinv3*A*X)*X'*A'*auxinv3*barY;  % dim_theta x 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Iterative algorithm FIEM
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nbr of independent runs, for boxplot
NbrMC = input('Enter the number of independent paths of FIEM (default value: 1e3):\n');
if isempty(NbrMC)==1
   NbrMC = 1e3;
end;

%% Nbr of iterations per FIEM run 
Kmax = input('Enter the number of iterations per path, Kmax (default value: 20 n):\n');
if isempty(Kmax)==1,
    Kmax = n*20;
end;


%% Definition of the learning rate gamma
fprintf('Definition of the learning rate \n');
WhichGamma = menu('Learning rate as in Fort et al. (recom: first choice)','(recommended) rate n^(2/3)','rate n^(1/2)');
% fix the values of mu and lambda
mu = input(sprintf('\t \t Value of mu in (0,1) - default value: 0.25:\n'));
ll = input(sprintf('\t \t Value of lambda (default value: 0.5):\n'));
if isempty(mu)==1,
    mu = 0.25;
end;
if isempty(ll)==1,
    ll = 0.5;
end;
if WhichGamma==1,
    Cstar_1 = findcstar_1((L*vmin/Lvdot)^2,n,ll,2*mu*vmin*L/Lvdot);
    gamma_gfm = sqrt(Cstar_1)/(n^(2/3)*L);
else 
    Cstar_2 = findcstar_2((L*vmin/Lvdot)^2,n,ll,Kmax,2*mu*vmin*L/Lvdot);
    gamma_gfm = sqrt(Cstar_2)/(n^(1/3)*Kmax^(1/3)*L);
end;
gamma_grid = gamma_gfm*ones(1,Kmax);


%% Initialization of the algorithm 
fprintf('Initialization of the FIEM paths: \n');
fprintf('\t the first time, choose at random; it will store the answer in a file \n');
fprintf('\t then, you can choose to read it by loading this file \n');
NewInit = menu('How to initialize the FIEM path','At random','Read it in a file');
if NewInit == 1,
    Sinit = randn(dim_theta,1); % dim_theta x 1
    save InitStat.mat Sinit
else
    load InitStat;
end;

%% compute T(S) when S = Sinit
mapSinit = auxinv2*Sinit;   % dim_theta x 1

%% Sample all the indices of the examples selected during the path
% First, the indices for the auxiliary quantities
fprintf('Selection of the examples called along the FIEM path: \n');
fprintf('\t when updating the auxiliary quantity \n');
fprintf('\t \t the first time, choose at random; it will store the answer in a file \n');
fprintf('\t \t then, you can choose to use the same by loading this file \n');
NewSample = menu('Selection of the examples along the path','At random','Read it in a file');
if NewSample == 1,
    RandomIndexImatrix =  ceil(n*rand(NbrMC,Kmax));
    save RandomIndex.mat RandomIndexImatrix
else
    load RandomIndex.mat
end;

% Then, the indices for the update mecanism of the statistics 
fprintf('\t when updating the statistics \n');
fprintf('\t \t the first time, choose at random; it will store the answer in a file \n');
fprintf('\t \t then, you can choose to read it by loading this file \n');
NewSampleFIEM = menu('Selection of the examples along the path','At random','Read it in a file');
if NewSampleFIEM == 1,
    RandomIndexJmatrix =  ceil(n*rand(NbrMC,Kmax));
    RandomIndexItildematrix = ceil(n*rand(NbrMC,Kmax));
    save RandomIndexFIEM.mat RandomIndexJmatrix 
else
    load RandomIndexFIEM.mat
end;




%% Loop on the independent runs
% Store the norm of the field H
FieldH = zeros(NbrMC,Kmax);
ExpFieldH = zeros(NbrMC,Kmax);
StoreCoeff = zeros(NbrMC,Kmax);

for nn=1:NbrMC,
    % display the number of independent runs
    fprintf('FIEM run, number: %f \n',nn);

    % store : a matrix with the successive Kmax values of the statistics
    % (size dim_theta), 
    S_FIEM = zeros(dim_theta,Kmax);  % (dim_theta) x Kmax
    % initializaton of the FIEM sequence
    S_FIEM(:,1) = Sinit; % dim_theta x 1
    % initialization of the  Smemory vector
    s1aux = auxPdt+(Pi2*Sinit)*ones(1,n);  % dim_theta x n
    % Store as many Smemory as stepsize sequences (number: 1)
    Smem_FIEM = s1aux;    % dim_theta x n
    % initialization of the recursive computation of the mean of Smemory
    s2aux = mean(s1aux,2); % dim_theta x 1
    % Store as many means as stepsize sequences (number: NbrGamma)
    auxmeanSmem_FIEM = s2aux; % dim_theta x 1

    
    %% FIEM loop
    %%%%%%%%%%%%
    % choice of the random indices, 
    RandomIndexI = RandomIndexImatrix(nn,:);
    RandomIndexJ = RandomIndexJmatrix(nn,:);
   
    for k = 1:(Kmax-1),
        % select the index I and J 
        I = RandomIndexI(k+1);
        J = RandomIndexJ(k+1);
        % Update component I of Smemory
        prevSmem = Smem_FIEM(:,I);   % dim_theta x 1
        auxiter = Pi2*S_FIEM(:,k);
        Smem_FIEM(:,I) = auxPdt(:,I)+auxiter; % dim_theta x 1
        % Update the mean value of Smemory
        auxmeanSmem_FIEM = auxmeanSmem_FIEM+(Smem_FIEM(:,I)-prevSmem)/n;   % (dim_theta) x 1
        % Update the statistics 
            % compute a multiplicative coefficient for the control variate
        Coeff = findlambda(Smem_FIEM,S_FIEM(:,k),n,auxPdt,Pi2);
            %Coeff = 1;     % Original FIEM
            %Coeff = 0;     % No control variate = Online EM 
        StoreCoeff(nn,k+1) = Coeff;
            % store the field H (to compare the strategies)
        auxField = auxPdt(:,J)+auxiter-S_FIEM(:,k)+Coeff*(auxmeanSmem_FIEM-Smem_FIEM(:,J)); % (dim_theta) x 1
        %[norm(auxPdt(:,J)+auxiter-S_FIEM(:,k)) norm((auxmeanSmem_FIEM-Smem_FIEM(:,J))) norm(auxField)]
        FieldH(nn,k+1) = norm(auxField)^2; % 1 x 1
        ExpFieldH(nn,k+1) = norm(auxEM+auxiter-S_FIEM(:,k))^2;  % 1 x 1
            % compute the statistis
        S_FIEM(:,k+1) = S_FIEM(:,k)+gamma_grid(k+1).*auxField;    % (dim_theta) x 1
    end; % FIEM loop
   
  
        % compute the theta-path 
        Theta_FIEM = auxinv2*S_FIEM;  % (dim_theta*) x Kmax
        % Compute the distance to theta_star
        DeltaTheta_FIEM = sqrt(sum((Theta_FIEM-theta_star).^2,1));  % 1 x Kmax
        
    %% DISPLAY
    if DisplayBool == 1,
        CntFig = 1;
        % Display on the same plot
        figure(CntFig);
        clf;
        burnin = 2;
        stop = Kmax;
        semilogy(burnin:stop,DeltaTheta_FIEM(1,burnin:stop)','Linewidth',2);
        title('FIEM theta-path : distance to the unique optimum');
        % Display the burn-in phase
        CntFig = CntFig+1;
        figure(CntFig);
        clf;
        burnin = 1;
        stop = 100;
        semilogy(burnin:stop,DeltaTheta_FIEM(1,burnin:stop)','Linewidth',2);
        title('theta-path : distance to theta_* (burn in)');
    end;
    
    % Store some boxplots
    GridIter = [1:1000:Kmax];
        % GridIter = [4000 6000 8000 10000 12000];
    GridGamma = 1;
    ErrorTheta(nn,:) = DeltaTheta_FIEM(1,GridIter);
    
    
  if DisplayBool == 1,
      CntFig = CntFig+1;
      figure(CntFig);
      clf;
      plot(2:Kmax,StoreCoeff(nn,2:end));
      hold on;
      plot([1 Kmax],[1,1],'r--');
      title('The optimal leverage coeff \lambda_k (see section 2.3.4.)'); 
      
      CntFig = CntFig+1;
      figure(CntFig);
      clf;
      semilogy(1:Kmax,FieldH(nn,:));
      title('Squared norm of the perturbed field H_k vs nbr iterations'); 
      
      CntFig = CntFig+1;
      figure(CntFig);
      clf;
      plot(1:Kmax,S_FIEM);
      title('The FIEM statistics S_k vs nbr iterations k'); 
      
  end;
    
    
    
    save StoreCoeffopt.mat StoreCoeff FieldH ExpFieldH ErrorTheta; 
end;    % Indep runs loop

toc

