%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 	SAMPLE THE DATA
%%
%% Define A, X, theta_true
%% Sample Z, Y
%%
%% save Y, A, X and theta_true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


%% number of examples 
n = input('Enter the size of the data set (the default value is 1e3):\n');
if isempty(n)==1
   n = 1e3;
end;

%% dimensions of the models
dim_theta = input('Enter the size of the parameter theta (default value: 20):\n');
if isempty(dim_theta)==1,
    dim_theta = 20;
end;

dim_Y = input('Enter the size of each observation Y (default value: 15):\n');
if isempty(dim_Y)==1,
    dim_Y = 15;
end;

dim_Z = input('Enter the size of each latent variable Z (default value: 10):\n');
if isempty(dim_Z)==1,
    dim_Z = 10;
end;

dim_A = [dim_Y,dim_Z];
dim_X = [dim_Z,dim_theta];


disp('For creating the data set, the true parameter theta is sparse (40 \% of the components are zero)');
%% theta_true, with 40% of the components fixed to zero
theta_true = 10*rand(dim_theta,1)-5;    % dim_theta x 1
erase = randsample(dim_theta,floor(0.4*dim_theta));
theta_true(erase) = zeros(1,length(erase));

%% the matrix A, as an AR(1)
rho = 0.8;
A = zeros(dim_A);
A(:,1) = sqrt(1-rho^2)*randn(dim_Y,1);

for jj=2:dim_Z,
    A(:,jj) =  rho*A(:,jj-1)+randn(dim_Y,1);  % dim_Y x 1
end;
disp(['The matrix A is obtained as a stationary AR(1) with variance 1 and coefficients ',num2str(rho)]);



%% the matrix X
rho = 0.9; 
X = zeros(dim_X);
X(:,1) = sqrt(1-rho^2)*randn(dim_Z,1);

for jj=2:dim_theta,
    X(:,jj) = sqrt(1-rho^2)*randn(dim_Z,1);   % 1 x dim_Z
end;
disp(['The matrix X is obtained as a stationary AR(1) with variance 1 and coefficients ',num2str(rho)]);


%% Sample the Zi
Zmatrix = X*theta_true + randn(dim_Z,n);    % dim_Z x n

%% Sample the Yi
Ymatrix = A*Zmatrix + randn(dim_Y,n);   % dim_Y  x n


%% Save the data
save Data.mat Ymatrix X A theta_true;

clear X A Zmatrix;
