%%%% Find Cstar
%% n : is the batch size
%% Kmax : the maximal number of iterations
%% ll : is a coefficient  in (0,1)
%% A : defines the quantity to be reached
%% Cinit : an initialization

function out = findcstar_2(Cinit,Kmax,n,ll,A)

Cup = ll*Kmax^(2/3)*n^(-1/3);
Cdown = 0; 
C = Cinit;

dist = A - sqrt(C)*((n*Kmax)^(-1/3)+C*(1/n+1/(1-ll)));
while abs(dist)>=1e-8
    if dist>0 
        Cdown = C;
        C = (C+Cup)/2;
    else
        Cup = C;
        C = (C+Cdown)/2;
    end;
    dist = A - sqrt(C)*((n*Kmax)^(-1/3)+C*(1/n+1/(1-ll)));
end;


out = C;
        
