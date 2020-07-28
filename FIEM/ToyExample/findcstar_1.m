%%%% Find Cstar
%% n : is the batch size
%% ll : is a coefficient  in (0,1)
%% A : defines the quantity to be reached
%% Cinit : an initialization
function out = findcstar_1(Cinit,n,ll,A)

Cup = ll*n^(1/3);
Cdown = 0; 
C = Cinit;

dist = A - sqrt(C)*(n^(-2/3)+(C*(ll-C*n^(-1/3))^(-1))*(1/n+1/(1-ll)));
while abs(dist)>=1e-8
    if dist>0 
        Cdown = C;
        C = (C+Cup)/2;
    else
        Cup = C;
        C = (C+Cdown)/2;
    end;
    dist = A - sqrt(C)*(n^(-2/3)+(C*(ll-C*n^(-1/3))^(-1))*(1/n+1/(1-ll)));
end;



out = C;
        
