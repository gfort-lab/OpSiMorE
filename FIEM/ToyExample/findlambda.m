function out = findlambda(Smem,S,n,auxPdt,Pi2)

%% the current value of Smemory, dim_theta x 1
%% the current value of the FIEM-statistic, dim_theta x 1
%% the number of examples, n
%% auxPdt, Pi2 : parameters


    bars = auxPdt+Pi2*S;  % dim_theta x n
    auxmean = mean(Smem,2);
    NumCoeff = mean(sum(bars.*(auxmean-Smem),1));
    DenomCoeff = mean(sum(Smem.^2,1))-sum(auxmean.^2);

     
out = -NumCoeff/DenomCoeff;