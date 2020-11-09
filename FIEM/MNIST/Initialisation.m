%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%  INITIALIZATION OF EM
%%
%% based on the paper 
%%  https://link.springer.com/content/pdf/10.1007/s10044-014-0441-3.pdf
%% by W. Kwedlo, Pattern Anal Applic (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VectMean InvCovMatrix] = Initialisation(X,L)
    [d n] = size(X);

    TrCov = trace(X*X');
    
    % Nbr of trials
    T = 20; % as advocated in the paper

    %% Choose one vector at random
    I = ceil(n*rand(1));
    %% Initialize the vector of means
    VectMean = zeros(d,L);
    VectMean(:,1) = X(:,I);
    %% Initialize the covariance matrix
    InvCovMatrix = zeros(d,d*L);
    InvCovMatrix(:,1:d) = SampleInvCov(d,TrCov,L);

    %% Loop in order to determine the other means
    OldSet = 1:n;
    NewSet = setdiff(OldSet,I);

    for mm = 2:L,
        % sample T vectors of features
        I = randsample(NewSet,T);

        % compute the distance of these selected samples to the Means
        for tt=1:T,
            x = X(:,I(tt));
            for mmm = 1:mm-1,
                InvCov = InvCovMatrix(:,(mmm-1)*d+1:mmm*d);
                dist2(tt,mmm) = (x-VectMean(:,mmm))'*InvCov*(x-VectMean(:,mmm));
            end;
            [value(tt) location(tt)]= min(dist2(tt,:));
        end;
       % for each vector, find the closest mean
       [vv indice] = max(diag(dist2(:,location)));
       % Store the New Mean
       VectMean(:,mm) = X(:,I(indice));

        % Store the New Cov
        InvCovMatrix(:,(mm-1)*d+1:mm*d) = SampleInvCov(d,TrCov,L);

        clear OldSet
        OldSet  = NewSet;
        clear NewSet
        NewSet = setdiff(OldSet,I(indice));
    end;


end



