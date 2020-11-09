function IC = SampleInvCov(d,TrCov,L)
        %% choose the eigenvalues 
            % pick at random d positive values in [0, A]
        VectEig = sort(TrCov*rand(d,1)/(10*d*L));
            % the values lower than max(eig)/10 are set to max(eig)/10
        VectEig = max(VectEig,VectEig(end)/10);
            % the sum of the eig is equal to A
        VectEig = VectEig/(sum(VectEig))*TrCov/(10*d*L);

        %% Obtain a d x d unitary matrix
        [Q  R] = qr(randn(d,d));

        %% Return the inverse of a covariance matrix
        IC = Q*diag(1./VectEig)*Q';
end