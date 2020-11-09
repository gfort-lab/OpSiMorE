%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD and VISUALIZE the MNIST data base
%%
%% 70 000 images of size 28x28 = 784
%% the first 60000 are the training set
%% the last 10000 are the test set
%%
%% the entries of the images are in the range [0, ..., 255]
%%
%% the mat file was downloaded from
%% https://www.kaggle.com/avnishnish/mnist-original#mnist-original.mat
%%
%% the description of the data is 
%% http://yann.lecun.com/exdb/mnist/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
CntFig = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD
%% this provides : 
%% data : size 784 x 70 000;
%% label : size 1 x 70 000

load mnist-original.mat;

%% TRAINING and TEST sets
data_train = double(data(:,1:60000)');
[n dinit] = size(data_train);
data_test = double(data(:,60001:70000)');
clear data;

label_train = label(1:60000);   % 
label_test = label(60001:70000);
clear label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some pixels are always zero; find them and remove them
%% there are 67 such pixels
minimal = min(data_train);
maximal = max(data_train);
keep = find(minimal~=maximal);
X = data_train(:,keep); % n x d
[n d] = size(X);   % d = 717 
clear minimal maximal keep


%%%%%%%%%%%%%%%%%%%%
%% DISPLAY AN IMAGE
%%%%%%%%%%%%%%%%%%%
% a picture of a "3"
nn = min(find(label_train==3));
    figure(CntFig);
    CntFig = CntFig+1;
    clf;
    imagesc(reshape(data_train(nn,:),28,28)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reduce the dimension by an ACP
%% center and normalize the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% center each column
Xcent = X - mean(X);
% standardize  each column
Xstd = Xcent/diag(std(Xcent));

% compute the svd of Xstd'*Xstd
[P Delta Q] = svd(Xstd'*Xstd);   % matrice = P Delta Q';
% compute the eigenvalues
eigen = diag(Delta); 

    % plot the eigenvalues
    figure(CntFig);
    clf;
    subplot(2,1,1);
    semilogy(eigen,'x-');
    subplot(2,1,2);
    semilogy(eigen(1:end-5),'x-');  % remove the lats 5 ones

%% Projection of the data along the principal components
Xred = Q'*Xstd';    % d x n

    % plot the n observations based on the first 3 directions
    CntFig=CntFig+1;
    figure(CntFig);
    clf;
    plot3(Xred(1,:),Xred(2,:),Xred(3,:),'.r');

    % distinguish per class
    vectcolor = ['.r','.g','.b','.y','.c','.k','.m','xr','xg','xb'];
    CntFig=CntFig+1;
    figure(CntFig);
    clf;
    for ll = [0:1:9],
        class = find(label_train==ll);
        plot3(Xred(1,class),Xred(2,class),Xred(3,class),[vectcolor(2*ll+1) vectcolor(2*ll+2)]);
        hold on;
    end;
    
    
    % distinguish per class
    vectcolor = ['.r','.g','.b','.y','.c','.k','.m','xr','xg','xb'];
    CntFig=CntFig+1;
    figure(CntFig);
    clf;
    for ll = [0:1:9],
        class = find(label_train==ll);
        plot(Xred(1,class),Xred(2,class),[vectcolor(2*ll+1) vectcolor(2*ll+2)]);
        hold on;
    end;
        
    
    % distinguish per class
    vectcolor = ['.r','.g','.b','.y','.c','.k','.m','xr','xg','xb'];
    CntFig=CntFig+1;
    figure(CntFig);
    clf;
    for ll = [0:1:9],
        class = find(label_train==ll);
        plot(class,Xred(1,class),[vectcolor(2*ll+1) vectcolor(2*ll+2)]);
        hold on;
    end;
    
save Data.mat Xred 