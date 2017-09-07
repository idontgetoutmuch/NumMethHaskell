% This file illustrates the usage of the function 'gmmVBEM.m'
% Refer to Bishop's book for notation and details 
% @book{bishop2006pra,
%   title={{Pattern recognition and machine learning}},
%   author={Bishop, C.M.},
%   year={2006},
%   publisher={Springer}
% }
% This function needs faithful.txt, gmmVEBM.m, MyEllipse.m
% Written by Emtiyaz, CS, UBC 
% June, 2007

%read data
fid = fopen('faithful.txt','r');
x = fscanf(fid,'%f',[2,272]);
fclose(fid);
train_data = x;
%standardize the data
train_data = train_data - repmat(mean(train_data,2),1,size(train_data,2));
train_data = train_data./repmat(var(train_data,[],2),1,size(train_data,2));
[dim N] = size(train_data);

%%%%%%%%%%%%%%%%%%%%
% GMM VBEM clustering
%%%%%%%%%%%%%%%%%%%%
% initialize with EM algorithm 
ncentres = 15;
mix = gmm(2, ncentres, 'full'); 
options = foptions;
options(14) = 10;
mix = gmminit(mix, train_data', options);
maxIter = 30; options(3) = 0.1; options(14) = maxIter;
[mix, options, errlog] =  gmmem(mix, train_data', options);
% intialize the priors
PriorPar.alpha = .001;
PriorPar.mu = zeros(dim,1);
PriorPar.beta = 1;
PriorPar.W = 200*eye(dim);
PriorPar.v = 20;

% set the options for VBEM
clear options;
options.maxIter = 100;
options.threshold = 1e-5;
options.displayFig = 1;
options.displayIter = 1;

% Call the function
[out] = gmmVBEM(train_data, mix, PriorPar, options);
%plot lower bound
figure
plot(out.L)


