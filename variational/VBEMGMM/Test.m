## Read data

id = fopen('faithful.txt','r')
x = fscanf(fid,'%f',[2,272]);
fclose(fid);
train_data = x;

## Standardize the data

train_data = train_data - repmat(mean(train_data,2),1,size(train_data,2));
train_data = train_data./repmat(var(train_data,[],2),1,size(train_data,2));
[dim N] = size(train_data);

## Intialize the priors

ncentres = 15
PriorPar.alpha = .001;
PriorPar.mu = zeros(dim,1);
PriorPar.beta = 1;
PriorPar.W = 200*eye(dim);
PriorPar.v = 20;

## The first iteration

x = train_data;
[D N] = size(x);

## Initialize variables

K = ncentres;
alpha0 = PriorPar.alpha;
m0 = PriorPar.mu;
beta0 = PriorPar.beta;
W0 = PriorPar.W;



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
