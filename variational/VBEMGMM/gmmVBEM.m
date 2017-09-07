function [out] = gmmVBEM3(x, mix, PriorPar, options)

% Variational Bayes EM algorithm for Gaussian Mixture Model
% This implementation is based on Bishop's Book
% Refer to Bishop's book for notation and details 
% @book{bishop2006pra,
%   title={{Pattern recognition and machine learning}},
%   author={Bishop, C.M.},
%   year={2006},
%   publisher={Springer}
% }
% Function uses inputs similar to Netlab
%%%%%%%%%%%%%%%%
% INPUT
% D is the dimension, N is the number of Data points
% x  - training data of size DxN 
% mix - gmm model initialize with netlab's gmmem function 
% PriorPar - structure containing priors
% options - options for maxIter, threshold etc. etc.
% FOR DETAILS OF ABOVE VARIABLE, REFER TO EXAMPLE FILE 
% out is the output vector (see the end of this file for details)
%%%%%%%%%%%%%%%%
% Written by Emtiyaz, CS, UBC 
% June, 2007

% initialize variables
[D N] = size(x);
K = mix.ncentres;
likIncr = options.threshold + eps;
logLambdaTilde = zeros(1,K);
E = zeros(N,K);
trSW = zeros(1,K);
xbarWxbar = zeros(1,K);
mWm = zeros(1,K);
trW0invW = zeros(1,K);

% priors
% prior over the mixing coefficients - CB 10.39
alpha0 = PriorPar.alpha;
% prior over gaussian means -  CB 10. 40
m0 = PriorPar.mu;
% prior over the gaussian variance - CB 10.40
beta0 = PriorPar.beta;
% wishart prior variables - CB 10.40
W0 = PriorPar.W;
W0inv = inv(W0);
v0 = PriorPar.v;

% Use 'responsibilities' from initialization to set sufficient statistics - 
% CB 10.51-10.53.
Nk = N*mix.priors';
xbar = mix.centres';
S = mix.covars;

% Use above sufficient statistics for M step update equations - CB 10.58, 
% 10.60-10.63. 
alpha = alpha0 + Nk;
beta = beta0 + Nk;
v = v0 + Nk;
m = ((beta0*m0)*ones(1,K) + (ones(D,1)*Nk').*xbar)./(ones(D,1)*beta');
W = zeros(D,D,K);
for k = 1:K
    mult1 = beta0.*Nk(k)/(beta0 + Nk(k));
    diff3 = xbar(:,k) - m0;
    W(:,:,k) = inv(W0inv + Nk(k)*S(:,:,k) + mult1*diff3*diff3');
end  
  

% Main loop of algorithm
for iter = 1:options.maxIter
  % Calculate r
  psiAlphaHat = psi(0,sum(alpha));
  logPiTilde = psi(0,alpha) - psiAlphaHat;
  const = D*log(2);
  for k = 1:K
    t1 = psi(0, 0.5*repmat(v(k)+1,D,1) - 0.5*[1:D]');
    logLambdaTilde(k) = sum(t1) + const  + log(det(W(:,:,k)));
    for n = 1:N
      % Calculate E
      diff = x(:,n) - m(:,k);
      E(n,k) = D/beta(k) + v(k)*diff'*W(:,:,k)*diff;
    end
  end
  logRho = repmat(logPiTilde' + 0.5*logLambdaTilde, N,1) - 0.5*E;
  logSumRho = logsumexp(logRho,2);
  logr = logRho - repmat(logSumRho, 1,K);
  r = exp(logr);
  
  % compute N(k)
  Nk = exp(logsumexp(logr,1))';
  % add a non-zero term for the components with zero responsibilities
  Nk = Nk + 1e-10;
  % compute xbar(k), S(k)
  for k=1:K
    xbar(:,k) = sum(repmat(r(:,k)',D,1).*x,2)/Nk(k);
    diff1 = x - repmat(xbar(:,k),1,N);
    diff2 = repmat(r(:,k)',D,1).*diff1;
    S(:,:,k) = (diff2*diff1')./Nk(k);
  end


  % compute Lower bound (refer to Bishop for these terms)
  % C(alpha0)
  logCalpha0 = gammaln(K*alpha0) - K*gammaln(alpha0);
  % B(lambda0)
  logB0 = (v0/2)*log(det(W0inv)) - (v0*D/2)*log(2) ...
          - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(v0+1 -[1:D])));
  % log(C(alpha))
  logCalpha = gammaln(sum(alpha)) - sum(gammaln(alpha));
  % Various other parameters for different terms
  H =0;
  for k = 1:K
    % sum(H(q(Lamba(k))))
    logBk = -(v(k)/2)*log(det(W(:,:,k))) - (v(k)*D/2)*log(2)...
            - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(v(k) + 1 - [1:D])));;
    H = H -logBk - 0.5*(v(k) -D-1)*logLambdaTilde(k) + 0.5*v(k)*D;
    % for Lt1 - third term
    trSW(k) = trace(v(k)*S(:,:,k)*W(:,:,k));
    diff = xbar(:,k) - m(:,k);
    xbarWxbar(k) = diff'*W(:,:,k)*diff;
    % for Lt4 - Fourth term
    diff = m(:,k) - m0;
    mWm(k) = diff'*W(:,:,k)*diff; 
    trW0invW(k) = trace(W0inv*W(:,:,k));
  end

  Lt1 = 0.5*sum(Nk.*(logLambdaTilde' - D./beta...
        - trSW' - v.*xbarWxbar' - D*log(2*pi)));
  Lt2 = sum(Nk.*logPiTilde);
  Lt3 = logCalpha0 + (alpha0 -1)*sum(logPiTilde);
  Lt41 = 0.5*sum(D*log(beta0/(2*pi)) + logLambdaTilde' - D*beta0./beta - beta0.*v.*mWm');
  Lt42 = K*logB0 + 0.5*(v0-D-1)*sum(logLambdaTilde) - 0.5*sum(v.*trW0invW');
  Lt4 = Lt41+Lt42;
  Lt5 = sum(sum(r.*logr));
  Lt6 = sum((alpha - 1).*logPiTilde) + logCalpha;
  Lt7 = 0.5*sum(logLambdaTilde' + D.*log(beta/(2*pi))) - 0.5*D*K - H;

  %Bishop's Lower Bound
  L(iter) = Lt1 + Lt2 + Lt3 + Lt4 - Lt5 - Lt6 - Lt7;

  % warning  if lower bound decreses
  if iter>2 && L(iter)<L(iter-1) 
    fprintf('Lower bound decreased by %f ', L(iter)-L(iter-1));
  end

  % Begin M step
  % compute new parameters
  alpha = alpha0 + Nk;
  beta = beta0 + Nk;
  v = v0 + Nk;
  m = (repmat(beta0.*m0,1,K) + repmat(Nk',D,1).*xbar)./repmat(beta',D,1);
  for k = 1:K
    mult1 = beta0.*Nk(k)/(beta0 + Nk(k));
    diff3 = xbar(:,k) - m0;
    W(:,:,k) = inv(W0inv + Nk(k)*S(:,:,k) + mult1*diff3*diff3');
  end

  %PLOT 
  if options.displayIter == 1
    fprintf('%d ',iter);
    fprintf('\n');
  end
  if options.displayFig == 1
    figure(3)
    clf
    plot(x(1,:),x(2,:),'o');
    hold on
    plot(m(1,:), m(2,:),'or','linewidth',2);
    weight = alpha/sum(alpha);
    for i = 1:K
      MyEllipse(inv(W(:,:,i))/(v(i)-D-1), m(:,i),'style','r','intensity',weight(i), 'facefill',.8);
      text(m(1,i), m(2,i), num2str(i),'BackgroundColor', [.7 .9 .7]);
    end
    pause(.01);
  end

  % check if the  likelihood increase is less than threshold
  if iter>1
    likIncr = abs((L(iter)-L(iter-1))/L(iter-1));
  end
  if likIncr < options.threshold
    break;
  end
end

out.alpha = alpha;
out.beta = beta;
out.m = m;
out.W = W;
out.v = v;
out.L = L;


