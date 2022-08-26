function [Px, model, BIC] = GaussianMM(x, K, epoch)
% ============================================================
%A function to apply the Gaussian mixture model using E-M algorithm and to calculate the BIC
%value
%
%  Notation used :-
%  Input:
%  x: N-by-D data matrix.
%  K: the number of components (Interger)
%  epoch: the maiximum number of iterations
%
%  Output:
%  PX: N-by-K matrix indicating the probability of each
%       component generating each point.
%  MODEL: a structure containing the parameters for a GMM:
%       MODEL.Mu: a K-by-D matrix.
%       MODEL.Sigma: a D-by-D-by-K matrix.
%       MODEL.Tau: a 1-by-K vector.
%  BIC: Bayesian information criterion for the specified K



% ============================================================
% @SourceCode Author: Pluskid (http://blog.pluskid.org)
% @Edited by : T.Yeo (2022 August)
    

    %% Generate Initial Centroids
    threshold = 1e-10;
    
    tol = 1e-6;
    llh = -inf(1,epoch);

    [N, D] = size(x);
 
    rndp = randperm(N);
    centroids = x(rndp(1:K),:);
  
    %% initial values
    [pMu,pTau,pSigma] = theta();
    Lprev = -inf; 
    
    %% EM Algorithm
    for iter = 2:epoch
        %% Estimation Step
        Px = real(Gaussian_prob());
        Px(isnan(Px))=0;
        Px(isnan(Px))=0;
       
        % calculate responsibility
        pGamma = Px .* repmat(pTau, N, 1); %numerator  = Tau(k) * N(xi | pMu(k), pSigma(k))
        pGamma = pGamma ./ repmat(sum(pGamma, 2), 1, K); %denominator = Tau(j) * N(xi | pMu(j), pSigma(j))sum over j
        
        %% Maximization Step - through Maximize likelihood Estimation
        
        Nk = sum(pGamma, 1);    %number of samples in each cluster
        Nk(isnan(Nk))=0;
        Nk(isinf(Nk))=0;
        
        % update pMu
        pMu = diag(1./Nk) * pGamma' * x;
        pTau = Nk/N;
        
        
        % update pSigma
        for kk = 1:K 
            Xshift = x-repmat(pMu(kk, :), N, 1);
            pSigma(:, :, kk) = (Xshift' * ...
                (diag(pGamma(:, kk)) * Xshift)) / Nk(kk);
        end
 
        % check for convergence
        L = sum(log(Px*pTau'));
        if L-Lprev < threshold || abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter))
            break;
        end
        Lprev = L;
        
    end
 
    %% Bayesian information criterion
    BIC=(K*log(N))-(2*L);

    %% Output
        model = [];
        model.Mu = pMu';
        model.Sigma = pSigma;
        model.Tau = pTau;

    %% Function Definition
    
    function [pMu,pTau,pSigma] = theta()
        pMu = centroids; 
        pTau = zeros(1, K); 
        pSigma = zeros(D, D, K); 
 
        % hard assign x to each centroids
        distmat = repmat(sum(x.*x, 2), 1, K) + ... 
            repmat(sum(pMu.*pMu, 2)', N, 1) - ...
            2*x*pMu';
        [~, labels] = min(distmat, [], 2);%Return the minimum from each row
 
        for k=1:K
            Xk = x(labels == k, :);
            pTau(k) = size(Xk, 1)/N;
            pSigma(:, :, k) = cov(Xk);
        end
    end
 
    function Px = Gaussian_prob() 
        %Gaussian posterior probability 
        %N(x|pMu,pSigma) = 1/((2pi)^(D/2))*(1/(abs(sigma))^0.5)*exp(-1/2*(x-pMu)'pSigma^(-1)*(x-pMu))
        Px = zeros(N, K);
        for k = 1:K
            Xshift = x-repmat(pMu(k, :), N, 1); %X-pMiu
            inv_pSigma = inv(pSigma(:, :, k));
            tmp = sum((Xshift*inv_pSigma) .* Xshift, 2);
            coef = (2*pi)^(-D/2) * sqrt(det(inv_pSigma));
            Px(:, k) = coef * exp(-0.5*tmp);
        end
    end
end