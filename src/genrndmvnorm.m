function X = genrndmvnorm(n,p,Mu,Sigma)
% rands = genrndmvnorm(n,p,Mu,Sigma)
%  Returns a data set from multivariate normal distribution with a given 
%  mean vector and a given covariance matrix.
%
%  Where
%  n --- scalar number of observations to be generated
%  p --- scalar number of variables
%  Mu --- (1xp) mean vector
%  Sigma --- (pxp) covariance matrix
%  rands --- (nxp) generated sample data matrix
%
%  Example: norm = genrndmvnorm(1000,2,[0,0],[2,0;0,1]);
%
%  Written by   Yuehui Fan
%  Updated by   J. Andrew Howe 20060523
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT
%
%  See Also PLOTRNDMVNORM, SRFCNTMVNORM, SRFCNTMIXMVNORM, GENRNDMIXMVNORM.

if (nargin ~= 4) || (not(isscalar(p))) || (not(isscalar(n))) || (p < 1) || ...
        (n < 1) || (sum(size(Mu) == [1,p]) ~= 2) || ((sum(size(Sigma) == [p,p]) ~= 2))
    % not 4 args, nonscalar p or n, nonpositive p or n, missized mu or Sigma
    help genrndmvnorm, return
end

T = chol(Sigma); Xtem = randn(n,p)*T; X = Xtem + ones(n,1)*Mu;


