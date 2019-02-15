function X = genrndmvpexp(n,p,Mu,Sigma,Beta)
% rands = genrndmvpexp(n,p,Mu,Sigma,Beta)
%  Returns a sample matrix from multivariate Power Exponential distribution
%  with a given mean vector, covariance matrix, and kurtosis parameter. For
%  repeatability, the states for both randn and rand must be set.
%
%  Where
%  n --- scalar number of sample sets to be generated
%  p --- scalar number of variables
%  Mu --- (px1) mean vector
%  Sigma --- (pxp) covarince matrix
%  Beta --- scalar kurtosis parameter of PE distribution
%  rands --- (nxp) sample variable matrix
%
%  Example: pexp = genrndmvpexp(500,2,[0;0],[1,0;0,1],2);
%
%  Updated by J. Andrew Howe 20061230
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT
%
%  See Also GENRNDMIXMVPEXP.

if (nargin ~= 5) || (not(isscalar(p))) || (not(isscalar(n))) || (not(isscalar(Beta))) ...
        || (p < 1) || (n < 1) || (sum(size(Mu) == [p,1]) ~= 2) || ((sum(size(Sigma) == [p,p]) ~= 2))
    % not 5 args, nonscalar p/n/beta, nonpositive p or n, missized mu or sigma
    help genrndmvpexp, return
end

% Generate n-by-p random points uniformly distributed on the surface of a
% hypersphere of p dimension
a = normrnd(0,1,n,p);
a_t = a';
a1 = sum(a_t.^2)';
a2 = sqrt(a1*ones(1,p));
u = a./a2;

% Generate Gamma(1/2,p/2beta)
g = gamrnd(p/2/Beta,2,n,1);
R = g.^(1/2/Beta);

% Cholesky Factorization of SIGMA
A = chol(Sigma);

% Generate the PE matrix
X = (Mu*ones(1,n) + (ones(p,1)*R').*(A'*u'))';
