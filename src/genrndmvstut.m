function X = genrndmvstut(n,p,Mu,Sigma,v)
% rands = genrndmvstut(n,p,Mu,Sigma,nu)
%  Returns a sample matrix from multivariate student's t distribution with
%  a given mean vector, covariance matrix, and degrees of freedom.
%
%  Where
%  n --- scalar number of sample sets to be generated
%  p --- scalar number of variables
%  Mu --- (1xp) mean vector
%  Sigma --- (pxp) covariance matrix: Sigma = cov(X) = (nu/(nu-2))*H^(-1)
%  nu --- scalar degrees of freedom for (2 < nu < 170-p); if nu > 170-p,
%     use genrndmvnorm.
%  rands --- (nxp) sample variable matrix
%
%  Example: stut = genrndmvstut(1000,2,[0,0],[2,0;0,1],3);
%
%  Created By Yuehui Fan     11-16-90
%  Updated By J. Andrew Howe 20060523
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT
%
%  See Also PLOTRNDMVSTUT, SRFCNTMVSTUT, SRFCNTMIXMVSTUT, GENRNDMIXMVSTUT.

if (nargin ~= 5) || (not(isscalar(p))) || (not(isscalar(n))) || (not(isscalar(v))) ...
        || (p < 1) || (n < 1) || (sum(size(Mu) == [1,p]) ~= 2) || ((sum(size(Sigma) == [p,p]) ~= 2))
    % not 5 args, nonscalar p/n/v, nonpositive p or n, missized mu or sigma
    help genrndmvstut, return
end


H = (v/(v-2))*inv(Sigma);
Sd = diag(Sigma)'; ex = (-(p+v)/2);
c1 = gamma((p+v)/2)*(det(H))^(1/2)/((v*pi)^(p/2)*gamma(v/2));
X = zeros(n,p); k = 1;

while k < n+1
      rnd = rand(1,p); tmp = 8*Sd.*rnd + Mu - 4*Sd;
      a = (tmp-Mu)*H*(tmp-Mu)';
      f1 = c1*(1+a/v)^ex; f2 = c1*rand;
      if f2 <= f1
         X(k,:) = tmp; k = k+1;
      end
end
