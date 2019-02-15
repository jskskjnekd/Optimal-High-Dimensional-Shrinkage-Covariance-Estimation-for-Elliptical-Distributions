function [beta_O_Ell, alpha_O_Ell] = get_Ell_beta_alpha(X)
%
%   X:  (n * p)


[n, p] = size(X);

S = get_S(X);
S_sgn = get_S_sgn(X);

gamma_hat = p* trace(S_sgn^2) - p/n;
k_j_hat = kurtosis(X) - 3;
kappa_hat = max(-2/(p+2), ((sum(k_j_hat))/(3*p)));
T = gamma_hat -1;


beta_O_Ell = max(0, T/(T + (1/n)*( kappa_hat*(2*gamma_hat + p)+gamma_hat+p ) ));
alpha_O_Ell = (1-beta_O_Ell)*trace(S)/p;



end

