function [NMSE] = simulate_mv_n_dist(p, n, varrho)
%
%  p:  number of variables
%  n:  number of observations



Mu = zeros(p,1);
SigmaMatrix = zeros(p, p);

for i = 1:p
    for j = 1:p
        SigmaMatrix(i, j) = varrho^(abs(i-j));
    end
end

X = genrndmvnorm(n, p, Mu', SigmaMatrix);
S = get_S(X);
[beta_O_Ell, alpha_O_Ell] = get_Ell_beta_alpha(X);
S_alpha_beta_Ell = beta_O_Ell * S + alpha_O_Ell * eye(p);

 NMSE = (norm(SigmaMatrix-S_alpha_beta_Ell, 'fro'))^2/(norm(SigmaMatrix, 'fro')^2);





end

