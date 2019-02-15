


clc
clear
close all




% -------------- Basic Parameters
p = 100;
n = p * 0.2;
Mu = ones(p,1);
SigmaMatrix = zeros(p, p);
varrho = 0.4;

for i = 1:p
    for j = 1:p
        SigmaMatrix(i, j) = varrho^(abs(i-j));
    end
end

X = genrndmvnorm(n, p, Mu', SigmaMatrix);


[beta_O_Ell, alpha_O_Ell] = get_Ell_beta_alpha(X)


