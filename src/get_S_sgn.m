function [S_sgn] = get_S_sgn(X)
%
%   X:  (n * p)

[n, p] = size(X);

S_sgn = zeros(p, p);


for i = 1:n
    S_sgn = S_sgn + (X(i, :)'*X(i, :))/((norm(X(i, :)))^2);
end

S_sgn = (1/n)*S_sgn;


end

