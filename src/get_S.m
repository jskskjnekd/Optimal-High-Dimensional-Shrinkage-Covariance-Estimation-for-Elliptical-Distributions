function [S] = get_S(X)
%
%   X:  (n * p)

[n, p] = size(X);

S = zeros(p, p);


for i = 1:n
    S = S + X(i, :)'*X(i, :);
end

S = (1/n)*S;




end

