


clc
clear
close all




% -------------- Basic Parameters
p = 100;

c_inv = [0.2:0.03:1.2];
NMSE_list = zeros(length(c_inv),1);
varrho = 0.4;

for i=1:length(c_inv)
    n = floor(p*c_inv(i));
    NMSE_list(i,1) = simulate_mv_t_dist(p, n, varrho);
end


figure
plot(c_inv, NMSE_list, 'r-o','linewidth',2)
xlabel('c^{-1}')
ylabel('NMSE')
title('Given \rho=0.4')


