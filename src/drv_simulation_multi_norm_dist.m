


clc
clear
close all




% -------------- Basic Parameters
p = 100;

num_simulations = 10;
c_inv = [0.2:0.02:1.2];



varrho_1 = 0.2;
C_1_list_1 = single_simulation_multi_norm(p, c_inv, varrho_1, num_simulations);
varrho_2 = 0.4;
C_1_list_2 = single_simulation_multi_norm(p, c_inv, varrho_2, num_simulations);
varrho_3 = 0.7;
C_1_list_3 = single_simulation_multi_norm(p, c_inv, varrho_3, num_simulations);

figure
loglog(c_inv, mean(C_1_list_1'), 'r-o')
hold on
loglog(c_inv, mean(C_1_list_2'), 'g-o')
hold on
loglog(c_inv, mean(C_1_list_3'), 'b-o')
xlabel('n/p')
ylabel('C_1')
legend('0.2', '0.4', '0.7')

