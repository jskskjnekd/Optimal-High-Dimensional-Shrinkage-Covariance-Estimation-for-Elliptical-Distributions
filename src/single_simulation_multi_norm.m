function [C_1_list] = single_simulation_multi_norm(p, c_inv, varrho, num_simulations)


C_1_list = zeros(length(c_inv),num_simulations);

for j=1:num_simulations
    for i=1:length(c_inv)
        n = floor(p*c_inv(i));
        C_1_list(i,j) = simulate_mv_n_dist_C1(p, n, varrho);
    end
end



end

