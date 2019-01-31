% iterate through arrays for xi, pi, T, pi_delta (attack model) 

% the following initial settings and attack settings are for studying "same green
% time ratio attacks"
xi1_arr = [0.5]; % the retaining ratios
xi2_arr = [0.5];

pi1_arr=[0.3]; % the initial green time ratios
pi2_arr=[0.3];

T_arr = [30, 60, 90, 120]; % array of cycle lengths
start_arr = [1, 2, 4, 10]; % array of attack start times
end_arr = [6, 6, 10, 15]; % array of attack end times (not inclusive)

k10_arr = 0:kj; % array of initial state ring 1 density values
k_arr = kj/2:kj; % array of initial state average network density values
T = 90; % cycle length
N = 25; % simulation length (in cycles)

pi1d_arr=[0.35, 0.4, 0.5]; % array of green time ratios during attack window
pi2d_arr =[0.35, 0.4, 0.5]; 

init_states = cat(1, [3,8], [8,8]); % 3,8 or 4,7 % these are DEOs of 
                                    %initial states that we want to test
                                    % [8, 8] is just fodder
target_states = cat(1, [3,8], [4,7]); 

% simulate for both single intersection and grid network
for i = 1:length(xi1_arr)
    for j = 1:length(pi1_arr)
        for l = 1:length(start_arr)
            if(start_arr(l) <= end_arr(l))
                for m = 1:length(pi1d_arr)
                    xi1 = xi1_arr(i);
                    xi2 = xi2_arr(i);
                    pi1 = pi1_arr(j);
                    pi2 = pi2_arr(j);
                    Tstart = start_arr(l);
                    Tend = end_arr(l);
                    pi1d = pi1d_arr(m);
                    pi2d = pi2d_arr(m);
                    [vuln_k1, vuln_k, vuln_xi, deo_01, deo_02, ...
                        deo_final1, deo_final2] = ...
                    GTR_attack_model(xi1, xi2, pi1, pi2, pi1d, pi2d, ...
                        Tstart, Tend, T, N, k10_arr, k_arr,...
                        init_states, target_states);
                    [gvuln_k1, gvuln_k, gvuln_xi, gdeo_01, gdeo_02, ...
                        gdeo_final1, gdeo_final2] = ...
                    GTR_grid_attack_model(xi1, xi2, pi1, pi2, pi1d, pi2d,...
                        Tstart, Tend, T, N, k10_arr, k_arr,...
                        init_states, target_states);
                end
            end
        end
    end
end
