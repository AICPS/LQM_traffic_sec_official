%% Description:
% The flow of the attack model is as follows:
% The initial condition must satisfy: k > kj /2 at least. 
% The only possible states that satisfy this constraint are:(2,6), (2,7),
% (4,7), (3,6), (3,7), (3,8) when xi > 0.5
% However, since we are concerned more about non gridlock states and states
% that are stationary, then we should only consider: (2,6), (3,7)
% However, whether or not an initial (k1,k) can be in one of these states
% depends on xi (retaining ratio). Therefore, we compute the region
% boundaries according to xi and determine which possible combination of
% regions can (k1,k) can be in. 
% Then, according to the region boundaries, the behavior of the current
% DEO, we can determine if a target region is reachable. 
% This can also be determined with the values of pi1 and pi2, lambda, and
% the boundaries of the target region. Any possible intermediary states
% must also have the same lambda as the initial state. 

function [vuln_k1, vuln_k, vuln_xi, deo_01, deo_02, deo_final1, deo_final2]=...
    GTR_grid_attack_model(xi1, xi2, gT1, gT2, gT1_d, gT2_d, ... 
    T_start, T_end, T, N, k10_arr, k_arr, init_states, target_states)
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here,'ppt'));
    width=650;
    height=350;
    nMax=1000; %maximum iteration number
    threshold=10^(-5); %stopping threshold
    %init_states = cat(1, [3,5], [1,1]);
%     init_states = cat(1, [4,8], [8,8]);
%     target_states = cat(1, [3,8], [4,7]);
    %disp(init_states);
    %% Desired Target State is (3,8) or (4,7)  
    %% Possible intermediary states include: (2,7) or (3,6) for xi > 0.5 
    %% fundamental diagram settings
    vf=60;
    kj=150;
    kc=30;
    C=vf*kc;
%     T = 90;
    %% network and signal settings
    L1 = 0.25;
    W = 4;
    H = 4;
    
    %% Algorithm Settings
    k10 =zeros(W*H,1); 
    k2 = zeros(W*H,1);
    k2d = zeros(W*H,1);
    k = zeros(W*H,1);
%     k10_arr = 0:kj;
%     k_arr = kj/2:kj;
    vuln_k1 = []; % (k_1, k) states that end in state with attack DEO after attack
    vuln_k = []; %
    vuln_xi = [];
    deo_01 = []; % first region R1 of DEO (R1, R2)
    deo_02 = []; % second region R2 of DEO (R1, R2)
    deo_final1 = []; % DEO after attack
    deo_final2 = []; 
    allq = zeros(length(k10_arr)*length(k_arr), N);
    allqd = zeros(length(k10_arr)*length(k_arr), N);
    allk1 = zeros(length(k10_arr)*length(k_arr), N);
    allk1d = zeros(length(k10_arr)*length(k_arr), N);
    deos = []; % keep track of the deo at every other period
    deos_d = []; % deos as a result of attack
    q_array = []; % array of average flows after each cycle without attack
    q_array_d = []; % array of average flow after each cycle with attack
    k1_array = []; % keeps track of k1 densities after each cycle
    k1_d_array = []; % keeps track of attacked k1 densities
    k1_prev = []; % keeps track of previous time step k1 densities
    k1_d_prev = []; % keeps track of previous time step k1 densities under attack
    k1_next = 0; % the next time step k1 densities
    k1_d_next = 0; % the next time step k1 densities under attack
    ctr = 1;
    
%     fname = strcat('gtrattack_gn_',num2str(init_state(1,1)),',',...
%         num2str(init_state(1,2)),'_', num2str(gT1), '_', num2str(gT2), '_', ...
%         num2str(gT1_d), '_', num2str(gT2_d), '_', num2str(T_start), '_', ...
%         num2str(T_end), '_', num2str(N), '_', num2str(T), '_', num2str(xi1));
%     write_nfname = strcat('./data/grid/',fname, '.xlsx');
%     if exist(write_nfname, 'file') 
%         delete(write_nfname); 
%     end
    % Initial k1 values
    for i = 1:length (k10_arr)
        k10 = k10_arr(i);
        % Initial k values
        for j = 1:length(k_arr)
            k = k_arr(j);
%                 fprintf('k1 is: %d\n', k10(1));
%                 fprintf('k is: %d\n', k(1));
            k2 = 2*k_arr(j) - k10; 
            k2d = k2;
            if(k2 <= kj && k2d <= kj)
                % Check if k10, k satisfies S0 boundaries with initial xi values
                for b = 1:length(init_states)
                    S0_1 = init_states(b,1);
                    S0_2 = init_states(b,2);
                    [r1, r2] = getDEO(xi1, xi2, gT1, gT2, k10, k,T);
                    fprintf('Xi is: %.2f, GTR is: %.2f, GTRD is: %.2f, (k1,k) is (%d, %d), DEO is (%d, %d)\n',...
                        xi1, gT1, gT1_d, k10, k, r1, r2);
                    if S0_1 == r1 && S0_2 == r2 %continue simulation; else skip
                        % N = specified number of cycles
                        if(S0_1 == 3 && S0_2 == 7 && k(1) ~= k10(1))
                            % do nothing
                        else
                            fname = strcat('gtrattack_gn_',num2str(r1),',',...
                                num2str(r2),'_', num2str(gT1), '_', num2str(gT2), '_', ...
                                num2str(gT1_d), '_', num2str(gT2_d), '_', num2str(T_start), '_', ...
                                num2str(T_end), '_', num2str(N), '_', num2str(T), '_', num2str(xi1));
                            write_nfname = strcat('./data/grid/',fname, '.xlsx');
                            if exist(write_nfname, 'file') 
                                delete(write_nfname); 
                            end
                            % uncomment the following for simulation
                            % results without attack
                            %if(max(k2) <= kj && max(k2d) <= kj) % make sure we don't test unreasonable states
%                           [qa, q, k1_next, k2_next] = LQM_grid(k(1),k10(1),xi1,xi2,T,N*T,gT1,gT2,vf,kj,kc,L1,W,H);
                            %q_array = [q_array q];
%                           xlswrite(write_nfname, qa, 4, strcat('A',num2str(ctr)));
%                           deo_final1 = [deo_final1 r1d];
%                           deo_final2 = [deo_final2 r2d];
                            % Simulate attack model on LQM grid network
                            [qda, q_d, k1_d_next, k2d_next] = ...
                                LQM_grid_attack(k(1),k10(1),xi1,xi2,T,N*T,T_start,T_end,...
                                 gT1, gT2, gT1_d,gT2_d,vf,kj,kc,L1,W,H);
                            q_array_d = [q_array_d q_d];
                            xlswrite(write_nfname, qda, 5, strcat('A',num2str(ctr)));
                            vuln_k1 = [vuln_k1 k10(1)];
                            vuln_k = [vuln_k k(1)];
                            vuln_xi = [vuln_xi xi1];
                            deo_01 = [deo_01 S0_1];
                            deo_02 = [deo_02 S0_2];
                            ctr = ctr+1;
                            %end
                        end
                    end
                end
            end
        end
    end
    
    
    disp('Writing now...');
%     disp(q_array);
    xlswrite(write_nfname, [vuln_k1' vuln_k' vuln_xi'], 1);
    xlswrite(write_nfname, [deo_01'  deo_02' ], 2);
    xlswrite(write_nfname, [q_array' q_array_d'], 4);

%     xlswrite(write_nfname, [deo_final1'  deo_final2'], 3);
%     xlswrite(write_nfname, allq, 4);
%     xlswrite(write_nfname, allqd, 5);

end
