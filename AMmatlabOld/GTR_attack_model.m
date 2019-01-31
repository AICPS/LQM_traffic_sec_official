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

function [vuln_k1, vuln_k, vuln_xi, deo_01, deo_02, deo_final1d, deo_final2d]=GTR_attack_model(gT1, gT2, gT1_d, gT2_d, T_start, T_end, N)
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here,'ppt'));
    width=650;
    height=350;

    nMax=1000; %maximum iteration number
    threshold=10^(-5); %stopping threshold
    %% Signal Settings
    %xi1_arr=[0.3 0.5 0.85];
    %xi2_arr=[0.3 0.5 0.85];
%     xi1_arr = [0.5, 0.6, 0.85];
%     xi2_arr = [0.5, 0.6, 0.85]; % for diff green time attack
    xi1_arr = [0.5];
    xi2_arr = [0.5]; % for same green time attack
    %gT1=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.6 0.6 0.6  0.4 0.4 0.4 ];
    %gT2=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.4 0.4 0.4  0.6 0.6 0.6 ];
    %gT1_arr=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5 ];
    %gT2_arr=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5 ];
    %% Algorithm:
    % For each value of xi in each range
    % For each Targeted Initial State S0:
    %  For each Initial Values k10, k s.t. k > kj/2
    %     Check if k10, k satisfies S0 boundaries with initial xi values
    %      Desired Target Attack States S_attack:
    %      Check if k satisfies the boundaries of the attack state
    %      If yes, determine how to change pi1, pi2 based on the lambda table for the current state
    %      Save this behavior as attack_behavior
    %      For each pi1, pi2 modifications that satisfy attack_behavior:
    %      Simulate and determine if the system reaches state S in S_attack; 
    %      Store S and pi1_d, pi2_d into arrays
    %% Targeted Initial State = (2,6)
%     init_states = cat(1,[2,6],[3,7]); % for diff GTattack
%     init_states = cat(1, [3,7], [4,8]);
    %init_states = cat(1, [3,5], [1,1]);
    init_states = cat(1, [4,8], [8,8]); % for same GTattack
    target_states = cat(1, [3,8], [4,7]);
    %disp(init_states);
    %% Desired Target State is (3,8) or (4,7)  
    %% Possible intermediary states include: (2,7) or (3,6) for xi > 0.5 
    %% fundamental diagram settings
    vf=60;
    kj=150;
    kc=30;
    C=vf*kc;
    T = 90;

    %% signal settings
    L1 = 0.25;
%     % for state-changing attacks
%     gT1 = 0.5;
%     gT2 = 0.5;
%     % for non state-changing attacks
%     gT1 = 0.4;
%     gT2 = 0.4;
    
    %% Algorithm Settings
    deo_period = 1;
    k10_arr = 0:kj;
    k_arr = kj/2:kj;
    vuln_k1 = [];
    vuln_k = [];
    vuln_xi = [];
    deo_01 = [];
    deo_02 = [];
    deo_final1 = [];
    deo_final2 = [];
    deo_final1d = [];
    deo_final2d = [];
    allk1_array = [];
    allk1d_array = [];
    allq_array = [];
    allqd_array = [];
    
    % Writing Settings
%     fname = strcat('GTR_dr_4-7_', num2str(gT1), '_', num2str(gT2), '_', num2str(gT1_d), '_', num2str(gT2_d), '_', num2str(T_start), '_', num2str(T_end), '_', num2str(N), '_', num2str(T), '_', num2str(xi1_arr[1]));
%     write_nfname = strcat('./data/double/',fname, '.xlsx');
%     if exist(write_nfname, 'file') 
%         delete(write_nfname); 
%     end
    
    %% ALGORITHM START
    % For each value of xi in each range
    for a = 1:length(xi1_arr)
        xi1 = xi1_arr(a);
        xi2 = xi2_arr(a);
        % Writing Settings
        fname = strcat('new_GTR_dr_4-8_', num2str(gT1), '_', num2str(gT2), '_', num2str(gT1_d), '_', num2str(gT2_d), '_', num2str(T_start), '_', num2str(T_end), '_', num2str(N), '_', num2str(T), '_', num2str(xi1));
        write_nfname = strcat('./data/double/',fname, '.xlsx');
        if exist(write_nfname, 'file') 
            delete(write_nfname); 
        end
        % For each Targeted Initial State S0:
        %  For each Initial Values k10, k s.t. k > kj/2

        % Initial k1 values

        for i = 1:length (k10_arr)
            k10 = k10_arr(i);
            % Initial k values
            for j = 1:length(k_arr)
                k = k_arr(j);
                fprintf('k1 is: %d\n', k10);
                fprintf('k is: %d\n', k);
                k_2 = 2*k_arr(j) - k10; 
                % Check if k10, k satisfies S0 boundaries with initial xi values
%                 disp(length(init_states));
                for b = 1:length(init_states)
%                     disp(length(init_states))
                    S0_1 = init_states(b,1);
                    S0_2 = init_states(b,2);
                    %fprintf('Current xi is: (%f)\n', xi1);
                    %fprintf('Desired Initial DEO is: (%d, %d)\n', S0_1, S0_2);
                    [r1, r2] = getDEO(xi1, xi2, gT1, gT2, k10, k,T);
                    fprintf('Current DEO is: (%d, %d)\n', r1, r2');
                    if S0_1 == r1 && S0_2 == r2 %continue simulation; else skip
                        % Desired Target Attack States S_attack:
                        % Check if k satisfies the boundaries of the attack state
                        deos = []; % keep track of the deo at every other period
                        deos_d = []; % deos as a result of attack
                        q_array = []; % array of average flows after each cycle without attack
                        q_array_d = []; % array of average flow after each cycle with attack
                        k1_array = [k10];
                        k1_d_array = [k10];
                        k1_prev = k10;
                        k1_d_prev = k10;
                        k1_next = 0;
                        k1_d_next = 0;
                        k = k_arr(j);
                        %fprintf('k1 is: %d\n', k10);
                        %fprintf('k is: %d\n', k);
                        if(S0_1 == 3 && S0_2 == 7 && k ~= k10)
                            % do nothing
                        else
                            % number of cycles
                            for n = 1:N 
                                k_2 = 2*k - k1_prev; 
                                k_2d = 2*k - k1_d_prev;
                                if(k_2 <= kj && k_2d <= kj) % make sure we don't test unreasonable states
                                    if (n > T_start && n <= T_end) % attack window
                                        [q, k1_next, k2] = flow_calculation(k,k1_prev,xi1,xi2,T,gT1,gT2, vf, kj, kc, L1);
                                        [q_d, k1_d_next, k2] = flow_calculation(k,k1_d_prev,xi1,xi2,T,gT1_d,gT2_d, vf, kj, kc, L1);
                                        q_array = [q_array q];
                                        q_array_d = [q_array_d q_d];
                                        k1_array = [k1_array k1_next];
                                        k1_d_array = [k1_d_array k1_d_next];
                                        k1_prev = k1_next;
                                        k1_d_prev = k1_d_next;
                                        if(mod(n, deo_period) == 0)
                                            [r1,r2] = getDEO(xi1, xi2, gT1, gT2, k1_prev, k,T);
                                            [r1d,r2d] = getDEO(xi1, xi2, gT1_d, gT2_d, k1_d_prev, k,T);
                                            deos = [deos  [r1,r2]];
                                            deos_d = [deos_d [r1d,r2d]];
                                        end
                                    else % no attack
                                        [q, k1_next, k2] = flow_calculation(k,k1_prev,xi1,xi2,T,gT1,gT2, vf, kj, kc, L1);
                                        [q_d, k1_d_next, k2] = flow_calculation(k,k1_d_prev,xi1,xi2,T,gT1,gT2, vf, kj, kc, L1);
                                        q_array = [q_array q];
                                        q_array_d = [q_array_d q_d];
                                        k1_array = [k1_array k1_next];
                                        k1_d_array = [k1_d_array k1_d_next];
                                        k1_prev = k1_next;
                                        k1_d_prev = k1_d_next;
                                        if(mod(n, deo_period) == 0)
                                            [r1,r2] = getDEO(xi1, xi2, gT1, gT2, k1_prev, k,T);
                                            [r1d,r2d] = getDEO(xi1, xi2, gT1, gT2, k1_d_prev, k,T);
                                            deos = [deos  [r1,r2]];
                                            deos_d = [deos_d [r1d,r2d]];
                                        end
                                    end
                                end
                            end
                            if(k_2 <= kj && k_2d <= kj)
                                k_2 = 2*k - k10; 
                                vuln_k1 = [vuln_k1 k10];
                                vuln_k = [vuln_k k];
                                vuln_xi = [vuln_xi xi1];
                                deo_01 = [deo_01 S0_1];
                                deo_02 = [deo_02 S0_2];
                                deo_final1 = [deo_final1 r1];
                                deo_final2 = [deo_final2 r2];
                                deo_final1d = [deo_final1d r1d];
                                deo_final2d = [deo_final2d r2d];
                                allq_array = [allq_array; q_array];
                                allqd_array = [allqd_array; q_array_d];
                                allk1_array = [allk1_array; k1_array];
                                allk1d_array = [allk1d_array; k1_d_array];
                            end
                        end
                    end
                end
            end
        end
    end
    disp('Writing now...');
    xlswrite(write_nfname, [vuln_k1' vuln_k' vuln_xi'], 1);
    xlswrite(write_nfname, [deo_01'  deo_02'], 2);
    xlswrite(write_nfname, [deo_final1'  deo_final2'], 3);
    xlswrite(write_nfname, [deo_final1d' deo_final2d'], 4);
    xlswrite(write_nfname, allq_array,5);
    xlswrite(write_nfname, allqd_array, 6);
    xlswrite(write_nfname, allk1_array, 7);
    xlswrite(write_nfname, allk1d_array, 8);
%     disp(allqd_array);
end
