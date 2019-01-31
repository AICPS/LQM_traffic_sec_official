function [cumQa, cumQ, k1_ret, k2_ret]=LQM_grid(k,k10,xi1,xi2,cycle,Tsim,gT1,gT2,vf,kj,kc,L1,W,H)

% need to make two matrixes
% one for the E/W links -> aka ring 1
% one for the N/S links -> aka ring 2

k20=2*k-k10;

%simulation settings
ddt=0.05;
dt=ddt/3600;%2 seconds simulation time step;in hours

%signal settings
NT=cycle*(1/ddt);%1 minute for each cycle -> total time * (1/0.05)
nLostTime=(1-gT1-gT2)*cycle/2/ddt; %total steps for the lost time in each phase
nGreen1=(gT1*cycle)/ddt;

% disp(NT);
% disp(nLostTime);
% disp(nGreen1);

N=Tsim*(1/ddt);
M=1;
%disp(N);

% two arrays corresponding to the average densities for each direction
% hlinks = zeros(N,W,1);
% vlinks = zeros(N,H,1);

%initialization
f1=zeros(N,W*H);
f2=zeros(N,W*H);
g1=zeros(N,W*H);
g2=zeros(N,W*H);
d1=zeros(N,W*H);
d2=zeros(N,W*H);
s1=zeros(N,W*H);
s2=zeros(N,W*H);
k1=zeros(N+1,W*H);
k2=zeros(N+1,W*H);
k1_tem=zeros(N,W*H);
k2_tem=zeros(N,W*H);
T=zeros(N+1,1);
cumQa = [];
%W= 4; %width - number of links -> horizontal
%H = 4; %height  - number of links -> vertical
% for each intersection, there are two inputs and two outputs of
% edges/links
% this will help us find which links these are
H_input = zeros(H, W); % horizontal links inputs
H_output = zeros(H, W); % horizontal links outputs
V_input = zeros(H, W); % vertical inputs
V_output = zeros(H, W); % vertical outputs

% fill in the link values for each direction
% horizontal direction
idx = 1;
i = 1;
j = 1;
for ctr = 1:W*H
    H_input (i,j)= idx; % idx corresponds to link id number
    idx = mod(idx,W*H)+1;
    H_output(i,j) = idx;
    if(mod(i,2) == 1) % odd rows -> east 
        if j == W % take it out for other version
            H_output (i,j) = H_output(i,1)-1;
        end
        j = j+1;
        if (j > W) % wrap around correctly
            j = W;
            i = mod(i,H)+1;
        end
    else % even rows -> west 
        if j == 1 % take it out for other version
            H_output (i,j) = H_output(i,W)-1;
        end
        j = j-1;
        if (j < 1) % wrap around correctly
            j = 1; 
            i = mod(i,H)+1;
        end
    end
end

%vertical direction
idx = 1;
i = 1;
j = 1;
for ctr = 1: W*H
    V_input (i,j) = idx; % id of input link
    idx = mod(idx,W*H) +1;
    V_output (i,j) = idx;
    if(mod(j,2) == 1) % odd columns -> south 
        if i == H % take it out for other version
            V_output (i,j) = V_output(1,j)-1;
        end
        i = i+1;
        if(i > H) % wrap around correctly
            i = H;
            j = mod(j,W)+1;
        end
    else % even columns -> north
        if i == 1 % take it out for other version
            V_output(i,j) = V_output(H,j)-1;
        end
        i = i-1;
        if(i<1) % wrap around correctly
            i = 1;
            j = mod(j,W)+1;
        end
    end
end

k1(1,:)=k10; % set the initial link densities
k2(1,:)=2*k-k10;
% k1pre = k1(1,:);
% k2pre = k2(1,:);
% need mapping of intersections with links
% time steps

%disp('START_LQMGRIDNETWORK');
for t = 1:N
    i = 1;
    j = 1;
    % go through intersections in horizontal direction and update the flows
    % in each
    %fprintf('Cycle is: %d\n', round(t/NT));
    %k1pre = k1(t,:);
    %k2pre = k2(t,:);
    %fprintf('Current time step is %d and k1pre, k2pre is: (%d, %d)\n', t, k1pre, k2pre);
    %determine the phase
    if(mod(t,NT)<=nGreen1)
        signal=1;
    elseif(mod(t,NT)>nGreen1+nLostTime && mod(t,NT)<=NT-nLostTime)
        signal=2;
    else
        signal=0;
    end

    % pre-calculate and set the initial demand and supplies for the links
    for ctr = 1:H*W
        %fprintf('Current coordinates are (%d, %d)\n', i, j);
        % get correct ids of horizontal and vertical links
        h_up = H_input(i,j);
        h_down = H_output(i,j);
        v_up = V_input(i,j);
        v_down = V_output(i,j); 

%         fprintf('Current inputs are (%d, %d)\n', h_up, v_up);
%         fprintf('Current outputs are (%d, %d)\n', h_down, v_down);
        % update demand and supply; update link densities corresponding
        % to current intersection coordinate (i,j)
        C=kc*vf;
        % demand from upstream horizontal of current intersection
        d1(t,h_up)=min(vf*k1(t,h_up),C);
        % demand from upstream vertical
        d2(t,v_up)=min(vf*k2(t,v_up),C);
        % supply from downstream horizontal
        s1(t,h_down)=min(C,C*(kj-k1(t,h_down))/(kj-kc));
        % supply from downstream vertical
        s2(t,v_down)=min(C,C*(kj-k2(t,v_down))/(kj-kc));
        k1_tem(t,h_up)=k1(t,h_up);
        k2_tem(t,v_up)=k2(t,v_up);
        T(t)=dt*t;

        if(signal==1)%ring 1 goes first
            %outflow horizontal upstream = Demand horizontal upstream,
            %demand horizontal downstream, supply vertical downstream
            g1(t,h_up)=min(d1(t,h_up),min(s1(t,h_down)/xi1,s2(t,v_down)/(1-xi1)));
            % outflow vertical upstream = demand vertical upstream, supply
            % vertical downstream , supply horizontal downstream
            g2(t,v_up)=0;
            % inflow horizontal downstream = outflow horizontal upstream,
            f1(t,h_down)=g1(t,h_up)*xi1+g2(t,v_up)*(1-xi2);
            % inflow vertical downstream = horizontal upstream ,
            % vertical upstream
            f2(t,v_down)=g1(t,h_up)*(1-xi1)+g2(t,v_up)*xi2;
        elseif(signal==2)%ring 2 has green time
            g1(t,h_up)=0;
            g2(t,v_up)=min(d2(t,v_up),min(s2(t, v_down)/xi2,s1(t,h_down)/(1-xi2)));
            f1(t,h_down)=g1(t,h_up)*xi1+g2(t,v_up)*(1-xi2);
            f2(t,v_down)=g1(t,h_up)*(1-xi1)+g2(t,v_up)*xi2;
        else %during the lost times SUPER IMPORTANT
            g1(t,h_up)=0;
            g2(t,v_up)=0;
            f1(t,h_down)=g1(t,h_up)*xi1+g2(t,v_up)*(1-xi2);
            f2(t,v_down)=g1(t,h_up)*(1-xi1)+g2(t,v_up)*xi2;
        end
        % change in average link density of i,j -> update link density
        % of i+1, j or i, j+1
%        fprintf('h down is: %d\n', h_down);
%        fprintf('v_down is: %d\n', v_down);

%         k1(t+1,h_down)= min(kj, max(k1(t,h_down)+dt*(f1(t,h_down)-g1(t,h_down))/L1, 0));
%         k2(t+1,v_down) = min(kj, max(k1(1,h_down) + k2(1,v_down)-k1(t+1, h_down),0));
        if(mod(i,2) == 1) % odd rows -> east 
            j = j+1;
            if (j > W) % wrap around correctly
                j = W;
                i = mod(i,H)+1;
            end
        else % even rows -> west 
            j = j-1;
            if (j < 1) % wrap around correctly
                j = 1; 
                i = mod(i,H)+1;
            end
        end
    end
    i = 1;
    j = 1;
    for ctr = 1:H*W
        %fprintf('Current coordinates are (%d, %d)\n', i, j);
        % get correct ids of horizontal and vertical links
        h_up = H_input(i,j);
        h_down = H_output(i,j);
        v_up = V_input(i,j);
        v_down = V_output(i,j); 
        k1(t+1,h_down)= min(kj, max(k1(t,h_down)+dt*(f1(t,h_down)-g1(t,h_down))/L1, 0));
        k2(t+1,v_down) = min(kj, max(k1(1,h_down) + k2(1,v_down)-k1(t+1, h_down),0));
        if(mod(i,2) == 1) % odd rows -> east 
            j = j+1;
            if (j > W) % wrap around correctly
                j = W;
                i = mod(i,H)+1;
            end
        else % even rows -> west 
            j = j-1;
            if (j < 1) % wrap around correctly
                j = 1; 
                i = mod(i,H)+1;
            end
        end
    end
    if(mod(t-1,NT)==0 && t-1 ~= 0)
%         fprintf('NT is: %.3f\n', NT);
%         fprintf('current iteration is: %.10f\n', t);
        cumQ = ((sum(sum(g1(t-(NT):t,:))/(H*W))/(NT)+ sum(sum(g2(t-(NT):t,:))/(H*W))/(NT)))/2;
%         fprintf('Flow is: %.10f\n', cumQ);
        cumQa = [cumQa cumQ];
    end

%     if(t > NT)
% %         cumQ = (sum(sum(g1(t-NT:t,:))))*dt/(cycle*H*W)*3600;
%         cumQ = ((sum(sum(g1(t-NT:t,:)*dt)/(H*W))/cycle + sum(sum(g2(t-NT:t,:)*dt)/(H*W))/cycle))/2*3600;
%         disp(cumQ);
%     end
%     if(t > NT)
%        cumQ = (sum(g1(t-NT:t,:)*dt) + sum(g2(t-NT:t,:)*dt))/(cycle*H*W)*3600;
%         fprintf('cumq at time t=%f2 is %f2\n', t, cumQ);
%     end
    %fprintf('NT is %f2\n', NT);
end

% disp('DONE_LQMGRIDNETWORK');
% disp('size of g1');
% disp(size(g1(N-NT:N,:)));
%disp(g1(:, :));
cumQ = ((sum(sum(g1(N-NT:N,:)*dt)/(H*W))/cycle + sum(sum(g2(N-NT:N,:)*dt)/(H*W))/cycle))/2*3600;
cumQa = [cumQa cumQ]; 
%cumQ = sum(sum(g1*dt))/(H*W*cycle)*3600;

k1_ret = k1(t+1);
k2_ret = k2(t+1);