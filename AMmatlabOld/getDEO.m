%% Description
% Can choose to return the first values that satisfy the expected DEO
% Or can choose to return the DEO of given initial k1, k values

function [region1,region2] = getDEO(in_xi1, in_xi2, in_gt1, in_gt2, in_k1, in_k, T)

%clc;

%% fundamental diagram settings
vf=60;
kj=150;
kc=30;
C=vf*kc;

%% signal settings
% T=10;
L=1;

%% Targeted initial states
k1_targ = [];
k_targ = [];

%% Iteration Settings
nMax=1000; %maximum iteration number
threshold=10^(-5); %stopping threshold

if(isempty(in_xi1)==0 && isempty(in_xi2)==0)
    % use the user input
    xi1 = in_xi1;
    xi2 = in_xi2;
else
    % don't use user input
    xi1=[0.3 0.5 0.85 0.3 0.3 0.6  0.4 0.6 0.85 0.3 0.5 0.85 0.3 0.5 0.85];
    xi2=[0.3 0.5 0.85 0.4 0.6 0.85 0.3 0.3 0.6  0.3 0.5 0.85 0.3 0.5 0.85];
end

if(isempty(in_gt1) ==0 && isempty(in_gt2)==0)
    % use the user input
    gT1 = in_gt1;
    gT2 = in_gt2;
else
    %don't use the user input
    gT1=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.6 0.6 0.6  0.4 0.4 0.4 ];
    gT2=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.4 0.4 0.4  0.6 0.6 0.6 ];
end


if(isempty(in_k) == 0)
    % use the user input
    k_arr = in_k;
else
    % don't use the user input
    k_arr = 0:1:kj;
end

for i=1:length(gT1)
    j=i;
    %roots
    k13=[];
    k14=[];
    %fprintf('gT1 is %f2 and gT2 is %f2\n', gT1(i), gT2(j));
    %fprintf('xi1 is %f2 and xi2 is %f2\n', xi1(i), xi2(j));
    %fprintf('length of k_arr is: %d\n', length(k_arr));
    for l=1:length(k_arr) %for each initial density
        %searching regions
        % valid bounds according to average density
        k = k_arr(l);
        kStart=max(2*k-kj,0);
        kEnd=min(kj,2*k);

        if(isempty(in_k1)==0)
            % use the user input
            k1_arr = in_k1;
        else
            % don't use the user input
            k1_arr = kStart:0.2:kEnd;
        end
        
        %start searching
        for m = 1:length(k1_arr)
            %     for k1=kStart:1:kEnd
            %start from the first value
            k10=k1_arr(m);
            % get the regions for (k1(0), k) according to the Poincare Maps 
            [region1, k11] = getRegion1(k,k10,T,gT1(i),xi1(j),kc,kj,vf,L);
            [region2, pk1] = getRegion2(k,k11,T,gT2(i),xi2(j),kc,kj,vf,L);
%             if region1 ==3 && region2 == 8
%                 fprintf('(%d, %d) is the DEO for (k1, k) = (%f2, %f2)\n', region1, region2, k10, k);
%                 k1_targ = [k1_targ k10];
%                 k_targ = [k_targ k];
%             end
        end
    end

end