function [cumQ, q1_d, q2_d, q_d, k1_d]=flow_calculation_GTR_ATTACK(k,k10,xi1,xi2,cycle,gT1,gT2, vf, kj, kc, L1, gT1_d, gT2_d, Tstart, Tend)

%simulation settings
ddt=0.05;
dt=ddt/3600;%2 seconds simulation time step;

%signal settings
NT=cycle*(1/ddt);%1 minute for each cycle
nLostTime=(1-gT1-gT2)*cycle/2/ddt; %total steps for the lost time in each phase
nGreen1=(gT1*cycle)/ddt;

N=NT;
M=1;

%initialization
f1=zeros(N,M);
f2=zeros(N,M);
g1=zeros(N,M);
g2=zeros(N,M);
q_d=zeros(N,M);
d1=zeros(N,M);
d2=zeros(N,M);
s1=zeros(N,M);
s2=zeros(N,M);

k1=zeros(N+1,M);
k2=zeros(N+1,M);


k1_tem=zeros(N,M);
k2_tem=zeros(N,M);

T=zeros(N+1,M);

k1(1,:)=k10;
k2(1,:)=2*k-k10;
        
j=1;

for i=1:N
    C=kc*vf;
    if(i >= Tstart*(1/ddt) && i <= Tend*(1/ddt))
        nLostTime=(1-gT1_d-gT2_d)*cycle/2/ddt; %total steps for the lost time in each phase
        nGreen1=(gT1_d*cycle)/ddt;
    else
        nLostTime=(1-gT1-gT2)*cycle/2/ddt; %total steps for the lost time in each phase
        nGreen1=(gT1*cycle)/ddt;
    end
    d1(i,j)=min(vf*k1(i,j),C);
    d2(i,j)=min(vf*k2(i,j),C);
    s1(i,j)=min(C,C*(kj-k1(i,j))/(kj-kc));
    s2(i,j)=min(C,C*(kj-k2(i,j))/(kj-kc));

    k1_tem(i,j)=k1(i,j);
    k2_tem(i,j)=k2(i,j);
    T(i+1,j)=dt*i; %increment simulation time by dt


    if(mod(i,NT)<=nGreen1) %in phase 1
        signal=1;
        disp('1');
    elseif(mod(i,NT)>nGreen1+nLostTime && mod(i,NT)<=NT-nLostTime)  % in phase 2
        signal=2;
        disp('2');
    else
        signal=0;
    end

    if(signal==1)%ring 1 goes first
        g1(i,j)=min(d1(i,j),min(s1(i,j)/xi1,s2(i,j)/(1-xi1)));
        g2(i,j)=0;
        if (g1(i,j) < 0 ) 
            fprintf('1 g1 is: %f4\n', g1(i,j));
        end
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    elseif(signal==2)%ring 2 has green time
        g1(i,j)=0;
        g2(i,j)=min(d2(i,j),min(s2(i,j)/xi2,s1(i,j)/(1-xi2)));
        if (g1(i,j) < 0 ) 
            fprintf('2 g1 is: %f4\n', g1(i,j));
        end
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    else %during the lost times SUPER IMPORTANT
        g1(i,j)=0;
        g2(i,j)=0;
        if (g1(i,j) < 0 ) 
            fprintf('3 g1 is: %f4\n', g1(i,j));
        end
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    end
    q_d(i,j) = ((sum(g1*dt)*3600  + sum(g2*dt)*3600))/2;
%     if (g1(i,j) < 0 ) 
%         fprintf('g1 is: %f4\n', g1(i,j));
%     end
    k1(i+1,j)=k1(i,j)+dt*(f1(i,j)-g1(i,j))/L1;
    k2(i+1,j)=k1(1,j)+k2(1,j)-k1(i+1,j);
end

% total outflow of g1 over all the cycles
cumQ=(sum(g1*dt)+ sum(g2*dt))/(2*cycle)*3600;
cumQ2 = sum(g1*dt)/cycle*3600;
fprintf('Cumulative flow (only g1) is: %f4\n', cumQ2);
%q_d = sum(g1*dt)/3600  + sum(g2*dt)/3600;
%q_d = q_d/2;
q1_d = g1;
q2_d = g2;
k1_d = k1;



