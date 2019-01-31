function [cumQ, k1, k2]=flow_calculation(k,k10,xi1,xi2,cycle,gT1,gT2, vf,kj, kc, L1)

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

N=NT;
M=1;


%initialization
f1=zeros(N,M);
f2=zeros(N,M);
g1=zeros(N,M);
g2=zeros(N,M);
d1=zeros(N,M);
d2=zeros(N,M);
s1=zeros(N,M);
s2=zeros(N,M);

k1=zeros(N+1,M);
k2=zeros(N+1,M);


k1_tem=zeros(N,M);
k2_tem=zeros(N,M);

T=zeros(N+1,M);

% fprintf('k: %f\n', k);
k1(1,:)=k10;
k2(1,:)=2*k-k10;
% fprintf('k1: %f\n', k1(1,:));
% fprintf('k2: %f\n', k2(1,:));
        
j=1;

for i=1:N
    C=kc*vf;
    d1(i,j)=min(vf*k1(i,j),C);
    d2(i,j)=min(vf*k2(i,j),C);
    s1(i,j)=min(C,C*(kj-k1(i,j))/(kj-kc));
    s2(i,j)=min(C,C*(kj-k2(i,j))/(kj-kc));
%     fprintf('k1: %.1f\n', k1(i,j));
%     fprintf('k2: %.1f\n', k2(i,j));
%     fprintf('d1: %.2f\n', d1(i,j));
%     fprintf('s1: %.2f\n', s1(i,j));
%     fprintf('d2: %.2f\n', d2(i,j));
%     fprintf('s2: %.2f\n', s2(i,j));
    k1_tem(i,j)=k1(i,j);
    k2_tem(i,j)=k2(i,j);
    T(i+1,j)=dt*i;
    
    if(mod(i,NT)<=nGreen1)
        signal=1;
    elseif(mod(i,NT)>nGreen1+nLostTime && mod(i,NT)<=NT-nLostTime)
        signal=2;
    else
        signal=0;
    end
%     fprintf('xi1 is %.2f\n', xi1);
%     fprintf('xi2 is %.2f\n', xi2);
    if(signal==1)%ring 1 goes first
        g1(i,j)=min(d1(i,j),min(s1(i,j)/xi1,s2(i,j)/(1-xi1)));
        g2(i,j)=0;
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
%         fprintf('phase 1\n');
%         fprintf('g1: %.2f\n', g1(i,j));
%         fprintf('f1: %.2f\n', f1(i,j));
%         fprintf('f2: %.2f\n', f2(i,j));
    elseif(signal==2)%ring 2 has green time
        g1(i,j)=0;
        g2(i,j)=min(d2(i,j),min(s2(i,j)/xi2,s1(i,j)/(1-xi2)));
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
%         fprintf('phase 2\n');
%         fprintf('g2: %.2f\n', g2(i,j));
%         fprintf('f1: %.2f\n', f1(i,j));
%         fprintf('f2: %.2f\n', f2(i,j));
    else %during the lost times SUPER IMPORTANT
        g1(i,j)=0;
        g2(i,j)=0;
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    end
    
    k1(i+1,j)=k1(i,j)+dt*(f1(i,j)-g1(i,j))/L1;
    k2(i+1,j)=k1(1,j)+k2(1,j)-k1(i+1,j);
end

% total outflow of g1 over the entire cycle 

%cumQ=sum(g1*dt)/cycle*3600; % dt is in terms of hours ; cycle is in terms
%of seconds? hours? 
cumQ = (sum(g1*dt)/cycle + sum(g2*dt)/cycle)/2*3600;
k1 = k1(i+1,j) ;
k2 = k2(i+1,j);
%fprintf('Cumulative flow (only g1) is: %f\n', cumQ);
% fprintf('NT is: %f\n', NT);%1 minute for each cycle
% fprintf('nlostime is: %f\n', nLostTime); %total steps for the lost time in each phase
% fprintf('ngreen is: %f\n', nGreen1);
% fprintf('N is: %f\n', N);