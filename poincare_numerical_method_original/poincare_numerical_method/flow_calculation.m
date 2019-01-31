function [cumQ]=flow_calculation(k,k10,xi1,xi2,cycle,gT1,gT2, vf,kj, kc, L1)

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
    d1(i,j)=min(vf*k1(i,j),C);
    d2(i,j)=min(vf*k2(i,j),C);
    s1(i,j)=min(C,C*(kj-k1(i,j))/(kj-kc));
    s2(i,j)=min(C,C*(kj-k2(i,j))/(kj-kc));
    
    k1_tem(i,j)=k1(i,j);
    k2_tem(i,j)=k2(i,j);
    T(i+1,j)=dt*i;
    
    
    if(mod(i,NT)<=nGreen1)
        signal=1;
    elseif(mod(i,NT)>nGreen1+nLostTime &&mod(i,NT)<=NT-nLostTime)
        signal=2;
    else
        signal=0;
    end
    
    if(signal==1)%ring 1 goes first
        g1(i,j)=min(d1(i,j),min(s1(i,j)/xi1,s2(i,j)/(1-xi1)));
        g2(i,j)=0;
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    elseif(signal==2)%ring 2 has green time
        g1(i,j)=0;
        g2(i,j)=min(d2(i,j),min(s2(i,j)/xi2,s1(i,j)/(1-xi2)));
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    else %during the lost times
        g1(i,j)=0;
        g2(i,j)=0;
        f1(i,j)=g1(i,j)*xi1+g2(i,j)*(1-xi2);
        f2(i,j)=g1(i,j)*(1-xi1)+g2(i,j)*xi2;
    end
    
    k1(i+1,j)=k1(i,j)+dt*(f1(i,j)-g1(i,j))/L1;
    k2(i+1,j)=k1(1,j)+k2(1,j)-k1(i+1,j);
end
       
cumQ=sum(g1*dt)/cycle*3600;
