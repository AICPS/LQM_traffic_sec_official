
clear
clc

%% fundamental diagram settings
vf=60;
kj=150;
kc=30;
C=vf*kc;

%% signal settings
T=100;
delta=0;
xi=0.3;
L=1;

xi_cri=(kj-2*kc)/(kj-kc)/2;
%% draw the function


% kAvg=[20 60 78 120];
% for i=1:4
%     k=kAvg(i);
%     
%     k11=[];
%     k12=[];
%     if(k>kj/2)
%         kStart=2*k-kj;
%         kEnd=kj;
%     else
%         kStart=0;
%         kEnd=2*k;
%     end
%     
%     for k1=kStart:1:kEnd
%         k11=[k11;k1];
%         if(xi>xi_cri)
%             k12_tmp=poincareMapGreater(k,k1,T,delta,xi,kc,kj,vf,L);
%         else
%             k12_tmp=poincareMapSmaller(k,k1,T,delta,xi,kc,kj,vf,L);
%         end
%         k12=[k12;k12_tmp];
%     end
%     figure(1)
%     subplot(2,2,i)
%     plot(k11,k12)
%     hold on
%     plot(k11,k11,'r')
%     
%     xlabel('k_1','FontSize',13)
%     ylabel('P(k_1)','FontSize',13)
%     title(strcat('k=',int2str(k),'vpm'),'FontSize',13)
%     grid on
%     
%     figure(2)
%     subplot(2,2,i)
%     plot(k11,k11-k12,'r','LineWidth',2)
%     hold on
%     plot(k11,zeros(size(k11)),'k')
%     
%     xlabel('k_1','FontSize',15)
%     ylabel('\Phi (k_1)','FontSize',15)
%     title(strcat('k=',int2str(k),'vpm'),'FontSize',13)
%     grid on
%     axis([kStart kEnd min(k11-k12)-0.1 max(k11-k12)+0.1 ])
% end

%% newton method

nMax=1000; %maximum iteration number
threshold=10^(-5); %stopping threshold

k13=[]; %initial values
k14=[];
for k=0:1:kj %for each initial density
    
    %searching regions
    kStart=max(2*k-kj,0);
    kEnd=min(kj,2*k);
  
    %start searching
    for k1=kStart:0.2:kEnd
%     for k1=kStart:1:kEnd
        %start from the first value
        k10=k1;
        
        %obtain the corresponding poincare map value
        [k11]=poincareMapPhaseOne(k,k1,T,g1,xi1,kc,kj,vf,L);
        [k12]=poincareMapPhaseTwo(k,k11,T,g2,xi2,kc,kj,vf,L)
        if(xi>xi_cri)
            p1=poincareMapGreater(k,k10,T,delta,xi,kc,kj,vf,L);
        else
            p1=poincareMapSmaller(k,k10,T,delta,xi,kc,kj,vf,L);
        end
        %obtain the corresponding root function value
        g1=k10-p1;
        
        if(g1==0) %it is a root
            k13=[k13;k];
            k14=[k14;k10];
        else %if it is not a root
            k11=p1; %find another point

            if(xi>xi_cri) %get the corresponding value
                p2=poincareMapGreater(k,k11,T,delta,xi,kc,kj,vf,L);
            else
                p2=poincareMapSmaller(k,k11,T,delta,xi,kc,kj,vf,L);
            end
            %obtain the corresponding root function value
            g2=k11-p2;
            
            for n=1:nMax %do the iterations
                if(abs(k11-k10)<threshold || g2==0) %satisfy the requirement
                    k13=[k13;k];
                    k14=[k14;k11];
                    break;
                end
                
                %update the values
 
                k1n=k11-g2*(k11-k10)/(g2-g1);              
                g1=g2;
                k10=k11;
                k11=k1n;
                
                if(xi>xi_cri)
                    p2=poincareMapGreater(k,k11,T,delta,xi,kc,kj,vf,L);
                else
                    p2=poincareMapSmaller(k,k11,T,delta,xi,kc,kj,vf,L);
                end
                g2=k11-p2;
                
            end
        end
    end    
end

figure(3)
    plot(zeros(kj/2+1,1),(0:1:kj/2),'k','LineWidth',2)
    hold on
    plot(zeros(kj/2+1,1)+kj,(0:1:kj/2)+kj/2,'k','LineWidth',2)
    
    plot((0:1:kj),(0:1:kj)/2,'k','LineWidth',2)
    plot((0:1:kj),(0:1:kj)/2+kj/2,'k','LineWidth',2)
    plot((0:1:kj),(0:1:kj), 'k','LineWidth',2)
    scatter(k14,k13,'b');
    axis([-5 kj+5 -5 kj+5])
    xlabel('k_1','FontSize',13)
    ylabel('k','FontSize',13)
    title(strcat('\xi=',num2str(xi),' & ','T=',int2str(T),'s', ' & ','\Delta=',int2str(delta),'s'),'FontSize',13)
    
%     xlswrite('stationary_state.xlsx', [k14,k13],'stationary_state')
    
figure(4)
    [row,col]=size(k14);
    avgf=zeros(row,col);
    avgk=k13;
    for i=1:row
        avgf(i,1)=flow_calculation(avgk(i,1),k14(i,1),xi,T, delta, vf,kj, kc, L);
    end
    scatter(avgk,avgf,'b');
    hold on
    xlabel('k','FontSize',13)
    ylabel('q','FontSize',13)
    title(strcat('\xi=',num2str(xi),' & ','T=',int2str(T),'s', ' & ','\Delta=',int2str(delta),'s'),'FontSize',13)

    axis([0 150 0 1000])
    grid on

%     xlswrite('stationary_state.xlsx', [avgk,avgf],'flow-rate')

