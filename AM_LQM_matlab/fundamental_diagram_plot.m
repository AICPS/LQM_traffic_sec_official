
clear
clc

%% fundamental diagram settings
vf=60;
kj=150;
kc=30;
C=vf*kc;

%% signal settings
T=90;
L=1;

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
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'ppt'));
width=650;
height=350;

nMax=1000; %maximum iteration number
threshold=10^(-5); %stopping threshold

xi1=[0.3 0.5 0.85 0.3 0.3 0.6  0.4 0.6 0.85 0.3 0.5 0.85 0.3 0.5 0.85];
xi2=[0.3 0.5 0.85 0.4 0.6 0.85 0.3 0.3 0.6  0.3 0.5 0.85 0.3 0.5 0.85];

gT1=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.6 0.6 0.6  0.4 0.4 0.4 ];
gT2=[0.5 0.5 0.5  0.5 0.5 0.5  0.5 0.5 0.5  0.4 0.4 0.4  0.6 0.6 0.6 ];

for i=1:length(gT1)
    j=i;
    %roots
    k13=[];
    k14=[];
    fprintf('gT1 is %f1 and gT2 is %f1\n', gT1(i), gT2(j));
    fprintf('xi1 is %f1 and xi2 is %f1\n', xi1(i), xi2(j));
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
            [k11]=poincareMapPhaseOne(k,k10,T,gT1(i),xi1(j),kc,kj,vf,L);
            [p1]=poincareMapPhaseTwo(k,k11,T,gT2(i),xi2(j),kc,kj,vf,L);
            
            %obtain the corresponding root function value
            g1=k10-p1;
            
            if(g1==0) %it is a root
                k13=[k13;k];
                k14=[k14;k10];
            else %if it is not a root
                k11=p1; %find another point
                
                [k12]=poincareMapPhaseOne(k,k11,T,gT1(i),xi1(j),kc,kj,vf,L);
                [p2]=poincareMapPhaseTwo(k,k12,T,gT2(i),xi2(j),kc,kj,vf,L);
                
                %obtain the corresponding root function value
                g2=k11-p2;
                
                for n=1:nMax %do the iterations
                    if(abs(k11-k10)<threshold || g2==0) %satisfy the requirement
                        k13=[k13;k];
                        k14=[k14;k11];
                        break;
                    end
                    
                    if(g1==g2 && g1~=0)% not the solution
                        break;
                    end
                    %update the values
                    
                    k1n=k11-g2*(k11-k10)/(g2-g1);
                    g1=g2;
                    k10=k11;
                    k11=k1n;
                    
                    [k12]=poincareMapPhaseOne(k,k11,T,gT1(i),xi1(j),kc,kj,vf,L);
                    [p2]=poincareMapPhaseTwo(k,k12,T,gT2(i),xi2(j),kc,kj,vf,L);
                    
                    g2=k11-p2;
                    
                end
            end
        end
    end
    
%     f=figure('Position', [450 400 width height],'Visible','on');
%     
%     plot(zeros(kj/2+1,1),(0:1:kj/2),'k','LineWidth',2)
%     hold on
%     plot(zeros(kj/2+1,1)+kj,(0:1:kj/2)+kj/2,'k','LineWidth',2)
%     
%     plot((0:1:kj),(0:1:kj)/2,'k','LineWidth',2)
%     plot((0:1:kj),(0:1:kj)/2+kj/2,'k','LineWidth',2)
%     plot((0:1:kj),(0:1:kj), '--k')
%     
%     %draw stationary states
%     h1=scatter(k14,k13,'b');
%     
%     %draw regions in phase one
%     xa1=[0:0.1:kc];
%     ya1=kj/2-(((1-xi1(j))*kj-(2-xi1(j))*kc))*xa1/2/kc;
%     h2=plot(xa1,ya1,'--c');
%     
%     ya2=[kc/2:0.1:(xi1(j)*kj+(2-xi1(j))*kc)/2];
%     xa2=kc*ones(size(ya2));
%     plot(xa2,ya2,'--c');
%     
%     ya3=[(kj-xi1(j)*(kj-kc))/2:0.1:(kj+kc)/2];
%     xa3=(kj-xi1(j)*(kj-kc))*ones(size(ya3));
%     plot(xa3,ya3,'--c');
%     
%     xa4=[kc:0.1:(kj-xi1(j)*(kj-kc))];
%     ya4=(xi1(j)*kj+(1-xi1(j))*kc+xa4)/2;
%     plot(xa4,ya4,'--c');
%     
%     xa5=[(kj-xi1(j)*(kj-kc)):0.1:kj];
%     ya5=(2*xi1(j)-1)/2/xi1(j)*kj+1/2/xi1(j)*xa5;
%     plot(xa5,ya5,'--c');
%     
%     %draw regions in phase one
%     xb1=[0:0.1:kj-(1-xi2(j))*(kj-kc)];
%     yb1=(xb1+kc)/2;
%     h3=plot(xb1,yb1,'-.m');
%     
%     xb2=[0:0.1:kj-(1-xi2(j))*(kj-kc)];
%     yb2=(xb1+(kj-xi2(j)*(kj-kc)))/2;
%     plot(xb2,yb2,'-.m');
%     
%     xb3=[kj-(1-xi2(j))*(kj-kc):0.1:kj];
%     yb3=xb3/2+kc*(kj-xb3)/2/(1-xi2(j))/(kj-kc);
%     plot(xb3,yb3,'-.m');
%     
%     xb4=[kj-(1-xi2(j))*(kj-kc):0.1:kj];
%     yb4=((1-2*xi2(j))*kj+xb4)/2/(1-xi2(j));
%     plot(xb4,yb4,'-.m');
%     
%     yb5=[(xi2(j)*kj+(2-xi2(j))*kc)/2:0.1:(kj+kc)/2];
%     xb5=kj-(1-xi2(j))*(kj-kc)*ones(size(yb5));
%     plot(xb5,yb5,'-.m');
%     
%     legend([h1 h2 h3],'Stationary states (fixed points)',...
%         'Region Boundaries in Phase 1','Region Boundaries in Phase 2','Location', 'NorthWest');
%     
%     grid on
%     set(gca,'XTick',[0:kc:kj]);
%     set(gca,'XLim',[0 kj]);
%     set(gca,'YTick',[0:kc:kj]);
%     set(gca,'YLim',[0 kj]);
%     xlabel('k_1','FontSize',13)
%     ylabel('k','FontSize',13)
%     saveas(f,strcat('PD_T', num2str(T),'_g1_',num2str(gT1(i)),'_g2_',num2str(gT2(i)),'_xi1_',num2str(xi1(j)),'_xi2_',num2str(xi2(j)),'.png'));
%     close(f)
    
    %double ring
%     f=figure('Position', [450 400 width height],'Visible','on');
%     [row,col]=size(k14);
%     avgf=zeros(row,col);
%     avgk=k13;
%     for t=1:row
%         fprintf('%d: \n', t);
%         [avgf(t,1),x,y]=flow_calculation(avgk(t),k14(t),xi1(j),xi2(j),T, gT1(i),gT2(i), vf,kj, kc, L);
%     end
%     scatter(avgk,avgf,'b');
%     hold on
%     xlabel('k','FontSize',13)
%     ylabel('q','FontSize',13)
%     axis([0 150 0 1000])
%     grid on
%     saveas(f,strcat('FD_T', num2str(T),'_g1_',num2str(gT1(i)),'_g2_',num2str(gT2(i)),'_xi1_',num2str(xi1(j)),'_xi2_',num2str(xi2(j)),'.png'));
%     close(f)
    
    %gridnetwork
    f=figure('Position', [450 400 width height],'Visible','on');
    [row,col]=size(k14);
    avgf=zeros(row,col);
    avgk=k13;
    for t=1:row
        %fprintf('%d: ', t);
        [avgf(t,1),x,y]=LQM_grid_network(avgk(t),k14(t),xi1(j),xi2(j),T, gT1(i),gT2(i), vf,kj, kc, L);
    end
    scatter(avgk,avgf,'b');
    hold on
    xlabel('k','FontSize',13)
    ylabel('q','FontSize',13)
    axis([0 150 0 1000])
    grid on
    saveas(f,strcat('FD_GN_T', num2str(T),'_g1_',num2str(gT1(i)),'_g2_',num2str(gT2(i)),'_xi1_',num2str(xi1(j)),'_xi2_',num2str(xi2(j)),'.png'));
    close(f)
    
end

