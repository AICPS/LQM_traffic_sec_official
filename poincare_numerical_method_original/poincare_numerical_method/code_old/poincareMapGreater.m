function [k12]=poincareMapGreater(k,k1,T,g1,g2,xi1,xi2,kc,kj,vf,L)

gamma1=(1-xi1)*vf/L;
gamma4=(1-xi2)*vf/L;

gamma2=(1-xi1)*vf*kc/L/xi1/(kj-kc);
gamma5=(1-xi2)*vf*kc/L/xi2/(kj-kc);

gamma3=vf*kc/L/(kj-kc);

%% signals
greenT1=(g1*T)/3600;
greenT2=(g2*T)/3600;

%% **************case 1*****************
if(k>(kj+kc)/2) 
    %% for phase 1
    bound1=2*xi1*k-(2*xi1-1)*kj;
    bound2=2*k*(1-xi2)-(1-2*xi2)*kj;
    
    if(k1<2*k-kj) %k1 out of bound
        k1=2*k-kj;
    end
    if(k1>kj)
        k1=kj;
    end
    
    if(k1>bound1) %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT1)+kj;
        if(k1_tmp>=bound1) %still in region 3
            k11=k1_tmp;
        else %go to region 4
            t1=log((kj-bound1)/(kj-k1))/gamma2;
            dt=greenT1-t1;
            k11=(kj-2*k+bound1)*exp(-gamma3*dt)-(kj-2*k);
        end
    else % in region 4
        k11=(kj-2*k+k1)*exp(-gamma3*greenT1)-(kj-2*k);
    end
    
    %% for phase 2    
    if(k11<bound2) %in region 7?
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT2)-(kj-2*k);
        if(k1_tmp<=bound2) %do not cross region 7
            k12=k1_tmp;
        else %go to region 8
            t2=log((kj-2*k+bound2)/(kj-2*k+k11))/gamma5;
            dt=greenT2-t2;
            k12=kj+(bound2-kj)*exp(-gamma3*dt);
        end 
    else %in region 8
        k12=kj+(k11-kj)*exp(-gamma3*greenT2);
    end
  
%% **************case 2*****************
elseif(k<=(kj+kc)/2 && k>kj/2) 
    if(k1<2*k-kj) %k1 out of bound
        k1=2*k-kj;
    end
    if(k1>kj)
        k1=kj;
    end
    bound1=2*k-xi1*kj-(1-xi1)*kc;
    bound2=kj-xi1*(kj-kc);
    
    bound3=2*k-kj+xi2*(kj-kc);
    bound4=kj-(1-xi2)*(kj-kc);
    %% for phase 1
    if(k1<=bound1)%in region 4
        k11=(kj-2*k+k1)*exp(-gamma3*greenT1)-(kj-2*k);
    elseif(k1<=bound2 && k1>bound1) %in region 2  
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % still in region 2
            k11=k1_tmp;
        else % go to region 4
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT1-t1;
            k11=(kj-2*k+bound1)*exp(-gamma3*dt)-(kj-2*k);
        end
    elseif(k1>bound2) %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT1)+kj;
        if(k1_tmp>=bound2) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound2)/(kj-k1))/gamma2;
            dt=greenT1-t1;
            k1_tmp1=bound2-gamma1*kc*dt;
            if(k1_tmp1>=bound1) % still in region 2
                k11=k1_tmp1;
            else % go to region 4
                t2=(bound2-bound1)/gamma1/kc;
                dt1=dt-t2;
                k11=(kj-2*k+bound1)*exp(-gamma3*dt1)-(kj-2*k);
            end
        end    
    end
    
    %% for phase 2
    if(k11>=bound4)%in region 8
        k12=kj+(k11-kj)*exp(-gamma3*greenT2);
    elseif(k11>=bound3 && k11<bound4) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound4) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound4-k11)/gamma4/kc;
            dt=greenT2-t1;
            k12=kj+(bound4-kj)*exp(-gamma3*dt);
        end       
    else %in region 7
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT2)-(kj-2*k);
        if(k1_tmp<=bound3) %still in region 7
            k12=k1_tmp;
        else %go to region 6
            t1=log((kj-2*k+bound3)/(kj-2*k+k11))/gamma5;
            dt=greenT2-t1;
            k1_tmp1=bound3+gamma4*kc*dt;
            if(k1_tmp1<=bound4) %still in region 6
                k12=k1_tmp1;
            else %go to region 8
                t2=(bound4-bound3)/gamma4/kc;
                dt1=dt-t2;
                k12=kj+(bound4-kj)*exp(-gamma3*dt1);
            end
        end 
    
    end

%% more complicated cases
elseif(k<=kj/2 && k>kc/2) 
    if(xi1>xi2)
        [k12]=Xi_1_greater_Xi_2_greater_half(k,k1,T,g1,g2,xi1,xi2,kc,kj,vf,L);
    else
        [k12]=Xi_1_smaller_Xi_2_greater_half(k,k1,T,g1,g2,xi1,xi2,kc,kj,vf,L);
    end

%% **************case 6*****************
elseif(k<=kc/2 && k>=0)
   if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k) % k1 out of bound
        k1=2*k;
    end
    
    %% for phase 1
    k11=k1*exp(-gamma1*greenT1); %region 1

    %% for phase 2
    k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2); %region 5
    
    
end