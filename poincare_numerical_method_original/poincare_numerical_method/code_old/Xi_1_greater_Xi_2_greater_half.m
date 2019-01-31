function [k12]=Xi_1_greater_Xi_2_greater_half(k,k1,T,g1,g2,xi1,xi2,kc,kj,vf,L)

%% consider the case xi1>xi2
gamma1=(1-xi1)*vf/L;
gamma4=(1-xi2)*vf/L;

gamma2=(1-xi1)*vf*kc/L/xi1/(kj-kc);
gamma5=(1-xi2)*vf*kc/L/xi2/(kj-kc);

gamma3=vf*kc/L/(kj-kc);

%% signals
greenT1=(g1*T)/3600;
greenT2=(g2*T)/3600;
%% **************case 1*****************
if(k<=kj/2 && k>(xi1*kj+(2-xi1)*kc)/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=(kj-2*k)*kc/((1-xi1)*kj-(2-xi1)*kc);
    bound2=2*k-xi1*kj-(1-xi1)*kc;
    bound3=kj-xi1*(kj-kc);
    
    bound4=2*k-kj+xi2*(kj-kc);
    bound5=kj-(1-xi2)*(kj-kc);
    bound6=(2*k*(1-xi2)*(kj-kc)-kc*kj)/((1-xi2)*(kj-kc)-kc);
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT1);
    elseif(k1<=bound2 && k1>bound1) %in region 4
        k1_tmp=(kj-2*k+k1)*exp(-gamma3*greenT1)-(kj-2*k);
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=-log((kj-2*k+bound1)/(kj-2*k+k1))/gamma3;
            dt=greenT1-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    elseif(k1<=bound3 && k1>bound2) %in region 2
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound2) % not go to region 4
            k11=k1_tmp;
        else % go to region 4
            t1=(k1-bound2)/gamma1/kc;
            dt=greenT1-t1;
            k1_tmp1=(kj-2*k+bound2)*exp(-gamma3*dt)-(kj-2*k);
            if(k1_tmp1>=bound1) % still in region 4
                k11=k1_tmp1;
            else %go to region 1
                t2=-log((kj-2*k+bound1)/(kj-2*k+bound2))/gamma3;
                dt1=dt-t2;
                k11=bound1*exp(-gamma1*dt1);
            end
        end
    elseif(k1>bound3) %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT1)+kj;
        if(k1_tmp>=bound3) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound3)/(kj-k1))/gamma2;
            dt=greenT1-t1;
            k1_tmp1=bound3-gamma1*kc*dt;
            if(k1_tmp1>=bound2) % still in region 2
                k11=k1_tmp1;
            else % go to region 4
                t2=(bound3-bound2)/gamma1/kc;
                dt1=dt-t2;
                k1_tmp2=(kj-2*k+bound2)*exp(-gamma3*dt1)-(kj-2*k);
                if(k1_tmp2>=bound1) %still in region 4
                    k11=k1_tmp2;
                else % go to region 1
                    t3=-log((kj-2*k+bound1)/(kj-2*k+bound2))/gamma3;
                    dt2=dt1-t3;
                    k11=bound1*exp(-gamma1*dt2);
                end
            end
        end
    end
    
    %% for phase 2
    if(k11>=bound6) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    elseif(k11>=bound5 && k11<bound6)%in region 8
        k1_tmp=kj+(k11-kj)*exp(-gamma3*greenT2);
        if(k1_tmp<bound6) % still in region 8
            k12=k1_tmp;
        else %go to region 5
            t1=-log((kj-bound6)/(kj-k11))/gamma3;
            dt=greenT2-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound6*exp(-gamma4*dt);
        end
    elseif(k11>=bound4 && k11<bound5) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound5) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound5-k11)/gamma4/kc;
            dt=greenT2-t1;
            k1_tmp1=kj+(bound5-kj)*exp(-gamma3*dt);
            if(k1_tmp1<=bound6) %still in region 8
                k12=k1_tmp1;
            else %go to region 5
                t2=-log((kj-bound6)/(kj-bound5))/gamma3;
                dt1=dt-t2;
                k12=2*k*(1-exp(-gamma4*dt1))+bound6*exp(-gamma4*dt1);
            end
        end
    else %in region 7
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT2)-(kj-2*k);
        if(k1_tmp<=bound4) %still in region 7
            k12=k1_tmp;
        else %go to region 6
            t1=log((kj-2*k+bound4)/(kj-2*k+k11))/gamma5;
            dt=greenT2-t1;
            k1_tmp1=bound4+gamma4*kc*dt;
            if(k1_tmp1<=bound5) %still in region 6
                k12=k1_tmp1;
            else %go to region 8
                t2=(bound5-bound4)/gamma4/kc;
                dt1=dt-t2;
                k1_tmp2=kj+(bound5-kj)*exp(-gamma3*dt1);
                if(k1_tmp2<=bound6) %still in region 8
                    k12=k1_tmp2;
                else %go to region 5
                    t3=-log((kj-bound6)/(kj-bound5))/gamma3;
                    dt2=dt1-t3;
                    k12=2*k*(1-exp(-gamma4*dt2))+bound6*exp(-gamma4*dt2);
                end
            end
        end
    end
    
%% **************case 2*****************
elseif(k<=(xi1*kj+(2-xi1)*kc)/2 && k>(xi2*kj+(2-xi2)*kc)/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=kc;
    bound2=kj-xi1*(kj-kc);
    
    bound4=2*k-kj+xi2*(kj-kc);
    bound5=kj-(1-xi2)*(kj-kc);
    bound6=(2*k*(1-xi2)*(kj-kc)-kc*kj)/((1-xi2)*(kj-kc)-kc);
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT1);
    elseif(k1<=bound2 && k1>bound1) %in region 2
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT1-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    else %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT1)+kj;
        if(k1_tmp>=bound2) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound2)/(kj-k1))/gamma2;
            dt=greenT1-t1;
            k1_tmp1=bound2-gamma1*kc*dt;
            if(k1_tmp1>=bound1) % still in region 2
                k11=k1_tmp1;
            else % go to region 1
                t2=(bound2-bound1)/gamma1/kc;
                dt1=dt-t2;
                k11=bound1*exp(-gamma1*dt1);
            end
        end
    end
    
   %% for phase 2
    if(k11>=bound6) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    elseif(k11>=bound5 && k11<bound6)%in region 8
        k1_tmp=kj+(k11-kj)*exp(-gamma3*greenT2);
        if(k1_tmp<bound6) % still in region 8
            k12=k1_tmp;
        else %go to region 5
            t1=-log((kj-bound6)/(kj-k11))/gamma3;
            dt=greenT2-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound6*exp(-gamma4*dt);
        end
    elseif(k11>=bound4 && k11<bound5) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound5) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound5-k11)/gamma4/kc;
            dt=greenT2-t1;
            k1_tmp1=kj+(bound5-kj)*exp(-gamma3*dt);
            if(k1_tmp1<=bound6) %still in region 8
                k12=k1_tmp1;
            else %go to region 5
                t2=-log((kj-bound6)/(kj-bound5))/gamma3;
                dt1=dt-t2;
                k12=2*k*(1-exp(-gamma4*dt1))+bound6*exp(-gamma4*dt1);
            end
        end
    else %in region 7
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT2)-(kj-2*k);
        if(k1_tmp<=bound4) %still in region 7
            k12=k1_tmp;
        else %go to region 6
            t1=log((kj-2*k+bound4)/(kj-2*k+k11))/gamma5;
            dt=greenT2-t1;
            k1_tmp1=bound4+gamma4*kc*dt;
            if(k1_tmp1<=bound5) %still in region 6
                k12=k1_tmp1;
            else %go to region 8
                t2=(bound5-bound4)/gamma4/kc;
                dt1=dt-t2;
                k1_tmp2=kj+(bound5-kj)*exp(-gamma3*dt1);
                if(k1_tmp2<=bound6) %still in region 8
                    k12=k1_tmp2;
                else %go to region 5
                    t3=-log((kj-bound6)/(kj-bound5))/gamma3;
                    dt2=dt1-t3;
                    k12=2*k*(1-exp(-gamma4*dt2))+bound6*exp(-gamma4*dt2);
                end
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
%% **************case 3*****************
elseif(k<=(xi2*kj+(2-xi2)*kc)/2 && k>(kj-xi2*(kj-kc))/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=kc;
    bound2=kj-xi1*(kj-kc);
    
    bound3=2*k-kj+xi2*(kj-kc);
    bound4=2*k-kc;
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT1);
    elseif(k1<=bound2 && k1>bound1) %in region 2
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT1-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    else %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT1)+kj;
        if(k1_tmp>=bound2) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound2)/(kj-k1))/gamma2;
            dt=greenT1-t1;
            k1_tmp1=bound2-gamma1*kc*dt;
            if(k1_tmp1>=bound1) % still in region 2
                k11=k1_tmp1;
            else % go to region 1
                t2=(bound2-bound1)/gamma1/kc;
                dt1=dt-t2;
                k11=bound1*exp(-gamma1*dt1);
            end
        end
    end
    
    %% for phase 2
    if(k11>=bound4) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    elseif(k11>=bound3 && k11<bound4) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound4) %still in region 6
            k12=k1_tmp;
        else %go to region 5
            t1=(bound4-k11)/gamma4/kc;
            dt=greenT2-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound4*exp(-gamma4*dt);
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
            else %go to region 5
                t2=(bound4-bound3)/gamma4/kc;
                dt1=dt-t2;
                k12=2*k*(1-exp(-gamma4*dt1))+bound4*exp(-gamma4*dt1);
            end
        end
    end
    
%% **************case 4*****************
elseif(k<=(kj-xi2*(kj-kc))/2 && k>(kj-xi1*(kj-kc))/2)
    
    
%% **************case 5*****************
elseif(k<=(kj-xi1*(kj-kc))/2 && k>kc/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=kc;
    bound2=2*k-kc;
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT1);
    else %in region 2
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT1-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    end
    
    %% for phase 2
    if(k11>=bound2) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    else %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound2) %still in region 6
            k12=k1_tmp;
        else %go to region 5
            t1=(bound2-k11)/gamma4/kc;
            dt=greenT2-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound2*exp(-gamma4*dt);
        end
    end
end
    