function [k12]=poincareMapSmaller(k,k1,T,delta,xi,kc,kj,vf,L)

gamma1=(1-xi)*vf/L;
gamma4=gamma1;

gamma2=(1-xi)*vf*kc/L/xi/(kj-kc);
gamma5=gamma2;

gamma3=vf*kc/L/(kj-kc);

%% signals
greenT=(T/2-delta)/3600;

%% **************case 2*****************
if(k>(kj+kc)/2) 
    %% for phase 1
    bound1=2*xi*k-(2*xi-1)*kj;
    bound2=2*k*(1-xi)-(1-2*xi)*kj;
    
    if(k1<2*k-kj) %k1 out of bound
        k1=2*k-kj;
    end
    if(k1>kj)
        k1=kj;
    end
    
    if(k1>bound1) %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT)+kj;
        if(k1_tmp>=bound1) %still in region 3
            k11=k1_tmp;
        else %go to region 4
            t1=log((kj-bound1)/(kj-k1))/gamma2;
            dt=greenT-t1;
            k11=(kj-2*k+bound1)*exp(-gamma3*dt)-(kj-2*k);
        end
    else % in region 4
        k11=(kj-2*k+k1)*exp(-gamma3*greenT)-(kj-2*k);
    end
    
    %% for phase 2    
    if(k11<bound2) %in region 7?
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT)-(kj-2*k);
        if(k1_tmp<=bound2) %do not cross region 7
            k12=k1_tmp;
        else %go to region 8
            t2=log((kj-2*k+bound2)/(kj-2*k+k11))/gamma5;
            dt=greenT-t2;
            k12=kj+(bound2-kj)*exp(-gamma3*dt);
        end 
    else %in region 8
        k12=kj+(k11-kj)*exp(-gamma3*greenT);
    end
  
%% **************case 2*****************
elseif(k<=(kj+kc)/2 && k>kj/2) 
    if(k1<2*k-kj) %k1 out of bound
        k1=2*k-kj;
    end
    if(k1>kj)
        k1=kj;
    end
    bound1=2*k-xi*kj-(1-xi)*kc;
    bound2=kj-xi*(kj-kc);
    bound3=2*k-kj+xi*(kj-kc);
    bound4=kj-(1-xi)*(kj-kc);
    %% for phase 1
    if(k1<=bound1)%in region 4
        k11=(kj-2*k+k1)*exp(-gamma3*greenT)-(kj-2*k);
    elseif(k1<=bound2 && k1>bound1) %in region 2  
        k1_tmp=k1-gamma1*kc*greenT;
        if(k1_tmp>=bound1) % still in region 2
            k11=k1_tmp;
        else % go to region 4
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT-t1;
            k11=(kj-2*k+bound1)*exp(-gamma3*dt)-(kj-2*k);
        end
    elseif(k1>bound2) %in region 3
        k1_tmp=(k1-kj)*exp(gamma2*greenT)+kj;
        if(k1_tmp>=bound2) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound2)/(kj-k1))/gamma2;
            dt=greenT-t1;
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
        k12=kj+(k11-kj)*exp(-gamma3*greenT);
    elseif(k11>=bound3 && k11<bound4) %in region 6
        k1_tmp=k11+gamma4*kc*greenT;
        if(k1_tmp<=bound4) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound4-k11)/gamma4/kc;
            dt=greenT-t1;
            k12=kj+(bound4-kj)*exp(-gamma3*dt);
        end       
    else %in region 7
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT)-(kj-2*k);
        if(k1_tmp<=bound3) %still in region 7
            k12=k1_tmp;
        else %go to region 6
            t1=log((kj-2*k+bound3)/(kj-2*k+k11))/gamma5;
            dt=greenT-t1;
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

%% **************case 3*****************
elseif(k<=kj/2 && k>(kj-xi*(kj-kc))/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=(kj-2*k)*kc/((1-xi)*kj-(2-xi)*kc);
    bound2=2*k-xi*kj-(1-xi)*kc;
    bound3=kj-xi*(kj-kc);
    
    bound4=2*k-kj+xi*(kj-kc);
    bound5=kj-(1-xi)*(kj-kc); 
    bound6=(2*k*(1-xi)*(kj-kc)-kc*kj)/((1-xi)*(kj-kc)-kc);
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT);
    elseif(k1<=bound2 && k1>bound1) %in region 4 
        k1_tmp=(kj-2*k+k1)*exp(-gamma3*greenT)-(kj-2*k);
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=-log((kj-2*k+bound1)/(kj-2*k+k1))/gamma3;
            dt=greenT-t1;
            k11=bound1*exp(-gamma1*dt);
        end  
    elseif(k1<=bound3 && k1>bound2) %in region 2  
        k1_tmp=k1-gamma1*kc*greenT;
        if(k1_tmp>=bound2) % not go to region 4
            k11=k1_tmp;
        else % go to region 4
            t1=(k1-bound2)/gamma1/kc;
            dt=greenT-t1;
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
        k1_tmp=(k1-kj)*exp(gamma2*greenT)+kj;
        if(k1_tmp>=bound3) %still in region 3
            k11=k1_tmp;
        else %go to region 2
            t1=log((kj-bound3)/(kj-k1))/gamma2;
            dt=greenT-t1;
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
        k12=2*k*(1-exp(-gamma4*greenT))+k11*exp(-gamma4*greenT);
    elseif(k11>=bound5 && k11<bound6)%in region 8
        k1_tmp=kj+(k11-kj)*exp(-gamma3*greenT);
        if(k1_tmp<bound6) % still in region 8
            k12=k1_tmp;
        else %go to region 5
            t1=-log((kj-bound6)/(kj-k11))/gamma3;
            dt=greenT-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound6*exp(-gamma4*dt);
        end
    elseif(k11>=bound4 && k11<bound5) %in region 6
        k1_tmp=k11+gamma4*kc*greenT;
        if(k1_tmp<=bound5) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound5-k11)/gamma4/kc;
            dt=greenT-t1;
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
        k1_tmp=(kj-2*k+k11)*exp(gamma5*greenT)-(kj-2*k);
        if(k1_tmp<=bound4) %still in region 7
            k12=k1_tmp;
        else %go to region 6
            t1=log((kj-2*k+bound4)/(kj-2*k+k11))/gamma5;
            dt=greenT-t1;
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
    
   
        %% **************case 4*****************
elseif( k<=(kj-xi*(kj-kc))/2 && k>(xi*kj+(2-xi)*kc)/2)
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    
    bound1=(kj-2*k)*kc/((1-xi)*kj-(2-xi)*kc);
    bound2=2*k-xi*kj-(1-xi)*kc;

    bound3=kj-(1-xi)*(kj-kc); 
    bound4=(2*k*(1-xi)*(kj-kc)-kc*kj)/((1-xi)*(kj-kc)-kc);
    
    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT); 
    elseif(k1<=bound2 && k1>bound1) %in region 4  
        k1_tmp=(kj-2*k+k1)*exp(-gamma3*greenT)-(kj-2*k);
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=-log((((kj-2*k+bound1))/(kj-2*k+k1)))/gamma3;
            dt=greenT-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    else %in region 2
        k1_tmp=k1-gamma1*kc*greenT;
        if(k1_tmp>=bound2) %still in region 2
            k11=k1_tmp;
        else %go to region 4
            t1=(k1-bound2)/gamma1/kc;
            dt=greenT-t1;
            k1_tmp1=(kj-2*k+bound2)*exp(-gamma3*dt)-(kj-2*k);
            if(k1_tmp1>=bound1) % still in region 4
                k11=k1_tmp1;
            else % go to region 1
                t2=-log((((kj-2*k+bound1))/(kj-2*k+bound2)))/gamma3;
                dt1=dt-t2;
                k11=bound1*exp(-gamma1*dt1);
            end
        end    
    end
    
    %% for phase 2
    if(k11>=bound4) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT))+k11*exp(-gamma4*greenT);
    elseif(k11>=bound3 && k11<bound4) %in region 8
        k1_tmp=kj+(k11-kj)*exp(-gamma3*greenT);
        if(k1_tmp<=bound4) %still in region 8
            k12=k1_tmp;
        else %go to region 5
            t1=-log((kj-bound4)/(kj-k11))/gamma3;
            dt=greenT-t1;            
            k12=2*k*(1-exp(-gamma4*dt))+bound4*exp(-gamma4*dt);
        end       
    else %in region 6
        k1_tmp=k11+gamma4*kc*greenT;
        if(k1_tmp<=bound3) %still in region 6
            k12=k1_tmp;
        else %go to region 8
            t1=(bound3-k11)/gamma4/kc;
            dt=greenT-t1;
            k1_tmp1=kj+(bound3-kj)*exp(-gamma3*dt);
            if(k1_tmp1<=bound4) %still in region 8
                k12=k1_tmp1;
            else %go to region 5
                t2=-log((kj-bound4)/(kj-bound3))/gamma3;
                dt1=dt-t2;                
                k12=2*k*(1-exp(-gamma4*dt1))+bound4*exp(-gamma4*dt1);
            end
        end 
    end


    %% **************case 5*****************
elseif(k<=(xi*kj+(2-xi)*kc)/2 && k>kc/2)
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
        k11=k1*exp(-gamma1*greenT); 
    else %in region 2  
        k1_tmp=k1-gamma1*kc*greenT;
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
        else % go to region 1
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    end
    
    %% for phase 2
    if(k11>=bound2) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT))+k11*exp(-gamma4*greenT);
    else %in region 6
        k1_tmp=k11+gamma4*kc*greenT;
        if(k1_tmp<=bound2) %still in region 6
            k12=k1_tmp;
        else %go to region 5
            t1=(bound2-k11)/gamma4/kc;
            dt=greenT-t1;            
            k12=2*k*(1-exp(-gamma4*dt))+bound2*exp(-gamma4*dt);
        end      
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
    k11=k1*exp(-gamma1*greenT); %region 1

    %% for phase 2
    k12=2*k*(1-exp(-gamma4*greenT))+k11*exp(-gamma4*greenT); %region 5
    
    
end