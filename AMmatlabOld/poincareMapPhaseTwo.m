function [k12]=poincareMapPhaseTwo(k,k11,T,g2,xi2,kc,kj,vf,L)

gamma4=(1-xi2)*vf/L;

gamma5=(1-xi2)*vf*kc/L/xi2/(kj-kc);

gamma3=vf*kc/L/(kj-kc);

%% signals
greenT2=(g2*T)/3600;

%% **************case 1*****************
if(k>(kj+kc)/2) 
    %% regions 7, 8
    bound2=2*k*(1-xi2)-(1-2*xi2)*kj;

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
    %% regions 6, 7, 8
    bound3=2*k-kj+xi2*(kj-kc);
    bound4=kj-(1-xi2)*(kj-kc);
        
    %% for phase 2
    if(k11>=bound4)%in region 8
        k12=kj+(k11-kj)*exp(-gamma3*greenT2);
    elseif(k11>=bound3 && k11<bound4) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound4) %still in region 6
            %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
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
                %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
                k12=k1_tmp1;
            else %go to region 8
                t2=(bound4-bound3)/gamma4/kc;
                dt1=dt-t2;
                k12=kj+(bound4-kj)*exp(-gamma3*dt1);
            end
        end 
    
    end

%% **************case 3*****************
elseif(k<=kj/2 && k>(xi2*kj+(2-xi2)*kc)/2)
    %% regions 5, 6, 7, 8
    
    bound4=2*k-kj+xi2*(kj-kc);
    bound5=kj-(1-xi2)*(kj-kc);
    bound6=(2*k*(1-xi2)*(kj-kc)-kc*kj)/((1-xi2)*(kj-kc)-kc);

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
            %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
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
                %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
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
elseif(k<=(xi2*kj+(2-xi2)*kc)/2 && k>(kj-xi2*(kj-kc))/2)
    %% regions 5, 6, 7
 
    bound3=2*k-kj+xi2*(kj-kc);
    bound4=2*k-kc;
    
    %% for phase 2
    if(k11>=bound4) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    elseif(k11>=bound3 && k11<bound4) %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound4) %still in region 6
            %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
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
                %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
            else %go to region 5
                t2=(bound4-bound3)/gamma4/kc;
                dt1=dt-t2;
                k12=2*k*(1-exp(-gamma4*dt1))+bound4*exp(-gamma4*dt1);
            end
        end
    end

    
%% **************case 5*****************
elseif(k<=(kj-xi2*(kj-kc))/2 && k>kc/2)
    %% regions 5, 6
    bound2=2*k-kc;
   
    %% for phase 2
    if(k11>=bound2) %in region 5
        k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2);
    else %in region 6
        k1_tmp=k11+gamma4*kc*greenT2;
        if(k1_tmp<=bound2) %still in region 6
            %fprintf('k1(t+T) (%d) is in region 6\n', k11); 
            k12=k1_tmp;
        else %go to region 5
            t1=(bound2-k11)/gamma4/kc;
            dt=greenT2-t1;
            k12=2*k*(1-exp(-gamma4*dt))+bound2*exp(-gamma4*dt);
        end
    end

%% **************case 6*****************
elseif(k<=kc/2 && k>=0)
    %% region 5
    %% for phase 2
    k12=2*k*(1-exp(-gamma4*greenT2))+k11*exp(-gamma4*greenT2); %region 5
    
    
end