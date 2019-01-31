function [k11]=poincareMapPhaseOne(k,k1,T,g1,xi1,kc,kj,vf,L)

gamma1=(1-xi1)*vf/L;

gamma2=(1-xi1)*vf*kc/L/xi1/(kj-kc);

gamma3=vf*kc/L/(kj-kc);

%% signals
greenT1=(g1*T)/3600;

%% **************case 1*****************
if(k>(kj+kc)/2) 
    %% regions 3 and 4
    %% for phase 1
    bound1=2*xi1*k-(2*xi1-1)*kj;
    
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
      
%% **************case 2*****************
elseif(k<=(kj+kc)/2 && k>kj/2) 
    %% regions 2, 3, 4 
    if(k1<2*k-kj) %k1 out of bound
        k1=2*k-kj;
    end
    if(k1>kj)
        k1=kj;
    end
    bound1=2*k-xi1*kj-(1-xi1)*kc;
    bound2=kj-xi1*(kj-kc);

    %% for phase 1
    if(k1<=bound1)%in region 4
        k11=(kj-2*k+k1)*exp(-gamma3*greenT1)-(kj-2*k);
    elseif(k1<=bound2 && k1>bound1) %in region 2  
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % still in region 2
            k11=k1_tmp;
            %fprintf('k1(0) (%d) is in region 2\n', k1); 
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
                %fprintf('k1(0) (%d) is in region 2\n', k1); 
                k11=k1_tmp1;
            else % go to region 4
                t2=(bound2-bound1)/gamma1/kc;
                dt1=dt-t2;
                k11=(kj-2*k+bound1)*exp(-gamma3*dt1)-(kj-2*k);
            end
        end    
    end
 
%% **************case 3*****************
elseif(k<=kj/2 && k>(xi1*kj+(2-xi1)*kc)/2)
    %% regions 1, 2, 3, 4
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=(kj-2*k)*kc/((1-xi1)*kj-(2-xi1)*kc);
    bound2=2*k-xi1*kj-(1-xi1)*kc;
    bound3=kj-xi1*(kj-kc);
      
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
            %fprintf('k1(0) (%d) is in region 2\n', k1); 
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
                %fprintf('k1(0) (%d) is in region 2\n', k1); 
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
    
%% **************case 4*****************
elseif(k<=(xi1*kj+(2-xi1)*kc)/2 && k>(kj-xi1*(kj-kc))/2)
    %% regions 1, 2, 3 
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=kc;
    bound2=kj-xi1*(kj-kc);

    %% for phase 1
    if(k1<=bound1)%in region 1
        k11=k1*exp(-gamma1*greenT1);
    elseif(k1<=bound2 && k1>bound1) %in region 2
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % not go to region 1
            k11=k1_tmp;
            %fprintf('k1(0) (%d) is in region 2\n', k1); 
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
                %fprintf('k1(0) (%d) is in region 2\n', k1); 
                k11=k1_tmp1;
            else % go to region 1
                t2=(bound2-bound1)/gamma1/kc;
                dt1=dt-t2;
                k11=bound1*exp(-gamma1*dt1);
            end
        end
    end
      
%% **************case 5*****************
elseif(k<=(kj-xi1*(kj-kc))/2 && k>kc/2)
    %% regions 1, 2
    if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k)
        k1=2*k;
    end
    bound1=kc;
    
    %% for phase 1
    if(k1<=bound1)%in region 1 the entire time
        k11=k1*exp(-gamma1*greenT1);
    else %in region 2; must check if k1 goes out of bounds in this phase or not
        k1_tmp=k1-gamma1*kc*greenT1;
        if(k1_tmp>=bound1) % not go to region 1, k1 stays in region 2
            %fprintf('k1(0) (%f2) is in region 2\n', k1); 
            k11=k1_tmp;
        else % go to region 1 for the remaining amount of time in the phase
            t1=(k1-bound1)/gamma1/kc;
            dt=greenT1-t1;
            k11=bound1*exp(-gamma1*dt);
        end
    end
 
%% **************case 6*****************
elseif(k<=kc/2 && k>=0)
    %% region 1
   if(k1<0) %k1 out of bound
        k1=0;
    end
    if(k1>2*k) % k1 out of bound
        k1=2*k;
    end
    
    %% for phase 1
    k11=k1*exp(-gamma1*greenT1); %region 1

end