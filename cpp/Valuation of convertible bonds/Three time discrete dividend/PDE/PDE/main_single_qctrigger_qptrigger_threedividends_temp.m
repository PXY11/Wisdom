% function U = convertible_code(Nt,Ns,F,conversion_price,Smax,S0,risk_free_rate,rate_T,sigma,T,final_coupon_rate,q,p,R,eta,no_call_time,...
% no_put_time,no_convert_time,Bc_star,Bp_star,theta,coupon_time_rate)

%function deltaP=main_single_qctrigger_qptrigger_test(sigma)

%% ¡ª¡ª¡ª¡ª¡ªSolving pde-Penalty mehtod-CN difference-Newton iteration¡ª¡ª¡ª¡ª¡ª¡ª
clc;
clear;
tic;
format long;
%%

fid=fopen('MyFile.txt','w');

%%
X=5.73;
oas=0;
r0=0.035;
r=r0+oas;     %128011
no_convert_time=[5.5 6];
no_call_time=[5.5 6];
no_put_time=[2 6];
% for three different cases, T changes will influence the Shistory length
Shistory=[6.25 6.13 6.23 6.03 6.19 6.32 6.39 6.33 6.20 6.04 6.13 6.12 6.11 6.25 6.19 6.09 6.01 6.13 6.13 6.24 6.33 6.80 7.48 7.70 8.30 7.75 7.77 7.81 7.61];
F=100; % face value
q=0; % dividend
p=0; % default probobility
R=1; % recovery rate
eta=0; % when default the stock price becomes (1-eta)*S
theta=1; % implicitness parameter£¬0.5ÎªCN method,1Îªfully implicit or backward euler method
n1=15; m1=30; ncX=1.3; n2=30; m2=30; npX=0.7;
coupon_rate=[1.8 1.5 1.5 1 0.7 0.5]/100;
bonus_rate=4.2/100;

for case_T=3 % corresponding to three different cases
    
    % three different cases for Shistory
    if case_T == 1
        T=5.7;
    end
    if case_T == 2
        T=5.45;
    end
    if case_T == 3
        T=6;
    end
    
    for case_D = 1:5
        switch case_D
            case 1
                D1 = 0.3;
                D2 = 0.3;
                D3 = 0.3;
            case 2
                D1 = 0.4;
                D2 = 0.4;
                D3 = 0.4;
            case 3
                D1 = 0.5;
                D2 = 0.5;
                D3 = 0.5;
            case 4
                D1 = 1.0;
                D2 = 1.0;
                D3 = 1.0;
            case 5
                D1 = 1.5;
                D2 = 1.5;
                D3 = 1.5;
            case 6
                D1 = 2.0;
                D2 = 2.0;
                D3 = 2.0;
        end
        
        for case_Dtime=1:3
            
            switch case_Dtime
            case 1
                D1_time = 1.0;
                D2_time = 3.0;
                D3_time = 5.0;
            case 2
                D1_time = 3.4;
                D2_time = 3.5;
                D3_time = 3.6;
            case 3
                D1_time = 3.49;
                D2_time = 3.5;
                D3_time = 3.51;
            end
            
            for case_sigma=1:7
                switch case_sigma
                    case 1
                        sigma = 0.2;
                    case 2
                        sigma = 0.4;
                    case 3
                        sigma = 0.6;
                    case 4
                        sigma = 0.8;
                    case 5
                        sigma = 1.0;
                    case 6
                        sigma = 1.5;
                    case 7
                        sigma = 2.5;
                end
                
                for case_S0 = 1:6
                    switch case_S0
                        case 1
                            S0 = X*1.28;
                        case 2
                            S0 = X*1.30;
                        case 3
                            S0 = X*1.32;
                        case 4
                            S0 = X*0.68;
                        case 5
                            S0 = X*0.70;
                        case 6
                            S0 = X*0.72;
                    end
                    
                    for method=1:2
                        %%
                        if method == 2 % Probability method
                            tic;
                            fprintf(fid, 'Probability method \n');
                            fprintf(fid, 'S0 = %f\n', S0');
                            fprintf(fid, 'T = %f\n', T');
                            fprintf(fid, 'D1 = %f\n', D1');
                            fprintf(fid, 'D2 = %f\n', D2');
                            fprintf(fid, 'D3 = %f\n', D3');
                            fprintf(fid, 'D1_time = %f\n', D1_time');
                            fprintf(fid, 'D2_time = %f\n', D2_time');
                            fprintf(fid, 'D3_time = %f\n', D3_time');
                            fprintf(fid, 'sigma = %f\n', sigma');
                            % equivalent call trigger
                            if n1>0
                                ctrigger=X*ncX;
                                qctrigger1_fun=@(qctrigger) calculateqctrigger_threedividends_part1_fixedr_test2(S0, Shistory, ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1,D1_time);
                                qctrigger2_fun=@(qctrigger) calculateqctrigger_threedividends_part2_fixedr_test2(S0, Shistory, ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time);
                                qctrigger3_fun=@(qctrigger) calculateqctrigger_threedividends_part3_fixedr_test2(S0, Shistory, ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time,D2,D3_time);
                                qctrigger4_fun=@(qctrigger) calculateqctrigger_threedividends_part4_fixedr_test2(S0, Shistory, ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time,D2,D3_time,D3);
                                qctrigger0=ctrigger*1.08;
                                % fsolve
                                % qctrigger1=fsolve(qctrigger1_fun, qctrigger0);
                                % qctrigger2=fsolve(qctrigger2_fun, qctrigger0);
                                % qctrigger3=fsolve(qctrigger3_fun, qctrigger0);
                                % qctrigger4=fsolve(qctrigger4_fun, qctrigger0);
                                % fzero
                                options=optimset('Display','iter','TolFun',0.001);
                                qctrigger1=fzero(qctrigger1_fun,qctrigger0,options);
                                qctrigger2=fzero(qctrigger2_fun,qctrigger0,options);
                                qctrigger3=fzero(qctrigger3_fun,qctrigger0,options);
                                qctrigger4=fzero(qctrigger4_fun,qctrigger0,options);
                            else
                                qctrigger1=S0*5;
                                qctrigger2=S0*5;
                                qctrigger3=S0*5;
                                qctrigger4=S0*5;
                            end
                            % equivalent put trigger
                            if n2>0 && (T-D1_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger1_fun=@(qptrigger) calculateqptrigger_threedividends_part1_fixedr_test2(S0, Shistory, ptrigger,qptrigger,T,r,sigma,no_put_time,n2,m2,D1_time);
                                qptrigger0=ptrigger*0.9;
                                % fsolve
                                % qptrigger1=fsolve(qptrigger1_fun, qptrigger0);
                                % fzero
                                options=optimset('Display','iter','TolFun',0.001);
                                qptrigger1=fzero(qptrigger1_fun,qptrigger0,options);
                                if qptrigger1==qptrigger0
                                    qptrigger1=0;
                                end
                            else
                                qptrigger1=0;
                            end
                            
                            if n2>0 && (T-D2_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger2_fun=@(qptrigger) calculateqptrigger_threedividends_part2_fixedr_test2(S0, Shistory, ptrigger,qptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time);
                                qptrigger0=ptrigger*0.9;
                                % fsolve
                                % qptrigger2=fsolve(qptrigger2_fun, qptrigger0);
                                % fzero
                                options=optimset('Display','iter','TolFun',0.001);
                                qptrigger2=fzero(qptrigger2_fun,qptrigger0,options);
                                if qptrigger2==qptrigger0
                                    qptrigger2=0;
                                end
                            else
                                qptrigger2=0;
                            end
                            
                            if n2>0 && (T-D3_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger3_fun=@(qptrigger) calculateqptrigger_threedividends_part3_fixedr_test2(S0, Shistory, ptrigger,qptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time,D2,D3_time);
                                qptrigger0=ptrigger*0.9;
                                % fsolve
                                % qptrigger3=fsolve(qptrigger3_fun, qptrigger0);
                                % fzero
                                options=optimset('Display','iter','TolFun',0.001);
                                qptrigger3=fzero(qptrigger3_fun,qptrigger0,options);
                                if qptrigger3==qptrigger0
                                    qptrigger3=0;
                                end
                            else
                                qptrigger3=0;
                            end
                            
                            if n2>0
                                ptrigger=X*npX;
                                qptrigger4_fun=@(qptrigger) calculateqptrigger_threedividends_part4_fixedr_test2(S0, Shistory, ptrigger,qptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time,D2,D3_time,D3);
                                qptrigger0=ptrigger*0.9;
                                % fsolve
                                % qptrigger4=fsolve(qptrigger4_fun, qptrigger0);
                                % fzero
                                options=optimset('Display','iter','TolFun',0.001);
                                qptrigger4=fzero(qptrigger4_fun,qptrigger0,options);
                                if qptrigger4==qptrigger0
                                    qptrigger4=0;
                                end
                            else
                                qptrigger4=0;
                            end
                        end
                        if method == 1 % Mean method
                            tic;
                            fprintf(fid, 'Mean method \n');
                            fprintf(fid, 'S0 = %f\n', S0');
                            fprintf(fid, 'T = %f\n', T');
                            fprintf(fid, 'D1 = %f\n', D1');
                            fprintf(fid, 'D2 = %f\n', D2');
                            fprintf(fid, 'D3 = %f\n', D3');
                            fprintf(fid, 'D1_time = %f\n', D1_time');
                            fprintf(fid, 'D2_time = %f\n', D2_time');
                            fprintf(fid, 'D3_time = %f\n', D3_time');
                            fprintf(fid, 'sigma = %f\n', sigma');
                            % equivalent call trigger
                            if n1>0
                                ctrigger=X*ncX;
                                qctrigger1=calculateqctrigger_threedividends_part1_fixedr(S0, Shistory, ctrigger,T,r,sigma,no_call_time,n1,m1,D1_time);
                                qctrigger2=calculateqctrigger_threedividends_part2_fixedr(S0, Shistory, ctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time);
                                qctrigger3=calculateqctrigger_threedividends_part3_fixedr(S0, Shistory, ctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time,D2,D3_time);
                                qctrigger4=calculateqctrigger_threedividends_part4_fixedr(S0, Shistory, ctrigger,T,r,sigma,no_call_time,n1,m1,ncX,D1_time,D1,D2_time,D2,D3_time,D3);
                            else
                                qctrigger1=S0*5;
                                qctrigger2=S0*5;
                                qctrigger3=S0*5;
                                qctrigger4=S0*5;
                            end
                            % equivalent put trigger
                            if n2>0 && (T-D1_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger1=calculateqptrigger_threedividends_part1_fixedr(S0, Shistory, ptrigger,T,r,sigma,no_put_time,n2,m2,D1_time);
                            else
                                qptrigger1=0;
                            end
                            
                            if n2>0 && (T-D2_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger2=calculateqptrigger_threedividends_part2_fixedr(S0, Shistory, ptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time);
                            else
                                qptrigger2=0;
                            end
                            
                            if n2>0 && (T-D3_time<no_put_time(1))
                                ptrigger=X*npX;
                                qptrigger3=calculateqptrigger_threedividends_part3_fixedr(S0, Shistory, ptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time,D2,D3_time);
                            else
                                qptrigger3=0;
                            end
                            
                            if n2>0
                                ptrigger=X*npX;
                                qptrigger4=calculateqptrigger_threedividends_part4_fixedr(S0, Shistory, ptrigger,T,r,sigma,no_put_time,n2,m2,npX,D1_time,D1,D2_time,D2,D3_time,D3);
                            else
                                qptrigger4=0;
                            end
                        end
                        %%
                        %¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªcomputation from parameter¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
                        Smax=max(max(max(qctrigger1,qctrigger2),qctrigger3),qctrigger4);%upbound of stock price
                        %Smax=qctrigger1;
                        Ns=500;
                        ds=Smax/Ns;
                        Nt=500;
                        dt=T/Nt;%time direction
                        Large=1/min(ds^2,dt^1.5);
                        rhopenltycall=Large;
                        rhopenltyput=Large;
                        tol=1/Large;
                        
                        %%
                        %¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªboundary conditions¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
                        %i=(1:Ns+1)';%stock price, [1,2,...,Ns+1],u(x,y)=u(T-xdt,(y-1)dh)
                        S=(0:Ns)'*ds;%stock price, [0,h,2h,...,Ns*h]
                        Bc_T=F+coupon_rate(1)*F;%Bc in the final time,using dirty price
                        if n2>0
                            Bp_T=F+coupon_rate(1)*F;
                        else
                            Bp_T=0;
                        end
                        k_T=F/(X-D1-D2-D3);%dividend influence
                        u(1:Ns+1,:)=max(k_T*S,F+coupon_rate(1)*F+bonus_rate*F);%from t=0 update u to t=T£»this is initial value of u when t=0,t=T-t,backward
                        B(1:Ns+1,:)=F+coupon_rate(1)*F+bonus_rate*F;%from t=0 update B to t=T£»this is initial value of B when t=0,t=T-t,backward
                    
                        
                        %¡ª¡ª¡ª¡ª¡ª¡ªsolve PDE by difference method, combined with penalty method¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
                        
                        flag=0;%if not convergent, break will stop the whole loop
                        for n=1:Nt %t=T-n*dt,n=1 means t=T-dt,n=N tmeans t=0 %want to u^n¡úu^(n+1)
                            count=0;
                            count1=0;
                            count2=0;
                            count3=0;
                            count4=0;
                            count5=0;
                            count6=0;
                            if flag==1
                                break;
                            end
                            
                            t=n*dt;%time to maturaty
                            
                            S=(0:Ns)'*ds+D3*(t>T-D3_time)+D2*(t>T-D2_time)+D1*(t>T-D1_time);%dividend influence to stock price
                            
                            
                            qctrigger_new=qctrigger1*(t>T-D1_time)+qctrigger2*(t>T-D2_time)*(t<=T-D1_time)+qctrigger3*(t>T-D3_time)*(t<=T-D2_time)+qctrigger4*(t<=T-D3_time);
                            qptrigger_new=qptrigger1*(t>T-D1_time)+qptrigger2*(t>T-D2_time)*(t<=T-D1_time)+qptrigger3*(t>T-D3_time)*(t<=T-D2_time)+qptrigger4*(t<=T-D3_time);
                            
                            K=0;  %coupon payment

                            alpha=((sigma^2*S(2:Ns,:).^2)/(2*ds^2)-(r+p*eta-q)*S(2:Ns,:)./(2*ds))*dt;%Ns-1, namely alpha_1~alpha_(Ns-1) means ds~(Ns-1)*ds
                            beta=((sigma^2*S(2:Ns,:).^2)/(2*ds^2)+(r+p*eta-q)*S(2:Ns,:)./(2*ds))*dt;
                            alpha_MuMB=[alpha;0;0];%take Mu and MB down diagonals as alpha_MuMB
                            beta_MuMB=[0;0;beta];%take Mu and MB up diagonals as beta_MuMB
                            if t<no_convert_time(1)
                                gamma_Mu=[-(r+p)*dt;-(alpha+beta+(r+p)*dt);0];%take Mu diagonals as gamma_Mu
                            else
                                gamma_Mu=[-(r+p)*dt;-(alpha+beta+(r+p)*dt);-(p*(1-eta)+q)*dt];%take Mu diagonals as gamma_Mu
                                %gamma_Mu=[-(r+p)*dt;-(alpha+beta+(r+p)*dt);0];
                            end
                            Mu=spdiags([alpha_MuMB,gamma_Mu,beta_MuMB],(-1:1),Ns+1,Ns+1);%this is Mu
                            gamma_MB=[-(r+p*(1-R))*dt;-(alpha+beta+(r+p*(1-R))*dt);0];%take MB diagonals as gamma_MB
                            MB=spdiags([alpha_MuMB,gamma_MB,beta_MuMB],(-1:1),Ns+1,Ns+1);%This is MB
                            
                            if ceil(t+dt+1e-10)<=length(coupon_rate)
                                AccI=(ceil(t+dt-1e-10)-t)*coupon_rate(ceil(t+dt-1e-10))*F;
                            end
                            
                            k_n=0*(t>=no_convert_time(1))+F/X*(t<no_convert_time(1))*(t>T-D1_time)+F/(X-D1)*(t<=T-D1_time)*(t>T-D2_time)+F/(X-D1-D2)*(t<=T-D2_time)*(t>T-D3_time)+F/(X-D1-D2-D3)*(t<=T-D3_time);
                            
                            if (n1>0) && (t<no_call_time(1))
                                Bc_n=F+AccI; %dirty price
                            else
                                Bc_n=100000;
                            end
                            if (n2>0) && (t<no_put_time(1)-n2/252)
                                Bp_n=F+AccI;
                            else
                                Bp_n=0;
                            end
                            u_n=u;%now value of u is n-1 step, will converte to n step following
                            B=(eye(Ns+1)-theta*MB)\((eye(Ns+1)+(1-theta)*MB)*B+K);%similar to eqation (3.45)£¬without constratint, B's eqation,
                            % solution B_n is n step iteration initial value
                            B(B>Bc_n)=Bc_n;%¡úexplict B<=Bc_n
                            Muu=p*dt*max(k_n*(1-eta)*S,R*B);%similar to (3.44) the right term
                            u=(speye(Ns+1)-theta*Mu)\((speye(Ns+1)+(1-theta)*Mu)*u_n+[Muu(1:Ns,:);0]+K);
                            %similar to eqation (3.44)£¬without constratint, u's eqation, solution u_n is n step iteration initial value
                            up=zeros(length(S),1);
                            low=zeros(length(S),1);
                            for j=1:length(S)
                                if S(j)>=qctrigger_new
                                    if t<no_convert_time(1)
                                        u(j)=k_n*qctrigger_new;
                                        up(j)=k_n*qctrigger_new;
                                        low(j)=k_n*qctrigger_new;
                                    else
                                        u(j)=F/X*qctrigger_new;
                                        up(j)=F/X*qctrigger_new;
                                        low(j)=F/X*qctrigger_new;
                                    end
                                elseif S(j)<qptrigger_new
                                    u(j)=max(u(j),max(k_n*S(j),Bp_n));
                                    up(j)=100000;
                                    low(j)=max(k_n*S(j),Bp_n);
                                else
                                    u(j)=max(u(j),k_n*S(j));
                                    %u(j)=max(u(j),Bp_n);
                                    up(j)=100000;
                                    low(j)=k_n*S(j);
                                    %low(j)=Bp_n;
                                end
                            end
                            
                            %% old method
                            P1=rhopenltyput*spdiags(u<low,0,Ns+1,Ns+1);
                            P2=rhopenltycall*spdiags(u>up,0,Ns+1,Ns+1);
                            while true
                                count=count+1;
                                u_old=u;%u_old is u^0_n,i.e.,n step initial value
                                P1_old=P1;
                                P2_old=P2;%fix n,every interation of k, u,P1,P2 will change value, so add 'old'
                                u=(speye(Ns+1)-theta*Mu+P1_old+P2_old)\((speye(Ns+1)+(1-theta)*Mu)*u_n+[Muu(1:Ns,:);0]+K+P1_old*low+P2_old*up);
                                %similar to (4.28)£¬solution is u_n^(k+1)
                                P1=rhopenltyput*spdiags(u<low,0,Ns+1,Ns+1);
                                P2=rhopenltycall*spdiags(u>up,0,Ns+1,Ns+1);
                                %         if (isequal(P1, P1_old) && isequal(P2, P2_old))
                                %             count1=count1+1;
                                %             break;
                                %         else
                                if max(abs(u-u_old)./max(1,abs(u)))<tol
                                    count2=count2+1;
                                    break;
                                elseif count>1e4
                                    count3=count3+1;
                                    disp(['Iteration probably does not converge, current n is ' num2str(n)]);
                                    flag=1;
                                    break;
                                else
                                    count4=count4+1;
                                end
                            end
                            err=max(abs(u-u_old)./max(1,abs(u)));
                            disp(['current n=' num2str(n) ', count=' num2str(count) ', error=' num2str(err)])
                            
                            %% new method
                            %     uu=u;
                            %     P1_old=rhopenltyput*speye(Ns+1);
                            %     P2_old=rhopenltycall*speye(Ns+1);
                            %     while true
                            %         count=count+1;
                            %         u_old=u;%u_old is u^0_n,i.e.,n step initial value
                            %         P1=rhopenltyput*spdiags(u+uu<2*low,0,Ns+1,Ns+1);
                            %         P2=rhopenltycall*spdiags(u+uu>2*up,0,Ns+1,Ns+1);
                            %         u=(speye(Ns+1)-theta*Mu+P1+P2)\((speye(Ns+1)+(1-theta)*Mu)*u_n+[Muu(1:Ns,:);0]-P1*(uu-2*low)-P2*(uu-2*up));
                            %         %similar to (4.28)£¬solution is u_n^(k+1)
                            %         P1=rhopenltyput*spdiags(u+uu<=2*low,0,Ns+1,Ns+1);
                            %         P2=rhopenltycall*spdiags(u+uu>=2*up,0,Ns+1,Ns+1);
                            %         if (isequal(P1, P1_old) && isequal(P2, P2_old))
                            %             count1=count1+1;
                            %             break;
                            %         elseif max(abs(u-u_old)./max(1,abs(u)))<tol
                            %             count2=count2+1;
                            %                 break;
                            %         elseif count>1e4
                            %             count4=count4+1;
                            %             disp(['Probably does not converge, current n is ' num2str(n)]);
                            %             flag=1;
                            %             break;
                            %         else
                            %             count3=count3+1;
                            %         end
                            %         P1_old=P1;
                            %         P2_old=P2;
                            %     end
                            %     err=max(abs(u-u_old)./max(1,abs(u)));
                            %     disp(['current n=' num2str(n) ', count=' num2str(count) ', error=' num2str(err)])
                            
                            B=min(B,u);
                            if (floor(t+dt+1e-10)-floor(t+1e-10))==1 && floor(t+dt+1e-10)+1<=length(coupon_rate)  %coupon rate payment
                                count5=count5+1;
                                for j=1:length(S)
                                    if S(j)<qctrigger_new
                                        u(j)=u(j)+coupon_rate(floor(t+dt+1e-10)+1)*F;
                                        %B(j)=B(j)+coupon_rate(ceil(t+dt-1e-10))*F;
                                    end
                                end
                                %u=u+coupon_rate(ceil(t+dt-1e-10))*F;
                                B=B+coupon_rate(floor(t+dt+1e-10)+1)*F;
                            end
                        end
                        Uindex=Ns*(S0-S(1))/(S(end)-S(1))+1;
                        U=u(ceil(Uindex-1e-10))*(Uindex-floor(Uindex+1e-10))+u(floor(Uindex+1e-10))*(ceil(Uindex-1e-10)-Uindex)
                        disp(['Total elapsed time : ' num2str(toc) ', U = ' num2str(U)])
                        fprintf(fid, 'Total elapsed time: %f \n', toc);
                        fprintf(fid, 'Bond Value: %f \n', U);
                        %%
                    end % case_method - 2
                end % case_S0 - 6 
            end % case sigma - 7
        end % case Dtime - 3
    end % case D - 6
end % three different cases - 3
    
fclose(fid);



% usingtime=toc;
% Out=['The convertible value is ', num2str(U)];
% disp(Out);
% disp(['Total elapsed time : ' num2str(usingtime) ', Nt = ' num2str(Nt) ', Ns = ' num2str(Ns)]);
% end