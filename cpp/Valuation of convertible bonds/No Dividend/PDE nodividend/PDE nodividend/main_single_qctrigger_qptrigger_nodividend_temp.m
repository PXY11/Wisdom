% function U = convertible_code(Nt,Ns,F,conversion_price,Smax,S0,risk_free_rate,rate_T,sigma,T,final_coupon_rate,q,p,R,eta,no_call_time,...
% no_put_time,no_convert_time,Bc_star,Bp_star,theta,coupon_time_rate)
%% ¡ª¡ª¡ª¡ª¡ªSolving pde-Penalty mehtod-CN difference-Newton iteration¡ª¡ª¡ª¡ª¡ª¡ª
clc;
clear;
tic;
format long;
%%
% ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªparameter specification¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
% coupon_rate=[2 1.5 1 0.8 0.6]/100; bonus_rate=4/100;
% S0=6.95; X=7.26; T=2.4904; n1=15; m1=30; ncX=1.3; n2=30; m2=30; npX=0.7; sigma=0.3176; %sigma=0.54148; r=0.047;
% r=0.0345; oas=0.0;    %110030
% no_convert_time=[4.5 5];
% no_call_time=[4.5 5];
% no_put_time=[2 5];
% %T=4.7;
% Shistory=zeros(1,29);
% %S0=4.5;

% coupon_rate=[1.6 1.5 1.5 1 0.5 0.2]/100; bonus_rate=5.4/100; 
% S0=20.71; X=42.80; T=3.9562; n1=20; m1=30; ncX=1.3; sigma=0.2141; 
% r=0.035; oas=0.0;    %110031
% no_convert_time=[5.5 6];
% no_call_time=[5.5 6];
% no_put_time=[2 6];
% 
% coupon_rate=[2 1.6 1.5 1 0.5 0.2]/100;  bonus_rate=4/100;
% S0=8.17; X=7.46; T=4.5205; n1=15; m1=30; nX=1.3; sigma=0.2743;
% r=0.035; oas=0.0;     %110032
% 
% coupon_rate=[2 1.7 1.4 0.9 0.5 0.3]/100;  bonus_rate=6/100; 
% S0=9.06; X=8.93; T=4.5233; n1=15; m1=30; nX=1.3; sigma=0.2216;
% r=0.0345; oas=0.0;     %110033
% 
% coupon_rate=[2 1.6 0.8 0.6 0.4 0.2]/100;  bonus_rate=6/100;
% S0=21.16; X=18.65; T=4.5507; n1=20; m1=30; nX=1.3; sigma=0.2397;
% r=0.035; oas=0.0;     %110034
% 
% coupon_rate=[1.6 1.5 1.5 1 0.5 0.2]/100;  bonus_rate=5/100;
% S0=7.57; X=10.65; T=3.6; n1=15; m1=30; nX=1.3; sigma=0.4;
% r=0.035; oas=0.0;    %113008
% 
% coupon_rate=[1.6 1.5 1.5 1 0.5 0.2]/100;  bonus_rate=4.4/100;
% S0=26.14; X=21.53; T=4.5699; n1=15; m1=30; nX=1.3; sigma=0.2458;
% r=0.035; oas=0.0;    %113009
% 
% coupon_rate=[2 1.8 1.5 1 0.5 0.3]/100;  bonus_rate=6/100;
% S0=6.89; X=9.3; T=4.7205; n1=15; m1=30; nX=1.3; sigma=0.2164;
% r=0.035; oas=0.0;    %113010
% 
% coupon_rate=[2 1.8 1.5 1 0.5 0.2]/100;  bonus_rate=3/100;
% S0=4.04; X=4.36; T=5.7178; n1=15; m1=30; nX=1.3; sigma=0.1429;
% r=0.035; oas=0.0;    %113011
% 
% coupon_rate=[1.8 1.5 1.3 1 0.5 0.3]/100;  bonus_rate=3.2/100;
% S0=14.2; X=16.78; T=5.7370; n1=15; m1=30; nX=1.3; sigma=0.2595;
% r=0.035; oas=0.0;    %113012
% 
% coupon_rate=[1 1 1 1 1]/100;  bonus_rate=7/100;
% S0=17.1; X=17.8; T=3.8055; n1=10; m1=20; nX=1.3; sigma=0.2139;
% r=0.035; oas=0.0;   %120001
% 
% coupon_rate=[2 1.8 1.5 1 0.7 0.5]/100;  bonus_rate=6/100;
% S0=7.81; X=9.93; T=4.474; n1=15; m1=30; nX=1.3; sigma=0.2695;
% r=0.035; oas=0.0;   %123001
% 
% coupon_rate=[2 1.8 1.5 1 0.7 0.5]/100;   bonus_rate=8/100;
% S0=3.7; X=5.26; T=4.9425; n1=15; m1=30; nX=1.3; sigma=0.1943;
% r=0.035; oas=0.0;   %127003
% 
% coupon_rate=[2 1.8 1.5 1 0.7 0.5]/100;   bonus_rate=8/100;
% S0=8.09; X=8; T=5.9260; n1=20; m1=30; nX=1.3; sigma=0.2655;
% r=0.036; oas=0.0;    %127004
% 
% coupon_rate=[1.6 1.6 1.6 1 0.7 0.5]/100;   bonus_rate=6.4/100;
% S0=19.14; X=13.04; T=3.4575; n1=15; m1=30; nX=1.3; sigma=0.2659;
% r=0.035; oas=0.0;    %128009
% 
% no_convert_time=[5.5 6];
% no_call_time=[5.5 6];
% no_put_time=[2 6];

%Shistory=[18.10 17.73 17.79 ... 18.84];
% 
% coupon_rate=[1.6 1.6 1.6 1 0.7 0.5]/100;   bonus_rate=6.4/100;
% S0=8.94; X=9.32; T=4.5671; n1=15; m1=30; nX=1.3; sigma=0.3314;
% r=0.035; oas=0.0;    %128010
% 
coupon_rate=[1.8 1.5 1.5 1 0.7 0.5]/100;   bonus_rate=4.2/100;
S0=7.68; X=5.73; T=4.674; n1=15; m1=30; ncX=1.3; n2=30; m2=30; npX=0.7; sigma=0.4;
%T=5.674;
r0=0.035; oas=0.0;    %128011
r=r0+oas;
no_convert_time=[5.5 6];
no_call_time=[5.5 6];
no_put_time=[2 6];
T=5.9;
S0=15;
Shistory=[6.25 6.13 6.23 6.03 6.19 6.32 6.39 6.33 6.20 6.04 6.13 6.12 6.11 6.25 6.19 6.09 6.01 6.13 6.13 6.24 6.33 6.80 7.48 7.70 8.30 7.75 7.77 7.81 7.61];
% 
% coupon_rate=[1.6 1.3 1.3 1 0.7 0.5]/100;  bonus_rate=1.4/100;
% S0=4.95; X=7.79; T=4.8137; n1=15; m1=30; nX=1.3; sigma=0.2068; 
% r=0.035; oas=0.0;   %128012
% 
% coupon_rate=[2 1.8 1.5 1 0.6 0.4]/100;  bonus_rate=6/100;
% S0=7.29; X=10.06; T=5.0822; n1=15; m1=30; nX=1.3; sigma=0.3171;
% r=0.036; oas=0.0;    %128013
% 
% coupon_rate=[2 1.8 1.5 1 0.7 0.5]/100;  bonus_rate=6/100;
% S0=19.04; X=30.77; T=5.8027; n1=15; m1=30; nX=1.3; sigma=0.4635;
% r=0.035; oas=0.0;     %128014
% 
% coupon_rate=[1.8 1.5 1.3 1 0.5 0.3]/100;   bonus_rate=6.2/100;
% S0=12.59; X=12.97; T=5.9425; n1=15; m1=30; nX=1.3; sigma=0.329;
% r=0.037; oas=0.0;    %128015
% 
% coupon_rate=[1.5 1.5 1.5]/100;   bonus_rate=1.5/100;
% S0=50.58; X=42.79; T=0.4493; n1=0; m1=0; nX=1.3; sigma=0.2365; %wind sigma=0.241331
% r=0.034; oas=0.0;   %132001
% 
% coupon_rate=[1 1 1 1 1]/100;   bonus_rate=0;
% S0=41.3; X=56.02; T=2.9452; n1=10; m1=20; nX=1.35; sigma=0.1994; 
% r=0.035; oas=0.0;   %132002
% 
% coupon_rate=[1 1 1]/100;   bonus_rate=0;
% S0=11.77; X=16.75; T=1.326; n1=10; m1=30; nX=1.2; sigma=0.2294; 
% r=0.034; oas=0.0;   %132003 
% 
% coupon_rate=[1 1 1 1 1 1]/100;   bonus_rate=3/100;
% S0=3.81; X=6.76; T=4.3562; n1=15; m1=30; nX=1.3; sigma=0.2151;
% r=0.035; oas=0.0;   %132004
% 
% coupon_rate=[1.7 1.7 1.7 1.7 1.7]/100;   bonus_rate=7.5/100;
% S0=33.22; X=38.34; T=3.4466; n1=0; m1=0; nX=1.3; sigma=0.2065;
% r=0.035; oas=0.0;   %132005
% 
% coupon_rate=[1 1 1 1 1]/100;  bonus_rate=11/100;
% S0=13.96; X=16.19; T=3.9863; n1=15; m1=30; nX=1.2; sigma=0.4859;
% r=0.035; oas=0.0;   %132006
% 
% coupon_rate=[1 1 1 1 1]/100;   bonus_rate=5/100;
% S0=9.79; X=16; T=4.3425; n1=10; m1=20; nX=1.3; sigma=0.1242;
% r=0.035; oas=0.0;   %132007
% 
% coupon_rate=[1.7 1.7 1.7 1.7 1.7]/100;  bonus_rate=6/100;
% S0=6.27; X=10; T=4.8219; n1=15; m1=30; nX=1.3; sigma=0.2403;
% r=0.035; oas=0.0;   %132008 

%r=r+oas;

%%
F=100;%face value
q=0;% dividend
p=0;% default probobility
R=1;% recovery rate
eta=0;%when default the stock price becomes (1-eta)*S
theta=1;%implicitness parameter£¬0.5ÎªCN method,1Îªfully implicit or backward euler method

%%
% r=0.03;
% T=0.5;% time to maturaty
% T=1.5;
% %T=2.5;
% X=42.79;% convert price
% S0=52.75;%now stock price
% %S0=30;
% sigma=0.3034;% volatility
% %bonus_rate=0;
% bonus_rate=0.015;
% %no_call_time=[0 3];%no allow call time
% no_call_time=[2 3];
% %no_put_time=[0 3];%no allow put time
% no_put_time=[1 3];
% no_convert_time=[2 3];%no allow convert time
% coupon_rate=[0.015 0.015 0.015 0 0 0 0 0 0 0 0 0];%coupon rate
% %coupon_rate=[0 0 0 0 0 0 0 0 0];
% n1=15;
% %n2=30;
% n2=0;

%% equivalent call trigger
% if n1>0
%     ctrigger=X*ncX;     %call trigger in announcement
%     qctrigger_fun=@(qctrigger) calculateqctrigger_fixedr_test1(S0,ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1);
%     qctrigger0=ctrigger*1.08;     %initial value of equivalent call trigger
%     %options=optimset('Display','iter','TolX',0.0001);
%     options=optimset('Display','iter','TolFun',0.0001);
%     qctrigger=fzero(qctrigger_fun,qctrigger0,options);
% else
%     qctrigger=S0*5;
% end

% if n1>0
%     ctrigger=X*ncX;     %call trigger in announcement
%     qctrigger_fun=@(qctrigger) calculateqctrigger_fixedr_test2(S0,Shistory,ctrigger,qctrigger,T,r,sigma,no_call_time,n1,m1);
%     qctrigger0=ctrigger*1.08;     %initial value of equivalent call trigger
%     %options=optimset('Display','iter','TolX',0.0001);
%     options=optimset('Display','iter','TolFun',0.0001);
%     qctrigger=fzero(qctrigger_fun,qctrigger0,options);
% else
%     qctrigger=S0*5;
% end   %probability method

% if n1>0
%     ctrigger=X*ncX;
%     qctrigger_fun=@(qctrigger) calculateqctrigger_fixedr_test3(S0, ctrigger, qctrigger,T,r,F,X,coupon_rate,sigma,no_call_time,n1,m1);
%     qctrigger0=ctrigger*1.08;     %initial value of equivalent call trigger
%     options=optimset('Display','iter','TolFun',0.00001);
%     qctrigger=fzero(qctrigger_fun,qctrigger0,options);
% else
%     qctrigger=S0*5;
% end

% if n1>0
%     ctrigger=X*ncX;
%     qctrigger_fun=@(qctrigger) calculateqctrigger_fixedr_test4(S0, Shistory, ctrigger, qctrigger,T,r,F,X,coupon_rate,sigma,no_call_time,n1,m1);
%     qctrigger0=ctrigger*1.05;     %initial value of equivalent call trigger
%     options=optimset('Display','iter','TolFun',0.00001);
%     qctrigger=fzero(qctrigger_fun,qctrigger0,options);
% else
%     qctrigger=S0*5;
% end

if n1>0
    ctrigger=X*ncX;
    qctrigger=calculateqctrigger_fixedr_test7(S0, Shistory, ctrigger,T,r0,sigma,no_call_time,n1,m1);
else
    qctrigger=S0*5;
end   %mean method

%% equivalent put trigger
% if n2>0
%     ptrigger=X*npX;     %call trigger in announcement
%     qptrigger_fun=@(qptrigger) calculateqptrigger_fixedr(S0,ptrigger,qptrigger,T,r,sigma,n2,m2);
%     qptrigger0=ptrigger*0.95;     %initial value of equivalent call trigger
%     %options=optimset('Display','iter','TolX',0.0001);
%     options=optimset('Display','iter','TolFun',0.0001);
%     qptrigger=fzero(qptrigger_fun,qptrigger0,options);
% else
%     qptrigger=0;
% end

% if n2>0
%     ptrigger=X*npX;
%     qptrigger_fun=@(qptrigger) calculateqptrigger_fixedr_test2(S0, Shistory, ptrigger, qptrigger,T,r,sigma,no_put_time,n2,m2);
%     qptrigger0=ptrigger*0.95;     %initial value of equivalent call trigger
%     options=optimset('Display','iter','TolFun',0.0001);
%     qptrigger=fzero(qptrigger_fun,qptrigger0,options);
% else
%     qptrigger=0;
% end   %probability method

if n2>0
    ptrigger=X*npX;
    qptrigger=calculateqptrigger_fixedr_test3(S0, Shistory, ptrigger,T,r0,sigma,no_put_time,n2,m2);
else
    qptrigger=0;
end     %mean method

% qctrigger=X*1.3;
% qptrigger=X*0.7;

%qctrigger=10.07;

%%
D=0;%fixed dividend
dtime=0.3301;%dividend at time t=X.33

%%
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªcomputation from parameter¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
% ds=0.0556;
% Smax=ceil(qctrigger/ds)*ds;%upbound of stock price
% Ns=round(Smax/ds);%price direction

Smax=qctrigger;
Ns=1000;
ds=Smax/Ns;

Nt=500;
dt=T/Nt;%time direction
Large=1/min(ds^2,dt^1.5)
Large=1e6;
rhopenltycall=Large;
rhopenltyput=Large;
tol=1/Large

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
k_T=F/(X-D*ceil(T-dtime));%dividend influence
u(1:Ns+1,:)=max(k_T*S,F+coupon_rate(1)*F+bonus_rate*F);%from t=0 update u to t=T£»this is initial value of u when t=0,t=T-t,backward
%u=max(u,Bp_T); 
%u=min(u,max(Bc_T,k_T*qctrigger));%¡úexplict determine value of u by considering constraint

B(1:Ns+1,:)=F+coupon_rate(1)*F+bonus_rate*F;%from t=0 update B to t=T£»this is initial value of B when t=0,t=T-t,backward
%B(1:Ns+1,:)=F+coupon_time_rate(1)*F; %I am unsure which one is correct.(by Li Xin)

% umax(1:Nt)=k*Smax;%s=Smax,u value, when u_ss=0
% Bmax(1:Nt)=F+final_coupon_rate*F;%s=Smax,B value,B_ss=0;¡ª¡ª%this 
% two condition is including in the matrix following, so they are useless

%%
%¡ª¡ª¡ª¡ª¡ª¡ªsolve pde by difference method, combined penalty method¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
flag=0;%if not convergenc, break will stop the whole loop
%%
for n=1:Nt%t=T-n*dt,n=1 means t=T-dt,n=N tmeans t=0 %want to u^n¡úu^(n+1)
    count=0;%total iteration time
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
    
    S=(0:Ns)'*ds+D*ceil(t-dtime);%dividend influence to stock price
    
%     if (ceil(t+dt)-ceil(t))==1%coupon rate payment
%         count5=count5+1;
%         u=u+coupon_rate(ceil(t+dt+1e-6))*F;
%         B=B+coupon_rate(ceil(t+dt+1e-6))*F;
%     end

    K=0;  %coupon payment
%     if (floor(t+dt)-floor(t))==1
%         K=coupon_rate(floor(t+dt)+1)*F;
%     end
        
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
       
    if t<no_convert_time(1)
        k_n=F/(X-D*(ceil(T-dtime)-ceil(t-dtime)));%when T=X.dtime, stock pay dividends
    else
        k_n=0;
    end
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
%     u=max(u,max(Bp_n,k_n*S));
%     u=min(u,max(Bc_n,k_n*qctrigger));   %¡úexplict constraint u, initial value of u_n
%     low=max(Bp_n,k_n*S);
%     up=max(Bc_n,k_n*qctrigger)*ones(length(S),1);
%     u=max(u,low);
%     u=min(u,up);
    
    up=zeros(length(S),1);
    low=zeros(length(S),1);
    for j=1:length(S)
        if S(j)>=qctrigger
            u(j)=max(u(j),k_n*S(j));
            u(j)=min(u(j),max(k_n*S(j),Bc_n));
            up(j)=max(k_n*S(j),Bc_n);
            low(j)=k_n*S(j);
        elseif S(j)<qptrigger
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
        P1=rhopenltyput*spdiags(u<low,0,Ns+1,Ns+1);
        P2=rhopenltycall*spdiags(u>up,0,Ns+1,Ns+1);
        if (isequal(P1, P1_old) && isequal(P2, P2_old))
            count1=count1+1;
            break;
        elseif max(abs(u-u_old)./max(1,abs(u)))<tol
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
%     P1_old=rhopenltyput*spdiags(u<max(Bp_n,k_n*S),0,Ns+1,Ns+1);
%     P2_old=rhopenltycall*spdiags(u>k_n*qctrigger*ones(Ns+1,1),0,Ns+1,Ns+1);
%     while true
%         count=count+1;
%         u_old=u;%u_old is u^0_n,i.e.,n step initial value
%         P1=rhopenltyput*spdiags(u+uu<2*max(Bp_n,k_n*S),0,Ns+1,Ns+1);
%         P2=rhopenltycall*spdiags(u+uu>2*k_n*qctrigger,0,Ns+1,Ns+1);
%         %u=(speye(Ns+1)-theta*Mu+P1-P2)\((speye(Ns+1)+(1-theta)*Mu)*u_n+[Muu(1:Ns,:);0]+P1*max(Bp_n,k_n*S)-P2*k_n*qctrigger*ones(Ns+1,1));
%         u=(speye(Ns+1)-theta*Mu+P1-P2)\((speye(Ns+1)+(1-theta)*Mu)*u_n+[Muu(1:Ns,:);0]-P1*(uu-2*max(Bp_n,k_n*S))+P2*(uu-2*k_n*qctrigger*ones(Ns+1,1)));
%         %similar to (4.28)£¬solution is u_n^(k+1)
%         P1=rhopenltyput*spdiags(u+uu<2*max(Bp_n,k_n*S),0,Ns+1,Ns+1);
%         P2=rhopenltycall*spdiags(u+uu>2*k_n*qctrigger*ones(Ns+1,1),0,Ns+1,Ns+1);
%         if (isequal(P1, P1_old) && isequal(P2, P2_old))
%             count1=count1+1;
%             break;
%         elseif max(abs(u-u_old)./max(1,abs(u)))<=1e-3
%             count2=count2+1;
%                 break;
%             elseif count>Nt*40
%                     count4=count4+1;
%                     disp(['Probably does not converge, current n is ' num2str(n)]);
%                     err=max(abs(u-u_old)./max(1,abs(u)))
%                     flag=1;
%                     break;
%         else
%             count3=count3+1;
%         end        
%         P1_old=P1;
%         P2_old=P2;
%     end

    B=min(B,u);
    if (floor(t+dt+1e-10)-floor(t+1e-10))==1 && floor(t+dt+1e-10)+1<=length(coupon_rate)  %coupon rate payment
        count5=count5+1;
        for j=1:length(S)
            if S(j)<qctrigger
                u(j)=u(j)+coupon_rate(floor(t+dt+1e-10)+1)*F;
                %B(j)=B(j)+coupon_rate(ceil(t+dt-1e-10))*F;
            end
        end
        %u=u+coupon_rate(ceil(t+dt-1e-10))*F;
        B=B+coupon_rate(floor(t+dt+1e-10)+1)*F;
    end
end

%%
Uindex=Ns*(S0-S(1))/(S(end)-S(1))+1;
U=u(ceil(Uindex-1e-10))*(Uindex-floor(Uindex+1e-10))+u(floor(Uindex+1e-10))*(ceil(Uindex-1e-10)-Uindex)
disp(['Total elapsed time : ' num2str(toc) ', U = ' num2str(U)]);

% usingtime=toc;
% Out=['The convertible value is ', num2str(U)];
% disp(Out);
% disp(['Total elapsed time : ' num2str(usingtime) ', Nt = ' num2str(Nt) ', Ns = ' num2str(Ns)]);
% end