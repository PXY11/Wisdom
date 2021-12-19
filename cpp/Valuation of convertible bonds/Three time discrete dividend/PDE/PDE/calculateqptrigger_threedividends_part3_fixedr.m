function qptrigger = calculateqptrigger_threedividends_part3_fixedr(S0, Shistory, ptrigger,T,r,sigma,no_put_time,n,m,npX,D1_time,D1,D2_time,D2,D3_time)
%% parameter input
% clc;
% clear;
% S0=18.21;
% Bc=43.05*1.3;
% Bc_star=62.2146;
% T=4.08;
% u=-0.01;
% sigma=0.42;
% n=15;
% m=30;

N=5000;
N=20000;%numbers of MC
N=50000;
% d_t=1;
% M=floor(T*252);
d_t=1/252;
M=round(T/d_t);
%% calculate
total_left_Smark=0;
%total_left_Smark2=0;
countput=0;
%countcall2=0;
filter_a=1;
filter_b=ones(1,m);
rng('shuffle');
%sam1=random('norm',0,1,N,M);
sam1=randn(N,M);
sam2=-sam1;
sam=[sam1;sam2];
S=ones(1,M)*10000;
for i=1:2*N
    S(1)=S0*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,1));
    for j=2:D1_time/d_t-1
        S(j)=S(j-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,j));
    end
    Stmp1=S(D1_time/d_t-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,D1_time/d_t));
    S(D1_time/d_t)=Stmp1-D1;
    for j=D1_time/d_t+1:D2_time/d_t-1
        S(j)=S(j-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,j));
    end
    Stmp2=S(D2_time/d_t-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,D2_time/d_t));
    S(D2_time/d_t)=Stmp2-D2;
    for j=D2_time/d_t+1:D3_time/d_t-1
        S(j)=S(j-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,j));
    end
    %S=10*ones(1,M);
    SS=[Shistory,S0,S];
    
    ptrigger_sequence=zeros(1,M+m);
    ptrigger_sequence(1:D1_time/d_t-1+m)=ptrigger;
    ptrigger_sequence(D1_time/d_t:D2_time/d_t-1+m)=ptrigger-npX*D1;
    ptrigger_sequence(D2_time/d_t:M+m)=ptrigger-npX*(D1+D2);
    
    if T+(m-1)*d_t<=no_put_time(1)-d_t
        judge=(SS<ptrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m
            left_triggertime=(left_j-m)*d_t;
            if left_triggertime>=D2_time  %分红后
                left_Smark=SS(left_j);
                countput=countput+1;
                total_left_Smark=total_left_Smark+left_Smark;
%             else  %分红前
%                 left_Smark2=SS(left_j);
%                 countcall2=countcall2+1;
%                 total_left_Smark2=total_left_Smark2+left_Smark2;
            end
        elseif ~isempty(left_j)
            disp('1 the convertible bond was already put')
        end
    elseif T+(m-1)*d_t>=no_put_time(1) && T<=no_put_time(1)-d_t
        left_delta=T+(m-1)*d_t-no_put_time(1);
        SS=SS(floor(left_delta/d_t+1e-10)+2:end);
        ptrigger_sequence=ptrigger_sequence(floor(left_delta/d_t+1e-10)+2:end);
        judge=(SS<ptrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m-1-floor(left_delta/d_t+1e-10)
            left_triggertime=(left_j-m+1+floor(left_delta/d_t+1e-10))*d_t;
            if left_triggertime>=D2_time
                left_Smark=SS(left_j);
                countput=countput+1;
                total_left_Smark=total_left_Smark+left_Smark;
%             else
%                 left_Smark2=SS(left_j);
%                 countcall2=countcall2+1;
%                 total_left_Smark2=total_left_Smark2+left_Smark2;
            end
        elseif ~isempty(left_j)
            disp('2 the convertible bond was already put')
        end      
    else
        %T>no_call_time(1)-d_t
        left_delta=T-no_put_time(1)+d_t;
        SS=S(ceil(left_delta/d_t-1e-10):end);
        ptrigger_sequence=ptrigger_sequence(ceil(left_delta/d_t-1e-10)+m:end);
        judge=(SS<ptrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if ~isempty(left_j)
            left_triggertime=(left_j-1+ceil(left_delta/d_t-1e-10))*d_t;
            if left_triggertime>=D2_time
                left_Smark=SS(left_j);
                countput=countput+1;
                total_left_Smark=total_left_Smark+left_Smark;
%             else
%                 left_Smark2=SS(left_j);
%                 countcall2=countcall2+1;
%                 total_left_Smark2=total_left_Smark2+left_Smark2;
            end
        end
    end
end
if countput>0
    qptrigger=total_left_Smark/countput;
else
    qptrigger=0;
end
%qctrigger2=(total_left_Smark1+total_left_Smark2)/(countcall1+countcall2);
end