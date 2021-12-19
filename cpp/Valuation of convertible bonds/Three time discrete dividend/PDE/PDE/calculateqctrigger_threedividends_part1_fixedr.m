function qctrigger = calculateqctrigger_threedividends_part1_fixedr(S0, Shistory, ctrigger,T,r,sigma,no_call_time,n,m,D1_time)
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
countcall=0;
%countcall2=0;
filter_a=1;
filter_b=ones(1,m);
rng('shuffle');
%sam1=random('norm',0,1,N,M);
sam1=randn(N,M);
sam2=-sam1;
sam=[sam1;sam2];
%S=zeros(1,D1_time/d_t-1);
S=zeros(1,M);
for i=1:2*N
    S(1)=S0*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,1));
    for j=2:D1_time/d_t-1
        S(j)=S(j-1).*exp((r-sigma^2/2)*d_t+sigma*sqrt(d_t)*sam(i,j));
    end
    SS=[Shistory,S0,S];
    
    %ctrigger_sequence=zeros(1,D1_time/d_t-1+m);
    ctrigger_sequence=zeros(1,M+m);
    ctrigger_sequence(1:end)=ctrigger;

    if T+(m-1)*d_t<=no_call_time(1)-d_t
        judge=(SS>=ctrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m
            left_Smark=SS(left_j);
            countcall=countcall+1;
            total_left_Smark=total_left_Smark+left_Smark;
        elseif ~isempty(left_j)
            disp('1 the convertible bond was already called')
        end
    elseif T+(m-1)*d_t>=no_call_time(1) && T<=no_call_time(1)-d_t
        left_delta=T+(m-1)*d_t-no_call_time(1);
        SS=SS(floor(left_delta/d_t+1e-10)+2:end);
        ctrigger_sequence=ctrigger_sequence(floor(left_delta/d_t+1e-10)+2:end);
        judge=(SS>=ctrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m-1-floor(left_delta/d_t+1e-10)
            left_Smark=SS(left_j);
            countcall=countcall+1;
            total_left_Smark=total_left_Smark+left_Smark;
        elseif ~isempty(left_j)
            disp('2 the convertible bond was already called')
        end      
    else
        %T>no_call_time(1)-d_t
        left_delta=T-no_call_time(1)+d_t;
        SS=S(ceil(left_delta/d_t-1e-10):end);
        ctrigger_sequence=ctrigger_sequence(ceil(left_delta/d_t-1e-10)+m:end);
        judge=(SS>=ctrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if ~isempty(left_j)
            left_Smark=SS(left_j);
            countcall=countcall+1;
            total_left_Smark=total_left_Smark+left_Smark;
        end
    end
end
qctrigger=total_left_Smark/countcall;
%qctrigger2=(total_left_Smark1+total_left_Smark2)/(countcall1+countcall2);
end