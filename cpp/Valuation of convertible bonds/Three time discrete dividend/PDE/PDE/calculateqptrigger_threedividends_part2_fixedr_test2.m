function difference = calculateqptrigger_threedividends_part2_fixedr_test2(S0, Shistory, ptrigger, qptrigger,T,r,sigma,no_put_time,n,m,npX,D1_time,D1,D2_time)
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
total_left=0;
total_right=0;
filter_a=1;
filter_b=ones(1,m);
rng('shuffle');
sam1=randn(N,M);
sam2=-sam1;
sam=[sam1;sam2];
S=zeros(1,M);
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
    %S=10*ones(1,M);
    SS=[Shistory,S0,S];
    
    ptrigger_sequence=zeros(1,M+m);
    ptrigger_sequence(1:D1_time/d_t-1+m)=ptrigger;
    ptrigger_sequence(D1_time/d_t+m:end)=ptrigger-npX*D1;
    
    if T+(m-1)*d_t<=no_put_time(1)-d_t
        judge=(SS<ptrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m
            left_triggertime=(left_j-m)*d_t;
            if left_triggertime>=D1_time  %·Öºìºó
                total_left=total_left+1;
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
            if left_triggertime>=D1_time
                total_left=total_left+1;
            end
        elseif ~isempty(left_j)
            disp('2 the convertible bond was already put')
        end      
    else
        %T>no_put_time(1)-d_t
        left_delta=T-no_put_time(1)+d_t;
        SS=S(ceil(left_delta/d_t-1e-10):end);
        ptrigger_sequence=ptrigger_sequence(ceil(left_delta/d_t-1e-10)+m:end);
        judge=(SS<ptrigger_sequence);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if ~isempty(left_j)
            left_triggertime=(left_j-1+ceil(left_delta/d_t-1e-10))*d_t;
            if left_triggertime>=D1_time
                total_left=total_left+1;
            end
        end
    end
    
    if T<=no_put_time(1)-d_t
        SS=[S0,S];
        right_j=find(SS<qptrigger,1);
        if right_j>1
            right_triggertime=(right_j-1)*d_t;
            if right_triggertime>=D1_time
                total_right=total_right+1;
            end
        end
    else
        %T>no_call_time(1)-d_t
        right_delta=T-no_put_time(1)+d_t;
        SS=S(ceil(right_delta/d_t-1e-10):end);
        right_j=find(SS<qptrigger,1);
        if ~isempty(right_j)
            right_triggertime=(right_j-1+ceil(right_delta/d_t-1e-10))*d_t;
            if right_triggertime>=D1_time
                total_right=total_right+1;
            end
        end
    end
end
difference=(total_left-total_right)/N/2;
end