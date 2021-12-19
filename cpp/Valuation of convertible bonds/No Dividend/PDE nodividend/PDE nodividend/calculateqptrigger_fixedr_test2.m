function difference = calculateqptrigger_fixedr_test2(S0, Shistory, ptrigger, qptrigger,T,r,sigma,no_put_time,n,m)
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
sam1=random('norm',0,sqrt(d_t),N,M);
sam2=-sam1;
sam=[sam1;sam2];
W=cumsum(sam,2);
for i=1:2*N
    S=S0*exp((r-sigma^2/2)*(1:M)*d_t+sigma*W(i,:));
    SS=[Shistory,S0,S];
    
    if T+(m-1)*d_t<=no_put_time(1)-d_t
        judge=(SS<ptrigger);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m
            total_left=total_left+1;
        elseif ~isempty(left_j)
            disp('1 the convertible bond was already put')
        end
    elseif T+(m-1)*d_t>=no_put_time(1) && T<=no_put_time(1)-d_t
        left_delta=T+(m-1)*d_t-no_put_time(1);
        SS=SS(floor(left_delta/d_t+1e-10)+2:end);
        judge=(SS<ptrigger);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if left_j>m-1-floor(left_delta/d_t+1e-10)
            total_left=total_left+1;
        elseif ~isempty(left_j)
            disp('2 the convertible bond was already put')
        end      
    else
        %T>no_put_time(1)-d_t
        left_delta=T-no_put_time(1)+d_t;
        SS=S(ceil(left_delta/d_t-1e-10):end);
        judge=(SS<ptrigger);
        judge_filter=filter(filter_b,filter_a,judge);
        left_j=find(judge_filter>=n,1);
        if ~isempty(left_j)
            total_left=total_left+1;
        end
    end
    
    
    if T<=no_put_time(1)-d_t
        SS=[S0,S];
        right_j=find(SS<qptrigger,1);
        if right_j>1
            total_right=total_right+1;
        end
    else
        %T>no_put_time(1)-d_t
        right_delta=T-no_put_time(1)+d_t;
        SS=S(ceil(right_delta/d_t-1e-10):end);
        right_j=find(SS<qptrigger,1);
        if ~isempty(right_j)
            total_right=total_right+1;
        end
    end
end
difference=(total_left-total_right)/N/2;
end