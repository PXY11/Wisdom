clc   %clear command window
close all  %close all figures
clear  %clear all variables in workspace

%这里的call和put是条款里真实的call和put

%diary('C:\Users\LX\Desktop\log.txt');
%diary on;

% format long;

clc;
clear;
tic;
format long;
%%

fid=fopen('MyFile.txt','w');

%%
X_1=5.73;
F=100;
dt=1;
coupon_rate=[1.8 1.5 1.5 1 0.7 0.5]/100/252;   %coupon rates
bonus_rate=4.2/100;
%data=xlsread('data',1,'C13:R252');
% T=126;
% r=zeros(1,T);  %risk free rate
% r(1:4)=2.891104+0.132771*(0:3);
% r(5:62)=3.400345+0.002063*(0:57);
% r(63:126)=3.519187+0.000454*(0:63);

no_convert_time=[5.5 6]*252;
no_call_time=[5.5 6]*252;
no_put_time=[2 6]*252;
call_m=30;
call_n=15;
put_m=30;
put_n=30;
call_Shistory=[6.25 6.13 6.23 6.03 6.19 6.32 6.39 6.33 6.20 6.04 6.13 6.12 6.11 6.25 6.19 6.09 6.01 6.13 6.13 6.24 6.33 6.80 7.48 7.70 8.30 7.75 7.77 7.81 7.61];
put_Shistory=[6.25 6.13 6.23 6.03 6.19 6.32 6.39 6.33 6.20 6.04 6.13 6.12 6.11 6.25 6.19 6.09 6.01 6.13 6.13 6.24 6.33 6.80 7.48 7.70 8.30 7.75 7.77 7.81 7.61];

for case_T=1:3 % corresponding to three different cases
    
    % three different cases for Shistory
    if case_T == 1
        T=round(252*5.7);
        fprintf(fid, 'T = %f\n', round(T/252)');
    end
    if case_T == 2
        T=round(252*5.45);
        fprintf(fid, 'T = %f\n', round(T/252)');
    end
    if case_T == 3
        T=round(252*4.6);
        fprintf(fid, 'T = %f\n', round(T/252)');
    end
    
    r=zeros(1,T);  %risk free rate
    r(1:T)=3.5;
    r=r/100/252;
    fr=zeros(1,T);
    fr(1)=r(1);
    for j=2:T
        fr(j)=log(exp(r(j)*j*dt)/exp(r(j-1)*(j-1)*dt))/dt;   %forward rate
    end
    N=5;
    K=zeros(1,T);
    K(T)=F*coupon_rate(1)*252;
    K(T-(1:floor(T/252))*252+1)=F*coupon_rate(2:floor(T/252)+1)*252;
    AccI=zeros(1,T);
    for k=1:floor(T/252+1e-10)
        AccI(T-251-252*(k-1):T-252*(k-1))=F*coupon_rate(k)*(1:252);
    end
    Bcclean=F;  %clean call price
    Bpclean=F;  %clean put price
    Bc=Bcclean+AccI;
    Bp=Bpclean+AccI;
    filter_a=1;
    call_filter_b=ones(1,call_m);
    put_filter_b=ones(1,put_m);
    price=zeros(1,22);
    
    for case_D = 1:5
        switch case_D
            case 1
                dividend1 = 0.3;
                dividend2 = 0.3;
                dividend3 = 0.3;
                fprintf(fid, 'dividend1 = %f\n', dividend1');
                fprintf(fid, 'dividend2 = %f\n', dividend2');
                fprintf(fid, 'dividend3 = %f\n', dividend3');
            case 2
                dividend1 = 0.4;
                dividend2 = 0.4;
                dividend3 = 0.4;
                fprintf(fid, 'dividend1 = %f\n', dividend1');
                fprintf(fid, 'dividend2 = %f\n', dividend2');
                fprintf(fid, 'dividend3 = %f\n', dividend3');
            case 3
                dividend1 = 0.5;
                dividend2 = 0.5;
                dividend3 = 0.5;
                fprintf(fid, 'D1 = %f\n', dividend1');
                fprintf(fid, 'D2 = %f\n', dividend2');
                fprintf(fid, 'D3 = %f\n', dividend3');
            case 4
                dividend1 = 1.0;
                dividend2 = 1.0;
                dividend3 = 1.0;
                fprintf(fid, 'D1 = %f\n', dividend1');
                fprintf(fid, 'D2 = %f\n', dividend2');
                fprintf(fid, 'D3 = %f\n', dividend3');
            case 5
                dividend1 = 1.5;
                dividend2 = 1.5;
                dividend3 = 1.5;
                fprintf(fid, 'D1 = %f\n', dividend1');
                fprintf(fid, 'D2 = %f\n', dividend2');
                fprintf(fid, 'D3 = %f\n', dividend3');
            case 6
                dividend1 = 2.0;
                dividend2 = 2.0;
                dividend3 = 2.0;
                fprintf(fid, 'D1 = %f\n', dividend1');
                fprintf(fid, 'D2 = %f\n', dividend2');
                fprintf(fid, 'D3 = %f\n', dividend3');
        end
        
        for case_Dtime=1:3
            
            switch case_Dtime
            case 1
                dividend1_time = 1.0*252;
                dividend2_time = 3.0*252;
                dividend3_time = 5.0*252;
                fprintf(fid, 'dividend1_time = %f\n', dividend1_time/252');
                fprintf(fid, 'dividend2_time = %f\n', dividend2_time/252');
                fprintf(fid, 'dividend3_time = %f\n', dividend3_time/252');
            case 2
                dividend1_time = 3.4*252;
                dividend2_time = 3.5*252;
                dividend3_time = 3.6*252;
                fprintf(fid, 'dividend1_time = %f\n', dividend1_time/252');
                fprintf(fid, 'dividend2_time = %f\n', dividend2_time/252');
                fprintf(fid, 'dividend3_time = %f\n', dividend3_time/252');
            case 3
                dividend1_time = 3.49*252;
                dividend2_time = 3.5*252;
                dividend3_time = 3.51*252;
                fprintf(fid, 'dividend1_time = %f\n', dividend1_time/252');
                fprintf(fid, 'dividend2_time = %f\n', dividend2_time/252');
                fprintf(fid, 'dividend3_time = %f\n', dividend3_time/252');
            end
            
            for case_sigma=1:7
                switch case_sigma
                    case 1
                        sig = 0.2/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 2
                        sig = 0.4/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 3
                        sig = 0.6/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 4
                        sig = 0.8/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 5
                        sig = 1.0/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 6
                        sig = 1.5/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                    case 7
                        sig = 2.5/sqrt(252);
                        fprintf(fid, 'sigma = %f\n', sig * sqrt(252)');
                end
                
                for case_S0 = 1:6
                    switch case_S0
                        case 1
                            S0 = X_1*1.28;
                            fprintf(fid, 'S0 = %f\n', S0');
                        case 2
                            S0 = X_1*1.30;
                            fprintf(fid, 'S0 = %f\n', S0');
                        case 3
                            S0 = X_1*1.32;
                            fprintf(fid, 'S0 = %f\n', S0');
                        case 4
                            S0 = X_1*0.68;
                            fprintf(fid, 'S0 = %f\n', S0');
                        case 5
                            S0 = X_1*0.70;
                            fprintf(fid, 'S0 = %f\n', S0');
                        case 6
                            S0 = X_1*0.72;
                            fprintf(fid, 'S0 = %f\n', S0');
                    end
                    tic;
                    for jj=1:20
                        cash=zeros(2*N,T);
                        sample=zeros(2*N,T);
                        conversion_price=zeros(2*N,T);
                        %     call_X=zeros(2*N,T+call_m);
                        %     put_X=zeros(2*N,T+put_m);
                        %     call_trigger=zeros(2*N,T+call_m);
                        %     put_trigger=zeros(2*N,T+put_m);
                        call_triggertime=zeros(2*N,1);
                        rng('shuffle');
                        ran1=randn(N,T);
                        ran2=-ran1;
                        ran=[ran1;ran2];
                        call_mark=[];
                        put_mark=[];
                        count=0;
                        for i=1:2*N
                            sample(i,1)=S0*exp((fr(1)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,1));
                            for j=2:dividend1_time-1
                                sample(i,j)=sample(i,j-1).*exp((fr(j)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,j));
                            end
                            sampletmp1=sample(i,dividend1_time-1).*exp((fr(dividend1_time)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,dividend1_time));
                            sample(i,dividend1_time)=sampletmp1-dividend1;
                            for j=dividend1_time+1:dividend2_time-1
                                sample(i,j)=sample(i,j-1).*exp((fr(j)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,j));
                            end
                            sampletmp2=sample(i,dividend2_time-1).*exp((fr(dividend2_time)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,dividend2_time));
                            sample(i,dividend2_time)=sampletmp2-dividend2;
                            for j=dividend2_time+1:dividend3_time-1
                                sample(i,j)=sample(i,j-1).*exp((fr(j)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,j));
                            end
                            sampletmp3=sample(i,dividend3_time-1).*exp((fr(dividend3_time)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,dividend3_time));
                            sample(i,dividend3_time)=sampletmp3-dividend3;
                            for j=dividend3_time+1:T
                                sample(i,j)=sample(i,j-1).*exp((fr(j)-sig^2/2)*dt+sig*sqrt(dt)*ran(i,j));
                            end
                            
                            call_SS=[call_Shistory,S0,sample(i,:)];
                            put_SS=[put_Shistory,S0,sample(i,:)];
                            
                            call_conversion_price(1:dividend1_time-1+call_m)=5.73;   %conversion price
                            call_conversion_price(dividend1_time+call_m:dividend2_time-1+call_m)=5.73-dividend1;
                            call_conversion_price(dividend2_time+call_m:dividend3_time-1+call_m)=5.73-dividend1-dividend2;
                            call_conversion_price(dividend3_time+call_m:T+call_m)=5.73-dividend1-dividend2-dividend3;
                            
                            put_conversion_price(1:dividend1_time-1+put_m)=5.73;
                            put_conversion_price(dividend1_time+put_m:dividend2_time-1+put_m)=5.73-dividend1;
                            put_conversion_price(dividend2_time+put_m:dividend3_time-1+put_m)=5.73-dividend1-dividend2;
                            put_conversion_price(dividend3_time+put_m:T+put_m)=5.73-dividend1-dividend2-dividend3;
                            
                            conversion_price(i,1:dividend1_time-1)=5.73;
                            conversion_price(i,dividend1_time:dividend2_time-1)=5.73-dividend1;
                            conversion_price(i,dividend2_time:dividend3_time-1)=5.73-dividend1-dividend2;
                            conversion_price(i,dividend3_time:T)=5.73-dividend1-dividend2-dividend3;
                            
                            call_trigger=1.3*call_conversion_price;
                            put_trigger=0.7*put_conversion_price;
                            
                            call_SS=call_SS/conversion_price(1,1);
                            put_SS=put_SS/conversion_price(1,1);
                            sample(i,:)=sample(i,:)/conversion_price(1,1);
                            
                            
                            %%    call
                            if T+(call_m-1)*dt<=no_call_time(1)-dt
                                judge=(call_SS>=call_trigger/conversion_price(1,1));
                                judge_filter=filter(call_filter_b,filter_a,judge);
                                call_j=find(judge_filter>=call_n,1);
                                if call_j>call_m
                                    call_triggertime(i)=(call_j-call_m)*dt;
                                    call_payoff=F/conversion_price(i,call_triggertime(i))*call_SS(call_j);
                                    cash(i,call_triggertime(i))=call_payoff;
                                    call_mark=[call_mark i];
                                elseif ~isempty(call_j)
                                    disp('1 the convertible bond was already called')
                                end
                            elseif T+(call_m-1)*dt>=no_call_time(1) && T<=no_call_time(1)-dt
                                call_delta=T+(call_m-1)*dt-no_call_time(1);
                                call_SS=call_SS(floor(call_delta/dt+1e-10)+2:end);
                                judge=(call_SS>=call_trigger(floor(call_delta/dt+1e-10)+2:end)/conversion_price(1,1));
                                judge_filter=filter(call_filter_b,filter_a,judge);
                                call_j=find(judge_filter>=call_n,1);
                                if call_j>call_m-1-floor(call_delta/dt+1e-10)
                                    call_triggertime(i)=(call_j-call_m+1+floor(call_delta/dt+1e-10))*dt;
                                    call_payoff=F/conversion_price(i,call_triggertime(i))*call_SS(call_j);
                                    cash(i,call_triggertime(i))=call_payoff;
                                    call_mark=[call_mark i];
                                elseif ~isempty(call_j)
                                    disp('2 the convertible bond was already called')
                                end
                            else
                                %T>no_call_time(1)-dt
                                call_delta=T-no_call_time(1)+dt;
                                call_SS=sample(i,ceil(call_delta/dt-1e-10):end);
                                judge=(call_SS>=call_trigger(ceil(call_delta/dt-1e-10)+call_m:end)/conversion_price(1,1));
                                judge_filter=filter(call_filter_b,filter_a,judge);
                                call_j=find(judge_filter>=call_n,1);
                                if ~isempty(call_j)
                                    call_triggertime(i)=(call_j-1+ceil(call_delta/dt-1e-10))*dt;
                                    call_payoff=F/conversion_price(i,call_triggertime(i))*call_SS(call_j);
                                    cash(i,call_triggertime(i))=call_payoff;
                                    call_mark=[call_mark i];
                                end
                            end
                            %      put
                            if T+(put_m-1)*dt<=no_put_time(1)-dt
                                judge=(put_SS<put_trigger/conversion_price(1,1));
                                judge_filter=filter(put_filter_b,filter_a,judge);
                                put_j0=find(judge_filter>=put_n);
                                if put_j0(end)>put_m
                                    index=find(put_j0>put_m);
                                    put_triggertime=(put_j0(index)-put_m)*dt;
                                    count=count+1;
                                    put_mark(i,1:length(put_triggertime))=put_triggertime;
                                elseif ~isempty(put_j0)
                                    disp('1 the convertible bond was already put')
                                end
                            elseif T+(put_m-1)*dt>=no_put_time(1) && T<=no_put_time(1)-dt
                                put_delta=T+(put_m-1)*dt-no_put_time(1);
                                put_SS=put_SS(floor(put_delta/dt+1e-10)+2:end);
                                judge=(put_SS<put_trigger(floor(put_delta/dt+1e-10)+2:end)/conversion_price(1,1));
                                judge_filter=filter(put_filter_b,filter_a,judge);
                                put_j0=find(judge_filter>=put_n,1);
                                if put_j0(end)>put_m-1-floor(put_delta/d_t+1e-10)
                                    index=find(put_j0>put_m-1-floor(put_delta/d_t+1e-10));
                                    put_triggertime=(put_j0(index)-put_m+1+floor(put_delta/dt+1e-10))*dt;
                                    count=count+1;
                                    put_mark(i,1:length(put_triggertime))=put_triggertime;
                                elseif ~isempty(put_j0)
                                    disp('2 the convertible bond was already put')
                                end
                            else
                                %T>no_put_time(1)-dt
                                put_delta=T-no_put_time(1)+dt;
                                put_SS=sample(i,ceil(put_delta/dt-1e-10):end);
                                judge=(put_SS<put_trigger(ceil(put_delta/dt-1e-10)+put_m:end)/conversion_price(1,1));
                                judge_filter=filter(put_filter_b,filter_a,judge);
                                put_j0=find(judge_filter>=put_n,1);
                                if ~isempty(put_j0)
                                    put_triggertime=(put_j0-1+ceil(put_delta/dt-1e-10))*dt;
                                    count=count+1;
                                    put_mark(i,1:length(put_triggertime))=put_triggertime;
                                end
                            end
                        end
                        %%
                        nocall=setdiff(1:2*N,call_mark);
                        cash(nocall,T)=max(F./conversion_price(nocall,T).*sample(nocall,T),F*(1+bonus_rate)/conversion_price(1,1)+K(T)/conversion_price(1,1));
                        
                        for j=T-1:-1:1
                            %j=dividend_time;
                            afterj=find(call_triggertime>j);   %find all paths that are called after time j
                            paths=union(nocall,afterj);
                            
                            x=sample(paths,j);
                            y=zeros(size(x));
                            for ii=1:length(paths)
                                pathi=paths(ii);
                                mark=find(cash(pathi,:)>0);
                                Fr=log(exp(r(mark)*mark*dt)/exp(r(j)*j*dt))/dt/(mark-j);   %forward rate
                                y(ii)=cash(pathi,mark)*exp(-Fr*(mark-j)*dt);
                                for kk=mark-1:-1:j+1
                                    Fr=((1+r(kk)*kk*dt)/(1+r(j)*j*dt)-1)/dt/(kk-j);
                                    y(ii)=y(ii)+K(kk)/conversion_price(1,1)*exp(-Fr*(kk-j)*dt);
                                end
                                y(ii)=y(ii)+K(j)/conversion_price(1,1);
                            end
                            
                            %Legendre polynomials
                            %%X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8 (63*x.^5-70*x.^3+15*x)/8 (231*x.^6-315*x.^4+105*x.^2-5)/16 (429*x.^7-693*x.^5+315*x.^3-35*x)/16 (6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35)/128 (12155*x.^9-25740*x.^7+18018*x.^5-4620*x.^3+315*x)/128];
                            %X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8 (63*x.^5-70*x.^3+15*x)/8 (231*x.^6-315*x.^4+105*x.^2-5)/16 (429*x.^7-693*x.^5+315*x.^3-35*x)/16 (6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35)/128];
                            %X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8 (63*x.^5-70*x.^3+15*x)/8 (231*x.^6-315*x.^4+105*x.^2-5)/16 (429*x.^7-693*x.^5+315*x.^3-35*x)/16];
                            %X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8 (63*x.^5-70*x.^3+15*x)/8 (231*x.^6-315*x.^4+105*x.^2-5)/16];
                            X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8 (63*x.^5-70*x.^3+15*x)/8];
                            %X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2 (35*x.^4-30*x.^2+3)/8];
                            %X=[ones(size(x)) x (3*x.^2-1)/2 (5*x.^3-3*x)/2];
                            %X=[ones(size(x)) x (3*x.^2-1)/2];
                            %X=[ones(size(x)) x];
                            
                            b=X\y;
                            v=X*b;
                            tmp1=find(v<F./conversion_price(paths,j).*sample(paths,j)/conversion_price(1,1));
                            count1=length(tmp1);
                            if count1>0
                                disp([num2str(count1) ' paths convert at j=' num2str(j)])
                            end
                            cash(paths(tmp1),j)=F./conversion_price(paths(tmp1),j).*sample(paths(tmp1),j)/conversion_price(1,1);
                            cash(paths(tmp1),j+1:T)=0;
                            tmp2=find(v<Bp(j)/conversion_price(1,1));
                            %         count2=length(tmp2);
                            %         if count2>0
                            %             disp([num2str(count2) ' paths want to put at j=' num2str(j)])
                            %         end
                            count3=0;
                            for mm=1:length(tmp2)
                                if ismember(j,put_mark(tmp2(mm),:))
                                    cash(paths(tmp2(mm)),j)=Bp(j)/conversion_price(1,1);
                                    cash(paths(tmp2(mm)),j+1:T)=0;
                                    count3=count3+1;
                                end
                            end
                            if count3>0
                                disp([num2str(count3) ' paths actually put at j=' num2str(j)])
                            end
                        end
                        
                        flag=cash>0;
                        sum0=sum(flag,2);
                        mark1=find(sum0==0);
                        if ~isempty(mark1)
                            disp(mark1')
                            dis('paths have no cash flow')
                        end
                        mark2=find(sum0>1);
                        if ~isempty(mark2)
                            disp(mark2')
                            disp('paths have more than one cash flows')
                        end
                        
                        for i=1:2*N
                            mark3=find(cash(i,:)>0,1);
                            cash(i,1:mark3-1)=cash(i,1:mark3-1)+K(1:mark3-1)/conversion_price(1,1);
                        end
                        
                        C=zeros(2*N,1);
                        for j=T:-1:1
                            %dr(j)=r(j)+4/100/250;
                            
                            dr(j)=r(j);
                            C=C+cash(:,j)*exp(-dr(j)*j*dt);
                        end
                        V=sum(C)*conversion_price(1,1)/N/2;
                        disp([jj V]);
                        price(jj)=V;
                    end
                    price(21)=mean(price(1:20));
                    price(22)=std(price(1:20));
                    disp(price(21:22));
                    fprintf(fid, 'price = %f\n', price(21)');
                    fprintf(fid, 'std = %f\n', price(22)');
                    fprintf(fid, 'Total elapsed time: %f \n', toc);
                end % case_S0 - 6 
            end % case sigma - 7
        end % case Dtime - 3
    end % case D - 6
end % three different cases - 3
    
fclose(fid);

