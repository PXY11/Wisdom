#include<PDESolver.h>
#include<iostream>
#include<cmath>
using namespace std;


void PDESolver::set_iter(int iter_time)
{
    this->iteration = iter_time;
}

void PDESolver::show_iter()
{
    cout<<"iter time = "<<this->iteration<<endl;
}

Parameter* PDESolver::get_pde_param_ptr()
{
    this->pde_param_ptr = new Parameter(11);
    return this->pde_param_ptr;
}

double PDESolver::interp(vector<double> x,vector<double> y,double ind)
{
    
    int token = 0;
    int i;
    for(i =0;i<ind;i++)
        token ++;
    double res = y[i-1] + (x[i] - x[i-1])/(y[i]-y[i-1]);
    return res;
}


double PDESolver::solve()
{
    int count = 0;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    int count5 = 0;
    int count6 = 0;
    int flag = 0;
    for (int n=1; n<=this->pde_param_ptr->Nt; n++)
    {
        if(flag==1)
            {break;}
        double t = n*this->pde_param_ptr->dt;
        //interp 具体待实现
        double rate1 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        double rate2 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        double r = -(log(rate2/rate1)/this->pde_param_ptr->dt);

        if ((ceil(t-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dt-this->pde_param_ptr->dtime))==1)
        {
            count6 += 1;
            for(int i =0;i<this->pde_param_ptr->B.size();i++)
            {
                    this->pde_param_ptr->B[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->B,2
                    );
            }
            for(int i =0;i<this->pde_param_ptr->u.size();i++)
            {
                    this->pde_param_ptr->u[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->u,2
                    );
            }
        }

        vector<double> alpha;
        for (int i =2;i<=this->pde_param_ptr->Ns;i++)
        {
            double tmp_1 = pow(this->pde_param_ptr->sigma,2)*pow(this->pde_param_ptr->S[i],2)/(2*pow(this->pde_param_ptr->ds,2));
            double tmp_2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]/(2*this->pde_param_ptr->ds);
            double tmp_3 = (tmp_1 - tmp_2)*this->pde_param_ptr->dt;
            alpha.push_back(tmp_3);
        }

        vector<double> beta;
        for (int i =2;i<=this->pde_param_ptr->Ns;i++)
        {
            double tmp_1 = pow(this->pde_param_ptr->sigma,2)*pow(this->pde_param_ptr->S[i],2)/(2*pow(this->pde_param_ptr->ds,2));
            double tmp_2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]/(2*this->pde_param_ptr->ds);
            double tmp_3 = (tmp_1 + tmp_2)*this->pde_param_ptr->dt;
            beta.push_back(tmp_3);
        }
        vector<double> alpha_MuMB;
        vector<double> beta_MuMB;
        for (int i=0;i<alpha.size();i++)
        {
            alpha_MuMB.push_back(alpha[i]);
        }
        for (int i = alpha_MuMB.size();i<(this->pde_param_ptr->Ns+1);i++)
        {
            alpha_MuMB.push_back(0); //alpha_MuMB不够alpha长度的部分补0
        }

        for (int i=0;i<(this->pde_param_ptr->Ns+1-beta.size());i++)
        {
            beta_MuMB.push_back(0); //beta_MuMB前面不够长的部分补0
        }
        for (int i=beta_MuMB.size()-1;i<beta.size();i++)
        {
            beta_MuMB.push_back(beta[i]); 
        }
        

        vector<double> gamma_Mu;
        double gamma_Mu_first = -(r+this->pde_param_ptr->p)*this->pde_param_ptr->dt;
        gamma_Mu.push_back(gamma_Mu_first);
        for (int i=0;i<alpha.size();i++)
        {
            double tmp = alpha[i]+beta[i]+gamma_Mu_first*(-1);
            gamma_Mu.push_back(tmp);
        }
        gamma_Mu.push_back(0);

        vector<double> gamma_MB;
        double gamma_MB_first = -(r+this->pde_param_ptr->p)*(1-this->pde_param_ptr->R)*this->pde_param_ptr->dt;
        gamma_MB.push_back(gamma_Mu_first);
        for (int i=0;i<alpha.size();i++)
        {
            double tmp = alpha[i]+beta[i]+gamma_MB_first*(-1);
            gamma_MB.push_back(tmp);
        }
        gamma_MB.push_back(0);

        vector<double> Mu;
        vector<double> MB; //[alpha_MuMB,gamma_MB,beta_MuMB]

        double AccI = (ceil(t) - t)*this->pde_param_ptr->coupon_time_rate[ceil(t)]*this->pde_param_ptr->F;
        double Bc_n;
        if ((this->pde_param_ptr->no_call_time[0]<=t&&t<=this->pde_param_ptr->no_call_time[1]))
        {
            Bc_n = 1000;
        }
        else
        {
            Bc_n = this->pde_param_ptr->Bc_star + AccI;
        }

        double Bp_n;
        if ((this->pde_param_ptr->no_put_time[0]<=t)&&t<=this->pde_param_ptr->no_put_time[1])
        {
            Bp_n = 0.01;
        }
        else
        {
            Bp_n = this->pde_param_ptr->Bp_star + AccI;
        }

        vector<double> k_n;
        if (this->pde_param_ptr->no_convert_time[0]<=t && t<=this->pde_param_ptr->no_put_time[1])
        {
            for(int i =0;i<this->pde_param_ptr->S.size();i++)
            {
                k_n.push_back(0.01*1);
            }
        }
        else 
        {
            for(int i=0;i<this->pde_param_ptr->S.size();i++)
            {
                double tmp;
                tmp = this->pde_param_ptr->F/(this->pde_param_ptr->conversion_price-this->pde_param_ptr->q*this->pde_param_ptr->S[i]*(ceil(this->pde_param_ptr->T-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dtime)));
            }
        }

        vector<double> u_n;
        for(int i=0;i<this->pde_param_ptr->u.size();i++)
        {
            u_n.push_back(this->pde_param_ptr->u[i]);
        }

        //设置B
        for (int i=0;i<gamma_MB.size();i++)
        {
            double tmp = (1 - this->pde_param_ptr->theta*gamma_MB[i])/((1+(1-this->pde_param_ptr->theta)*this->pde_param_ptr->theta*gamma_MB[i])*this->pde_param_ptr->B[i]);
            this->pde_param_ptr->B[i] = tmp;
        }

        for (int i=0;i<this->pde_param_ptr->B.size();i++)
        {
            if(this->pde_param_ptr->B[i]>Bc_n)
            {
                this->pde_param_ptr->B[i] = Bc_n;
            }    
        }
        vector<double> Muu;
        for(int i=0;i<this->pde_param_ptr->S.size();i++)
        {
            double tmp_1 = max(k_n[i]*(1-this->pde_param_ptr->eta)*this->pde_param_ptr->S[i],this->pde_param_ptr->R*this->pde_param_ptr->B[i]);
            double tmp_2 = this->pde_param_ptr->p*this->pde_param_ptr->dt*tmp_1;
            Muu.push_back(tmp_2);
        }
        vector<double> u_old;
        for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
        {
            double tmp = 1 - this->pde_param_ptr->theta*Mu[i];
            double tmp2  =  (1 - (1-this->pde_param_ptr->theta)*Mu[i])*u_n[i];
            double tmp3  = Muu[i];
            double res = tmp/tmp2 + tmp3;
            u_old.push_back(res);
        }

        vector<double> P1_old,P2_old;
        for (int i = 0;i< this->pde_param_ptr->Ns+1;i++)
        {
            double tmp;
            double tmp0 = k_n[i]*this->pde_param_ptr->S[i]+ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->F];
            double tmp2 = this->pde_param_ptr->rhopenltyput*max(Bp_n,tmp);
            double tmp3 = -1*this->pde_param_ptr->rhopenltycall*max(Bp_n,tmp);
            P1_old.push_back(tmp2);
            P2_old.push_back(tmp3);

        } 
        vector<double> P1,P2;
        vector<double> u;
        while(1)
        {
            count += 1;
            
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp = 1 - this->pde_param_ptr->theta*Mu[i];
                double tmp2  =  (1 - (1-this->pde_param_ptr->theta)*Mu[i])*u_n[i];
                double tmp3  = Muu[i];
                double res = tmp/tmp2 + tmp3;
                u.push_back(res);
            }
            
            for (int i = 0;i< this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                double tmp0 = k_n[i]*this->pde_param_ptr->S[i]+ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->F];
                double tmp2 = this->pde_param_ptr->rhopenltyput*max(Bp_n,tmp);
                double tmp3 = -1*this->pde_param_ptr->rhopenltycall*max(Bp_n,tmp);
                P1.push_back(tmp2);
                P2.push_back(tmp3);

            } 
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp0 = this->pde_param_ptr->Ns+1 - this->pde_param_ptr->theta*Mu[i] + P1[i] - P2[i];
                double tmp1 = this->pde_param_ptr->Ns+1 - this->pde_param_ptr->theta*Mu[i];
                double tmp2 = P1[i]*max(Bp_n,(k_n[i]*this->pde_param_ptr->S[i]+ceil(t+this->pde_param_ptr->dt)-ceil(t))*this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->F]);
                double tmp3 = -P2[i]*max(Bp_n,(k_n[i]*this->pde_param_ptr->S[i]+ceil(t+this->pde_param_ptr->dt)-ceil(t))*this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->F]);
                double res = tmp0/tmp1+tmp2+tmp3;
                u[i] = res;
            }
            for (int i = 0;i< this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                double tmp0 = k_n[i]*this->pde_param_ptr->S[i]+ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)*this->pde_param_ptr->F];
                double tmp2 = this->pde_param_ptr->rhopenltyput*max(Bp_n,tmp);
                double tmp3 = -1*this->pde_param_ptr->rhopenltycall*max(Bp_n,tmp);
                P1.push_back(tmp2);
                P2.push_back(tmp3);

            } 
            bool eq0 = true;
            for(int i=0;i<P1.size();i++)
            {
                if(P1[i]!=P1_old[i])
                    eq0=false;
            } 
            bool eq1 = true;
            for(int i=0;i<P2.size();i++)
            {
                if(P2[i]!=P2_old[i])
                    eq1 = false;
            }          
            if (eq1 && eq0)
            {
                count1 += 1;
                break;
            }
            vector<double> u_u_old ;
            u_u_old.push_back(0);
            vector<double> abs_u;
            for (int i=1;i<u.size();i++)
            {   
                double res = abs(u[i]-u_old[i]);
                if (res>u_u_old[i-1])
                    u_u_old.push_back(res);
                abs_u.push_back(abs(u[i]));
            }
            double tmp1 = u_u_old[-1];
            if ((tmp1/abs_u[-1]) < 1/0.0001)
            {
                count2 += 1;
                break;
            }
            if (count > this->pde_param_ptr->Nt*40)
            {
                count4 +=1;
                cout<<"Probably does not converge, current n is "<<n<<endl;
                flag = 1;
            }
            else 
                count3 += 1;
            for(int i =0;i<P1_old.size();i++)
            {
                
                P1_old[i] = P1[i];
            }

            for(int i =0;i<P2_old.size();i++)
            {
                
                P2_old[i] = P2[i];
            }

            if(ceil(this->pde_param_ptr->dt+t)-ceil(t)==1)
            {
                count5 +=1;
                for(int i=0;i<u.size();i++)
                {
                    u[i] = u[i] + this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)]*this->pde_param_ptr->F;
                    this->pde_param_ptr->B[i] = this->pde_param_ptr->B[i] +this->pde_param_ptr->coupon_time_rate[ceil(t+this->pde_param_ptr->dt)]*this->pde_param_ptr->F;

                }
            }
        }
        double Uindex = this->pde_param_ptr->Ns*this->pde_param_ptr->S0 + 1;
        double U;
        if (Uindex >1)
            U = u[ceil(Uindex)]*(Uindex-floor(Uindex+0.001))+u[floor(Uindex)]*(ceil(Uindex+0.001)-Uindex);
        else 
            U = u[ceil(Uindex)];
        return U;
    }

}

double PDESolver::calBpstar()
{
    int N = 10000;
    double d_t = 1/252;
    double M = ceil(this->pde_param_ptr->T/d_t);
    double eps=0.0001;
    vector<double> left_sum,right_sum;
    for(int i=0;i<N;i++)
    {
        left_sum.push_back(0);
        right_sum.push_back(0);
    }
    double windowSize = 30;
    double n = 30;
    vector<double>filter_b;
    for (int i=0;i<windowSize;i++)
    {
        filter_b.push_back(1/windowSize*1);
    }
    double filter_a = 1;
    double rate1,rate2;
    vector<double> r;
    for(int i=1;i<M;i++)
    {
        rate1 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,i*d_t);
        rate2 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(i-1)*d_t);
        r.push_back(-log(rate2/rate1)/d_t);
    }
    vector<double> W;
    for(int i=0;i<M;i++)
    {
        W.push_back(sqrt(d_t));
    }
    for(int i=0;i<N;i++)
    {
        double S = this->pde_param_ptr->S0*exp(r[i]-this->pde_param_ptr->sigma*this->pde_param_ptr->sigma/2)*d_t+this->pde_param_ptr->sigma*W[i];
        vector<bool> judge;
        vector<double>judgefilter;
        for(int i=0;i<this->pde_param_ptr->S.size();i++)
        {
            judge.push_back(this->pde_param_ptr->S[i]<this->pde_param_ptr->Bp_star);
            judgefilter.push_back(judge[i]+eps);
            if (judgefilter[i] >= n/M)
                left_sum.push_back(1);
            if (judgefilter[i] ==1)
                right_sum.push_back(1);
        }

    }
    double sum = 0;
    for(int i=0;i<left_sum.size();i++)
    {
        sum += left_sum[i];
    }
    double pbp = sum/left_sum.size();
    sum = 0;
    for(int i=0;i<right_sum.size();i++)
    {
        sum += left_sum[i];
    }
    double pbpstar = sum/right_sum.size();
    double pdiff = pbp - pbpstar;
    return pdiff;

}



double PDESolver::calBcstar()
{
    int N = 10000;
    double d_t = 1/252;
    double M = ceil(this->pde_param_ptr->T/d_t);
    double eps=0.0001;
    vector<double> left_sum,right_sum;
    for(int i=0;i<N;i++)
    {
        left_sum.push_back(0);
        right_sum.push_back(0);
    }
    double windowSize = 30;
    double n = 30;
    vector<double>filter_b;
    for (int i=0;i<windowSize;i++)
    {
        filter_b.push_back(1/windowSize*1);
    }
    double filter_a = 1;
    double rate1,rate2;
    vector<double> r;
    for(int i=1;i<M;i++)
    {
        rate1 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,i*d_t);
        rate2 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(i-1)*d_t);
        r.push_back(-log(rate2/rate1)/d_t);
    }
    vector<double> W;
    for(int i=0;i<M;i++)
    {
        W.push_back(sqrt(d_t));
    }
    for(int i=0;i<N;i++)
    {
        double S = this->pde_param_ptr->S0*exp(r[i]-this->pde_param_ptr->sigma*this->pde_param_ptr->sigma/2)*d_t+this->pde_param_ptr->sigma*W[i];
        vector<bool> judge;
        vector<double>judgefilter;
        for(int i=0;i<this->pde_param_ptr->S.size();i++)
        {
            judge.push_back(this->pde_param_ptr->S[i]<this->pde_param_ptr->Bp_star);
            judgefilter.push_back(judge[i]+eps);
            if (judgefilter[i] >= n/M)
                left_sum.push_back(1);
            if (judgefilter[i] ==1)
                right_sum.push_back(1);
        }

    }
    double sum = 0;
    for(int i=0;i<left_sum.size();i++)
    {
        sum += left_sum[i];
    }
    double pbc = sum/left_sum.size();
    sum = 0;
    for(int i=0;i<right_sum.size();i++)
    {
        sum += left_sum[i];
    }
    double pbcstar = sum/right_sum.size();
    double pdiff = pbc - pbcstar;
    return pdiff;

}

