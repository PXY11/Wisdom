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
    //具体待实现
    return 1;
}


void PDESolver::solve()
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
        if(flag=1)
            break;
        double t = n*this->pde_param_ptr->dt;
        double rate1 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        double rate2 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        double r = -(log(rate2/rate1)/this->pde_param_ptr->dt);

        if ((ceil(t-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dt-this->pde_param_ptr->dtime))==1)
        {
            count6 += 1;
            for(int i =0;i<this->pde_param_ptr->B.size();i++)
            {
                    this->pde_param_ptr->B[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->B,1.998
                    );
            }
            for(int i =0;i<this->pde_param_ptr->u.size();i++)
            {
                    this->pde_param_ptr->u[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->u,1.998
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
            double tmp_3 = (tmp_1 - tmp_2)*this->pde_param_ptr->dt;
            beta.push_back(tmp_3);
        }
        vector<double> alpha_MuMB;
        vector<double> beta_MuMB;
        for (int i=0;i<alpha.size();i++)
        {
            alpha_MuMB.push_back(alpha[i]);
        }
        for (int i = alpha_MuMB.size();i<alpha.size();i++)
        {
            alpha_MuMB.push_back(0); //alpha_MuMB不够alpha长度的部分补0
        }

        for (int i=0;i<beta.size();i++)
        {
            beta_MuMB.push_back(beta[i]); 
        }
        for (int i = beta_MuMB.size();i<beta.size();i++)
        {
            beta_MuMB.push_back(0); //beta_MuMB不够beta长度的部分补0
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
        vector<double> u;


    }
}






