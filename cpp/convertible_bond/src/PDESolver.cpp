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
        t = n*this->pde_param_ptr->dt;
        rate1 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->dt);
        rate2 = this->interp(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->dt);
        r = -(log(rate2/rate1)/this->pde_param_ptr->dt);

        if (ceil(t-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dt-this->pde_param_ptr->dtime))==1)
        {
            count6 += 1;
            for(int i =0;i<this->pde_param_ptr->B.size();i++)
            {
                    this->pde_param_ptr->B[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->B,(1-this->pde_param_ptr->q)*this->pde_param_ptr->S
                    );
            }
            for(int i =0;i<this->pde_param_ptr->u.size();i++)
            {
                    this->pde_param_ptr->u[i] = this->interp(
                    this->pde_param_ptr->S,this->pde_param_ptr->u,(1-this->pde_param_ptr->q)*this->pde_param_ptr->S
                    );
            }
        }

        vector<double> alpha;
        for (int i =2;i<=this->pde_param_ptr->Ns;i++)
        {
            double tmp_1 = power(this->pde_param_ptr->sigma,2)*power(this->pde_param_ptr->S[i],2)/(2*power(this->pde_param_ptr->ds,2));
            double tmp_2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]/(2*this->pde_param_ptr->ds);
            double tmp_3 = (tmp_1 - tmp_2)*this->pde_param_ptr->dt;
            alpha.push_back(tmp_3);
        }

        vector<double> beta;
        for (int i =2;i<=this->pde_param_ptr->Ns;i++)
        {
            double tmp_1 = power(this->pde_param_ptr->sigma,2)*power(this->pde_param_ptr->S[i],2)/(2*power(this->pde_param_ptr->ds,2));
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
        alpha_MuMB.push_back(0);
        alpha_MuMB.push_back(0);
        for (int i=0;i<beta.size();i++)
        {
            beta_MuMB.push_back(beta[i]);
        }
        beta_MuMB.push_back(0);
        beta_MuMB.push_back(0);
        gamma_Mu=[-(r+p)*dt;-(alpha+beta+(r+p)*dt);0];
        vector<double> gamma_Mu;

    }
}