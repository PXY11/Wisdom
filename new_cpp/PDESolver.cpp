#include"PDESolver.h"
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

//Parameter* PDESolver::get_pde_param_ptr()
//{
//    this->pde_param_ptr = new Parameter(11);
//    return this->pde_param_ptr;
//}

double PDESolver::interpCB(vector<double> x,vector<double> y,double ind)
{
    //https://zhidao.baidu.com/question/619575742421839812.html MATLAB插值函数原理
    if(x.size()!=y.size()) 
    {
        cout<<"The two list have different length!"<<endl;
        return -1;
    }

    int i=0;
    while (i<ind) i++;
    double ratio = (ind - x[i])/(x[i+1] - x[i]);
    double res = y[i] + (y[i+1] - y[i])*ratio;

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

    double t; //循环中用到的全局变量 定义在循环外
    double rate1;
    double rate2;
    double r;
    vector<double> alpha;
    vector<double> beta;
    vector<double> alpha_MuMB;
    vector<double> beta_MuMB;
    vector<double> gamma_Mu;

    //初始化一个二维动态数组 方法参考 https://blog.csdn.net/samuelcoulee/article/details/8674388
    //用于存储MATLAB源代码中的Mu，其是一个稀疏矩阵 
    int height = int(this->pde_param_ptr->Ns)+1;  //height = 2401
    int width = int(this->pde_param_ptr->Ns)+1;  //width=2401
    double** Mu = new double *[height];
    for(int i=0;i<height;i++) Mu[i] = new double[width];
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            Mu[i][j] = 0;
        }
    }
    
    vector<double> gamma_MB;

    // int height = int(this->pde_param_ptr->Ns)+1;  //height = 2401
    // int width = int(this->pde_param_ptr->Ns)+1;  //width=2401
    double** MB = new double *[height];
    for(int i=0;i<height;i++) MB[i] = new double[width];
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            MB[i][j] = 0;
        }
    }    

    double AccI;
    double Bc_n;
    double Bp_n;
    vector<double> k_n;
    vector<double> u_n;
    vector<double> Muu;
    // for(int n=0;n<this->pde_param_ptr->Nt;n++)
    for(int n=1;n<2;n++)
    {
        if(flag==1) break;
        // double t = n*this->pde_param_ptr->dt;
        t = n*this->pde_param_ptr->dt;
        cout<<"t = "<<t<<endl;
        // double rate1 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        rate1 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        cout<<"rate1 =  "<<rate1<<endl;
        // double rate2 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        rate2 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        cout<<"rate2 = "<<rate2<<endl;
        // double r = -log((rate2/rate1))/this->pde_param_ptr->dt;
        r = -log((rate2/rate1))/this->pde_param_ptr->dt;
        cout<<"r = "<<r<<endl;

        if ((ceil(t-this->pde_param_ptr->dtime) - ceil(t-this->pde_param_ptr->dt-this->pde_param_ptr->dtime)) == 1)
        {//此段代码对应MATLAB源码中的84至88行
            count6++;
            vector<double> tmpB = this->pde_param_ptr->B;
            tmpB.clear();
            for(int i=0;i<this->pde_param_ptr->B.size();i++) tmpB.push_back(this->interpCB(this->pde_param_ptr->S,this->pde_param_ptr->B,(1-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]));
            this->pde_param_ptr->B = tmpB;

            vector<double>tmpu = this->pde_param_ptr->u;
            tmpu.clear();
            for(int i=0;i<this->pde_param_ptr->u.size();i++) tmpu.push_back(this->interpCB(this->pde_param_ptr->S,this->pde_param_ptr->u,(1-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]));
            this->pde_param_ptr->u = tmpu;
        }

        alpha.clear();
        for(int i=1;i<this->pde_param_ptr->S.size()-1;i++)//循环赋值给alpha 在MATLAB源代码92行中是以向量化操作实现的
        { //alpha shape 2399 x 1
            double tmp1 = (this->pde_param_ptr->sigma*this->pde_param_ptr->sigma*this->pde_param_ptr->S[i]*this->pde_param_ptr->S[i])/(2*this->pde_param_ptr->ds*this->pde_param_ptr->ds);
            double tmp2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]/(2*this->pde_param_ptr->ds);
            double tmp3 = this->pde_param_ptr->dt;
            alpha.push_back((tmp1-tmp2)*tmp3);
        }
        cout<<"The first 4 and last 2 values of alpha = "<<alpha[0]<<" "<<alpha[1]<<" "<<alpha[2]<<" "<<alpha[3]<<" "<<alpha[2397]<<" "<<alpha[2398]<<endl;
        
        beta.clear();
        for(int i=1;i<this->pde_param_ptr->S.size()-1;i++)//循环赋值给alpha 在MATLAB源代码92行中是以向量化操作实现的
        { //alpha shape 2399 x 1
            double tmp1 = (this->pde_param_ptr->sigma*this->pde_param_ptr->sigma*this->pde_param_ptr->S[i]*this->pde_param_ptr->S[i])/(2*this->pde_param_ptr->ds*this->pde_param_ptr->ds);
            double tmp2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*this->pde_param_ptr->S[i]/(2*this->pde_param_ptr->ds);
            double tmp3 = this->pde_param_ptr->dt;
            beta.push_back((tmp1+tmp2)*tmp3);
        }
        cout<<"The first 4 and last 2 values of beta = "<<beta[0]<<" "<<beta[1]<<" "<<beta[2]<<" "<<beta[3]<<" "<<beta[2397]<<" "<<beta[2398]<<endl;
    
        alpha_MuMB.clear();
        alpha_MuMB = alpha;
        alpha_MuMB.push_back(0);
        alpha_MuMB.push_back(0);
        cout<<"The length of alpha_MuMB is:"<<alpha_MuMB.size()<<endl;
        cout<<"The first 2 and last 2 values of alpha_MuMB = "<<alpha_MuMB[0]<<" "<<alpha_MuMB[1]<<" "<<alpha_MuMB[2399]<<" "<<alpha_MuMB[2400]<<endl;

        beta_MuMB.clear();
        beta_MuMB.push_back(0);
        beta_MuMB.push_back(0);
        beta_MuMB.insert(beta_MuMB.end(),beta.begin(),beta.end());
        cout<<"The length of beta_MuMB is:"<<beta_MuMB.size()<<endl;
        cout<<"The first 2 and last 2 values of beta_MuMB = "<<beta_MuMB[0]<<" "<<beta_MuMB[1]<<" "<<beta_MuMB[2399]<<" "<<beta_MuMB[2400]<<endl;

        gamma_Mu.clear();
        gamma_Mu.push_back(-(r+this->pde_param_ptr->p)*this->pde_param_ptr->dt);
        for(int i=0;i<alpha.size();i++)
        {
            double tmp = -(alpha[i] + beta[i] + (r+this->pde_param_ptr->p)*this->pde_param_ptr->dt);
            gamma_Mu.push_back(tmp);
        }
        gamma_Mu.push_back(0);
        cout<<"The length of gamma_Mu is:"<<gamma_Mu.size()<<endl;
        cout<<"The first 3 and last 3 values of gamma_Mu = "<<gamma_Mu[0]<<" "<<gamma_Mu[1]<<" "<<gamma_Mu[2]<<" "<<gamma_Mu[2399]<<" "<<gamma_Mu[2400]<<endl;
    
        //此段代码用于给Mu赋值，源代码102行创建稀疏矩阵Mu，此处通过循环先给主对角线赋值
        //上副对角线的值是beta_MuMB
        //下副对角线的值是alpha_MuMB
        for(int i=0;i<height;i++)
        {
            Mu[i][i] = gamma_Mu[i];
        }
        for(int i=1;i<height;i++)
        {
            Mu[i-1][i] = beta_MuMB[i];
        }
        for(int i=1;i<height;i++)
        {
            Mu[i][i-1] = alpha_MuMB[i-1];
        }
        cout<<"row 1 col 1 of Mu:"<<Mu[0][0]<<endl;
        cout<<"row 2 col 2 of Mu:"<<Mu[1][1]<<endl;
        cout<<"row 1 col 2 of Mu:"<<Mu[0][1]<<endl;
        cout<<"row 2 col 3 of Mu:"<<Mu[1][2]<<endl;
        cout<<"row 3 col 4 of Mu:"<<Mu[2][3]<<endl;
        cout<<"row 2 col 1 of Mu:"<<Mu[1][0]<<endl;
        cout<<"row 3 col 2 of Mu:"<<Mu[2][1]<<endl;
        cout<<"row 2401 col 2401 of Mu:"<<Mu[2400][2400]<<endl;
        cout<<"row 2400 col 2401 of Mu:"<<Mu[2399][2400]<<endl;
        cout<<"row 2401 col 2400 of Mu:"<<Mu[2400][2399]<<endl;

        gamma_MB.clear();
        gamma_MB.push_back(-(r+this->pde_param_ptr->p*(1-this->pde_param_ptr->R))*this->pde_param_ptr->dt);
        for(int i=0;i<alpha.size();i++)
        {
            double tmp = -(alpha[i]+beta[i]+  (r+this->pde_param_ptr->p*(1-this->pde_param_ptr->R))*this->pde_param_ptr->dt);
            gamma_MB.push_back(tmp);
        }  
        gamma_MB.push_back(0);

        cout<<"The length of gamma_MB is:"<<gamma_MB.size()<<endl;
        cout<<"The first 2 and last 2 values of gamma_MB = "<<gamma_MB[0]<<" "<<gamma_MB[1]<<" "<<gamma_MB[2399]<<" "<<gamma_MB[2400]<<endl;
        cout<<"$$$"<<gamma_MB[4]<<endl;

        for(int i=0;i<height;i++)
        {
            MB[i][i] = gamma_MB[i];
        }
        for(int i=1;i<height;i++)
        {
            MB[i-1][i] = beta_MuMB[i];
        }
        for(int i=1;i<height;i++)
        {
            MB[i][i-1] = alpha_MuMB[i-1];
        }
        cout<<"row 1 col 1 of MB:"<<MB[0][0]<<endl;
        cout<<"row 2 col 2 of MB:"<<MB[1][1]<<endl;
        cout<<"row 1 col 2 of MB:"<<MB[0][1]<<endl;
        cout<<"row 2 col 3 of MB:"<<MB[1][2]<<endl;
        cout<<"row 3 col 4 of MB:"<<MB[2][3]<<endl;
        cout<<"row 2 col 1 of MB:"<<MB[1][0]<<endl;
        cout<<"row 3 col 2 of MB:"<<MB[2][1]<<endl;
        cout<<"row 2401 col 2401 of MB:"<<MB[2400][2400]<<endl;
        cout<<"row 2400 col 2401 of MB:"<<MB[2399][2400]<<endl;
        cout<<"row 2401 col 2400 of MB:"<<MB[2400][2399]<<endl;
        //至此，Mu，MB两个稀疏矩阵已经构建好了，并且数据与MATLAB核对无误
        AccI = (ceil(t) - t)*this->pde_param_ptr->coupon_time_rate[ceil(t)]*this->pde_param_ptr->F;
        cout<<"AccI = "<<AccI<<endl;

        if(this->pde_param_ptr->no_call_time[0]<=t&&t<=this->pde_param_ptr->no_call_time[1])
            Bc_n = 1000;
        else 
            Bc_n=this->pde_param_ptr->Bc_star+AccI;
        cout<<"Bc_n = "<<Bc_n<<endl;
        if(this->pde_param_ptr->no_put_time[0]<=t&&t<=this->pde_param_ptr->no_put_time[1])
            Bp_n=0.01;
        else
            Bp_n=this->pde_param_ptr->Bp_star+AccI;
        cout<<"Bp_n = "<<Bp_n<<endl;

        if(this->pde_param_ptr->no_convert_time[0]<=t&&t<=this->pde_param_ptr->no_convert_time[1])
        {
            for(int i=0;i<this->pde_param_ptr->S.size();i++) 
            {   
                 k_n.push_back(0.01*1); 
            }
        }
        else
        {
            for(int i=0;i<this->pde_param_ptr->S.size();i++)
            {
                double tmp1 = this->pde_param_ptr->F;
                double tmp2 = this->pde_param_ptr->conversion_price-this->pde_param_ptr->q*this->pde_param_ptr->S[i]*(ceil(this->pde_param_ptr->T-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dtime));
                k_n.push_back(tmp1/tmp2);
            }
        }
        cout<<"The length of k_n is:"<<k_n.size()<<endl;
        cout<<"The first 2 value and last 2 value of k_n ="<<k_n[0]<<" "<<k_n[1]<<" "<<k_n[2399]<<" "<<k_n[2400]<<endl;

        // 将u的值拷贝一份 给u_n  方法参考 https://blog.csdn.net/lengyuezuixue/article/details/78480040
        u_n.assign(this->pde_param_ptr->u.begin(),this->pde_param_ptr->u.end());
        cout<<"The length of u_n is:"<<u_n.size()<<endl;
        cout<<"The first 2 value and last 2 value of u_n ="<<u_n[0]<<" "<<u_n[1]<<" "<<u_n[2399]<<" "<<u_n[2400]<<endl;


        cout<<"The length of B ="<<this->pde_param_ptr->B.size()<<endl;
        cout<<"The first 2 value and last 2 value of B ="<<this->pde_param_ptr->B[0]<<" "<<this->pde_param_ptr->B[1]<<" "<<this->pde_param_ptr->B[2399]<<" "<<this->pde_param_ptr->B[2400]<<endl;
        cout<<"***************"<<endl;
        for(int i=0;i<this->pde_param_ptr->B.size();i++)
        {
            double tmp1 = 1 - this->pde_param_ptr->theta*MB[i][i];
            double tmp2 = (1 + (1 - this->pde_param_ptr->theta)*MB[i][i])*this->pde_param_ptr->B[i];
            if(i==0) cout<<"$$$ "<<tmp1/tmp2<<endl;
            this->pde_param_ptr->B[i] = tmp2/tmp1;
            
        }
        cout<<"The length of B ="<<this->pde_param_ptr->B.size()<<endl;
        cout<<"The first 2 value and last 2 value of B ="<<this->pde_param_ptr->B[0]<<" "<<this->pde_param_ptr->B[1]<<" "<<this->pde_param_ptr->B[2399]<<" "<<this->pde_param_ptr->B[2400]<<endl;

        for(int i=0;i<this->pde_param_ptr->B.size();i++)
        {
            if(this->pde_param_ptr->B[i]>Bc_n) this->pde_param_ptr->B[i] = Bc_n;
        }


        for(int i=0;i<this->pde_param_ptr->S.size();i++)
        {
            double tmp1 = max(k_n[i]*(1-this->pde_param_ptr->eta)*this->pde_param_ptr->S[i],this->pde_param_ptr->R*this->pde_param_ptr->B[i]);
            double tmp2 = this->pde_param_ptr->p*this->pde_param_ptr->dt;
            Muu.push_back(tmp1*tmp2);
        }
        cout<<"The length of Muu ="<<Muu.size()<<endl;
        cout<<"The first 2 value and last 2 value of Muu ="<<Muu[0]<<" "<<Muu[1]<<" "<<Muu[2399]<<" "<<Muu[2400]<<endl;

        cout<<"The length of u ="<<this->pde_param_ptr->u.size()<<endl;
        cout<<"The first 2 value and last 2 value of u ="<<this->pde_param_ptr->u[0]<<" "<<this->pde_param_ptr->u[1]<<" "<<this->pde_param_ptr->u[2399]<<" "<<this->pde_param_ptr->u[2400]<<endl;
        cout<<"$$$$$$$$$$$$$$$$$"<<endl;
        for(int i=0;i<u_n.size();i++)
        {
            double tmp1 = 1 - this->pde_param_ptr->theta*Mu[i][i];
            double tmp2 = (1 + (1-this->pde_param_ptr->theta)*Mu[i][i])*u_n[i];
            double tmp3;
            if (i <= this->pde_param_ptr->Ns-1)
                tmp3 = Muu[i];
            else
                tmp3 = 0;

            this->pde_param_ptr->u[i] = (tmp2 + tmp3)/tmp1;
            if(i<4) cout<<"tmp1 tmp2 tmp3 u_n[i] Mu[i] " <<tmp1<<" "<<tmp2<<" "<<tmp3<<" "<<u_n[i]<<" "<<Mu[i][i]<<endl;
            if(i==2398) cout<<"tmp1 tmp2 tmp3 u_n[i] Mu[i] " <<tmp1<<" "<<tmp2<<" "<<tmp3<<" "<<u_n[i]<<" "<<Mu[i][i]<<endl;
            
        }
        cout<<"The length of u ="<<this->pde_param_ptr->u.size()<<endl;
        cout<<"The first 3 value and last 3 value of u ="<<this->pde_param_ptr->u[0]<<" "<<this->pde_param_ptr->u[1]<<" "<<this->pde_param_ptr->u[2]<<" "<<this->pde_param_ptr->u[2398]<<" "<<this->pde_param_ptr->u[2399]<<" "<<this->pde_param_ptr->u[2400]<<endl;
        cout<<"u[4]"<<this->pde_param_ptr->u[4]<<endl;
        cout<<"u[5]"<<this->pde_param_ptr->u[5]<<endl;
    }

    return 1;
}