#include"PDESolver.h"
#include<iostream>
#include<cmath>
#include <Eigen/Dense>
 
using namespace Eigen;
using namespace std;


void PDESolver::set_iter(int iter_time)
{
    this->iteration = iter_time;
}

void PDESolver::show_iter()
{
    cout<<"iter time = "<<this->iteration<<endl;
}
// double PDESolver::interpCB(vector<double> x,vector<double> y,double ind)
// {
//     //https://zhidao.baidu.com/question/619575742421839812.html MATLAB插值函数原理
//     if(x.size()!=y.size()) 
//     {
//         cout<<"The two list have different length!"<<endl;
//         return -1;
//     }

//     int i=0;
//     while (i<ind) i++;
//     double ratio = (ind - x[i])/(x[i+1] - x[i]);
//     double res = y[i] + (y[i+1] - y[i])*ratio;
//     return res;
// }
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
    double t;
    double rate1;
    double rate2;
    double r;
    for(int n=1;n<2;n++)
    {
        if(flag==1) break;
        t = n*this->pde_param_ptr->dt;
        cout<<"t = "<<t<<endl;
        // rate1 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        // cout<<"rate1 =  "<<rate1<<endl;
        // rate2 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        // cout<<"rate2 = "<<rate2<<endl;
        // r = -log((rate2/rate1))/this->pde_param_ptr->dt;
        // cout<<"r = "<<r<<endl;
    }

    return 1;
}