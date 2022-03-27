#include<iostream>
#include<vector>
#include"Parameter.h"
#include"PDESolver.h"
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

int main()
{


    double T = 0.5;
    double Ns = 2400;
    double F = 100;
    double conversion_price = 42.79;
    double  S0 = 52.75;
    double sigma = 0.3034;
    double final_coupon_rate = 0.015;
    double q = 0;
    double p =  0;
    double R = 1;
    double eta = 0;

    MatrixXd no_call_time(2,1);
    no_call_time<<0,3;
    MatrixXd no_put_time(2,1);
    no_put_time<<0,3;
    MatrixXd no_convert_time(2,1);
    no_convert_time<<2,3;

    double Bc_star = 100000;
    double Bp_star =  0.01;
    double theta = 1;
    
    MatrixXd coupon_time_rate(12,1);
    coupon_time_rate<<0.015, 0.015, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    // cout<<"coupon time rate \n"<<coupon_time_rate<<endl;
    MatrixXd risk_free_rate(6,1);
    risk_free_rate<<1, 0.999691000000000, 0.999233000000000, 0.99081900000000, 0.982017000000000, 0.973488000000000;
    // cout<<"risk free rate \n"<<risk_free_rate<<endl;
    MatrixXd rate_T(6,1);
    rate_T<<0.500000000000000, 0.491666666666667, 0.480555555555556, 0.244444444444444, -0.00555555555555556, -0.250000000000000;
    cout<<"rate_T[1,0] = "<<rate_T(1,0)<<endl;
    double dtime= 0.33;


    Parameter paras(T,Ns,F,conversion_price,S0,sigma,final_coupon_rate,\
    q,p,R,eta,Bc_star,Bp_star,theta,dtime,no_call_time,no_put_time,no_convert_time,coupon_time_rate,risk_free_rate,rate_T);
    
    PDESolver pdesolver(&paras,9);
    pdesolver.set_iter(9);
    pdesolver.show_iter();

    pdesolver.pde_param_ptr->show();
    cout<<"\n*********show complete*********"<<endl;

    double res = pdesolver.solve();
    cout<<"Final result = "<<res<<endl;


    return 0;
}