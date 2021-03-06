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

    //mac的编译器版本问题，不能直接用列表赋值，用新的赋值方法 https://blog.csdn.net/qwer7512090/article/details/105272667
    // vector<int> no_call_time = { 0, 3 };
    // vector<int> no_put_time = { 0, 3 };
    // vector<int> no_convert_time = { 2, 3 };
    // int nocalltime[2] = {0,3};
    // vector<int> no_call_time(nocalltime,nocalltime+2);
    // int noputtime[2] = {0,3};
    // vector<int> no_put_time(noputtime,noputtime+2);
    // int noconvertime[2] = {2,3};
    // vector<int> no_convert_time(noconvertime,noconvertime+2);

    MatrixXd no_call_time(2,1);
    no_call_time<<0,3;
    MatrixXd no_put_time(2,1);
    no_put_time<<0,3;
    MatrixXd no_convert_time(2,1);
    no_convert_time<<2,3;

    double Bc_star = 100000;
    double Bp_star =  0.01;
    double theta = 1;

    //mac的编译器版本问题，不能直接用列表赋值，用新的赋值方法 https://blog.csdn.net/qwer7512090/article/details/105272667
    // vector<double> coupon_time_rate = { 0.015, 0.015, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // vector<double> risk_free_rate = { 1, 0.999691000000000, 0.999233000000000, 0.99081900000000, 0.982017000000000, 0.973488000000000 };
    // vector<double> rate_T = { 0.500000000000000, 0.491666666666667, 0.480555555555556, 0.244444444444444, -0.00555555555555556, -0.250000000000000 };
    // double coupontimerate[12] = { 0.015, 0.015, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // vector<double> coupon_time_rate(coupontimerate,coupontimerate+12);
    // double riskfreerate[6] = { 1, 0.999691000000000, 0.999233000000000, 0.99081900000000, 0.982017000000000, 0.973488000000000 };
    // vector<double> risk_free_rate(riskfreerate,riskfreerate+6);
    // double rateT[6] = { 0.500000000000000, 0.491666666666667, 0.480555555555556, 0.244444444444444, -0.00555555555555556, -0.250000000000000 };
    // vector<double> rate_T(rateT,rateT+6);
    
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



    Parameter paras(T,Ns,F,conversion_price,S0,sigma,final_coupon_rate,q,p,R,eta,Bc_star,Bp_star,theta,dtime,no_call_time,no_put_time,no_convert_time,coupon_time_rate,risk_free_rate,rate_T);
    
    PDESolver pdesolver(&paras,9);
    pdesolver.set_iter(9);
    pdesolver.show_iter();

    // cout<<"*******First value of B"<<endl;
    // cout<<paras.B[0]<<endl;  
    // cout<<pdesolver.pde_param_ptr->B[0]<<endl;
    // cout<<"***********************"<<endl;
    pdesolver.pde_param_ptr->show();
    cout<<"\n*********show complete*********"<<endl;


    /* 测试convertible_code.m 功能*/
    // pdesolver.pde_param_ptr->setBoundaryParam(); 
    double res = pdesolver.solve();
    cout<<"Final res = "<<res<<endl;


    //测试矩阵求逆功能
    // double matsrc[2][2] = {{1,0},{3,4}};
    // double matdes[2][2] = {{0,0},{0,0}};
    // bool ft = pdesolver.GetMatrixInverse(matsrc,2,matdes);
    // cout<<"ft = "<<ft<<endl;

    // double left[2][2] = {{1,1},{2,3}};
    // double right[2][2] = {{3,4},{5,6}};
    // pdesolver.MatMul(left,right);

    
    // pdesolver.pde_param_ptr->show();  //用于输出pde parameters
    // cout<<"rhopenltycall="<<pdesolver.pde_param_ptr->rhopenltycall<<endl;
    //cout<<pdesolver.pde_param_ptr->dt<<endl;


    return 0;
}