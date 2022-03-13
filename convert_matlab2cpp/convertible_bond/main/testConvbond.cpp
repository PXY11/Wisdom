#include<iostream>
#include<vector>
#include<Parameter.h>
// #include<PDESolver.h>
//#include<dataStructure.h>
using namespace std;

int main(){


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

    //通过数组地址初始化
    // int a[5]={1,2,3,4,5}
    // vector<int> vec_i(a,a+5);
    int nocalltime[2] = {0,3};
    // vector<int> no_call_time = { 0, 3};
    vector<int> no_call_time(nocalltime,nocalltime+2);
    // for (int i=0;i<no_call_time.size();i++)
    // {
    //     cout<<no_call_time[i]<<endl;
    // }
    int noputtime[2] = {0,3};
    // vector<int> no_put_time  { 0, 3 };
    vector<int> no_put_time(noputtime,noputtime+2);
    // vector<int> no_convert_time { 2, 3 };
    int noconvertime[2] = {2,3};
    vector<int> no_convert_time(noconvertime,noconvertime+2);

    double Bc_star = 100000;
    double Bp_star =  0.01;
    double theta = 1;
    // vector<double> coupon_time_rate  { 0.015, 0.015, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // vector<double> risk_free_rate  { 1, 0.999691000000000, 0.999233000000000, 0.99081900000000, 0.982017000000000, 0.973488000000000 };
    // vector<double> rate_T  { 0.500000000000000, 0.491666666666667, 0.480555555555556, 0.244444444444444, -0.00555555555555556, -0.250000000000000 };
    double coupontimerate[12] = { 0.015, 0.015, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<double> coupon_time_rate(coupontimerate,coupontimerate+12);
    double riskfreerate[6] = { 1, 0.999691000000000, 0.999233000000000, 0.99081900000000, 0.982017000000000, 0.973488000000000 };
    vector<double> risk_free_rate(riskfreerate,riskfreerate+6);
    double rateT[6] = { 0.500000000000000, 0.491666666666667, 0.480555555555556, 0.244444444444444, -0.00555555555555556, -0.250000000000000 };
    vector<double> rate_T(rateT,rateT+6);

    double dtime= 0.33;

    cout<<no_call_time[0]<<" "<<no_put_time[0]<<" "<<no_convert_time[0]<<" "<<coupon_time_rate[0]<<" "<<risk_free_rate[0]<<" "<<rate_T[0];
    // PDESolver pdesolver(1);
    // pdesolver.set_iter(9);
    // pdesolver.show_iter();
    Parameter* param_ptr;
    // param_ptr = pdesolver.get_pde_param_ptr();

    /* 测试convertible_code.m 功能*/
    // pdesolver.pde_param_ptr->readParam();
    // pdesolver.pde_param_ptr->calParam();
    // pdesolver.pde_param_ptr->setBoundaryParam();
    // double res = pdesolver.solve();
    // cout<<"Final res = "<<res;
    // pdesolver.pde_param_ptr->show();  //用于输出pde parameters
    // cout<<"rhopenltycall="<<pdesolver.pde_param_ptr->rhopenltycall<<endl;
    //cout<<pdesolver.pde_param_ptr->dt<<endl;
    

    return 0;
}