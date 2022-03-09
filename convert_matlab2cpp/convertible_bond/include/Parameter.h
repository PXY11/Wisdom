#ifndef PARAMETER_H
#define PARAMETER_H
#include <string>
#include <vector>
using namespace std;
class Parameter
{

public:
    Parameter(int ver){this->version=ver;}
    void show();
    void readParam();
    void calParam();
    void setBoundaryParam();
    ~Parameter(){}

//private:
    int version;
    double T; // time to maturity
    double Nt; //calculated time fraction
    double Ns; // price fraction
    double F; // face value
    double conversion_price; // convert price
    double S0; // present stock price
    double Smax; //calculated upbound of stock price
    double sigma; //volatility
    double final_coupon_rate; //coupon rate in the final time
    double q; // dividend
    double p; // default probability
    double R; // recovery rate
    double eta; // when default te stock price becomes (1*eta)*S
    vector<int> no_call_time; // no allow call time
    vector<int> no_put_time; // no allow put time
    vector<int> no_convert_time; // no allow convert time
    double Bc_star; // callable price
    double Bp_star; // puttable price
    double theta; // implicititness parameter 
    vector<double> coupon_time_rate; //coupon rate
    vector<double> risk_free_rate; 
    vector<double> rate_T;
    double dtime; // divided at time t=X.33
    double k; //calculated conversion ratio
    double ds; //calculated price direction
    double dt; //calculated time direction 
    double rhopenltycall; //calculated
    double rhopenltyput; //calculated

    vector<double> i; //boundary condition stock price [1,2...,Ns+1]
    vector<double> S; //boundary condition stock price [0,h,2h...Ns*h]
    double Bc_T; //boundary condition Bc in te final time using dirty price
    double Bp_T; //boundary condition 
    vector<double> k_T; //boundary condition
    vector<double> u; //boundary condition 
    vector<double> B; //boundary condition


};  

#endif