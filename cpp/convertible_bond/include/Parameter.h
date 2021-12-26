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

private:
    int version;
    double T;
    double Nt; //calculated
    double Ns;
    double F;
    double conversion_price;
    double S0;
    double Smax; //calculated
    double sigma;
    double final_coupon_rate;
    double q;
    double p;
    double R;
    double eta;
    vector<int> no_call_time;
    vector<int> no_put_time;
    vector<int> no_convert_time;
    double Bc_star;
    double Bp_star;
    double theta;
    vector<double> coupon_time_rate;
    vector<double> risk_free_rate;
    vector<double> rate_T;
    double dtime;
    double k; //calculated
    double ds; //calculated
    double dt; //calculated
    double rhopenltycall; //calculated
    double rhopenltyput; //calculated

    vector<double> i; //boundary condition
    vector<double> S; //boundary condition
    double Bc_T; //boundary condition
    double Bp_T; //boundary condition
    vector<double> k_T; //boundary condition
    vector<double> u; //boundary condition
    vector<double> B; //boundary condition


};  

#endif