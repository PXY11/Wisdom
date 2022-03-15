#ifndef PARAMETER_H
#define PARAMETER_H
#include <string>
#include <vector>
#include<cmath>
using namespace std;
class Parameter
{

public:
    Parameter(int ver){this->version=ver;}
    Parameter(double _T, double _Ns, double _F, double _conversion_price, double _S0,  double _sigma, double _final_coupon_rate, double _q, double _p, double _R, double _eta,
              double _Bc_star, double _Bp_star, double _theta, double _dtime,
              vector<int> _no_call_time, vector<int> _no_put_time, vector<int> _no_convert_time, vector<double> _coupon_time_rate, vector<double> _risk_free_rate, vector<double> _rate_T)
    {
        T = _T;
        Ns = _Ns;
        F = _F;
        conversion_price = _conversion_price;
        S0 = _S0;
        sigma = _sigma;
        final_coupon_rate = _final_coupon_rate;
        q = _q;
        p = _p;
        R = _R;
        eta = _eta;
        Bc_star = _Bc_star;
        Bp_star = _Bp_star;
        theta = _theta;
        dtime = _dtime;
        no_call_time = _no_call_time;
        no_put_time = _no_put_time;
        no_convert_time = _no_convert_time;
        coupon_time_rate = _coupon_time_rate;
        risk_free_rate = _risk_free_rate; //risk_free_rate在matlab文件35行之后 变成了逆序存储
        risk_free_rate.clear();
        for(int i=_risk_free_rate.size()-1;i>=0;i--)
        {risk_free_rate.push_back(_risk_free_rate[i]);}
        rate_T = _rate_T; //rate_T在matlab文件 34行 [rate_T,index]=unique(rate_T); 之后，变成了逆序存储
        rate_T.clear(); // 转换rate_T为逆序存储
        for(int i=_rate_T.size()-1;i>=0;i--) // 转换rate_T为逆序存储
        {rate_T.push_back(_rate_T[i]);} // 转换rate_T为逆序存储

        //---calParam
        Nt = ceil(500 * T);
        Smax = 8 * S0;
        k = F / conversion_price;
        ds = Smax / Ns;
        dt = T / Nt;
        rhopenltycall = 1000000 / (dt *dt);
        rhopenltyput = 1000000 / (dt * dt);

        //boundary condition
            //����i
        for (int i = 1; i < this->Ns + 2; i++)
        {
            this->i.push_back(i);
        }

        //����S
        for (int i = 0; i < this->Ns + 1; i++)
        {
            double tmp = i * this->ds;
            this->S.push_back(tmp);
        }

        this->Bc_T = this->Bc_star + this->final_coupon_rate * this->F;
        this->Bp_T = this->Bp_star + this->final_coupon_rate * this->F;

        //����k_T
        for (int i = 0; i < this->S.size(); i++)
        {
            double tmp = this->F / (this->conversion_price - this->q * ceil(this->T - this->dtime));
            this->k_T.push_back(tmp);
        }

        //����u
        for (int i = 0; i < this->S.size(); i++)
        {
            double tmp_1 = this->k_T[i] * this->S[i] + this->coupon_time_rate[0] * this->F;
            double tmp_2 = this->F + this->final_coupon_rate * this->F + this->coupon_time_rate[0] * this->F;
            double tmp_3 = max(tmp_1, tmp_2);
            u.push_back(tmp_3);
            //cout<<i+1<<" "<<tmp_3<<endl;
        }

        //����B
        for (int i = 0; i < this->Ns + 1; i++)
        {
            double tmp = this->F + this->final_coupon_rate * this->F;
            B.push_back(tmp);
        }
    }
    void show();
    //void readParam();
    //void calParam();
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