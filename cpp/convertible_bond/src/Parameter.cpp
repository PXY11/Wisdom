#include <Parameter.h>
#include <iostream>
#include <string>
#include <json/json.h>
#include <iostream>
#include <fstream>
#include <vector>
#include<math.h>
using namespace std;

void Parameter::show()
{
    cout<<"The value of T is:"<<this->T<<endl;
    cout<<"The value of Nt is:"<<this->Nt<<endl;
    cout<<"The value of Ns is:"<<this->Ns<<endl;
    cout<<"The value of F is:"<<this->F<<endl;
    cout<<"The value of conversion_price is:"<<this->conversion_price<<endl;
    cout<<"The value of S0 is:"<<this->S0<<endl;
    cout<<"The value of Smax is:"<<this->Smax<<endl;
    cout<<"The value of sigma is:"<<this->sigma<<endl;
    cout<<"The value of final_coupon_rate is:"<<this->final_coupon_rate<<endl;
    cout<<"The value of q is:"<<this->q<<endl;
    cout<<"The value of p is:"<<this->p<<endl;
    cout<<"The value of R is:"<<this->R<<endl;
    cout<<"The value of eta is:"<<this->eta<<endl;
    cout<<"The length of no_call_time is:"<<this->no_call_time.size()<<endl;
    cout<<"The length of no_put_time is:"<<this->no_put_time.size()<<endl;
    cout<<"The length of no_convert_time is:"<<this->no_convert_time.size()<<endl;
    cout<<"The value of Bc_star is:"<<this->Bc_star<<endl;
    cout<<"The value of Bp_star is:"<<this->Bp_star<<endl;
    cout<<"The value of theta is:"<<this->theta<<endl;
    cout<<"The value of dtime is:"<<this->dtime<<endl;
    cout<<"The value of k is:"<<this->k<<endl;
    cout<<"The value of ds is:"<<this->ds<<endl;
    cout<<"The value of dt is:"<<this->dt<<endl;
    cout<<"The value of rhopenltycall is:"<<this->rhopenltycall<<endl;
    cout<<"The value of rhopenltyput is:"<<this->rhopenltyput<<endl;
    cout<<"The length of coupon_time_rate is:"<<this->coupon_time_rate.size()<<endl;
    cout<<"The length of risk_free_rate is:"<<this->risk_free_rate.size()<<endl;
    cout<<"The length of rate_T is:"<<this->rate_T.size()<<endl;
    cout<<"The length of i is:"<<this->i.size()<<endl;
}

void Parameter::readParam()
{
    Json::Reader reader;
    Json::Value root;

    ifstream in("param_setting/CBparam.json",ios::binary);

    if(!in.is_open()){
        cout<<"Error opening"<<endl;
        return;
    }

    if(reader.parse(in,root)){

        this->T = root["T"].asDouble();
        this->Ns = root["Ns"].asDouble();
        this->F = root["F"].asDouble();
        this->conversion_price = root["conversion_price"].asDouble();
        this->S0 = root["S0"].asDouble();
        this->Smax = root["Smax"].asDouble();
        this->sigma = root["sigma"].asDouble();
        this->final_coupon_rate = root["final_coupon_rate"].asDouble();
        this->q = root["q"].asDouble();
        this->p = root["p"].asDouble();
        this->R = root["R"].asDouble();
        this->eta = root["eta"].asDouble();
        this->dtime = root["dtime"].asDouble();
        this->Bc_star = root["Bc_star"].asDouble();
        this->Bp_star = root["Bp_star"].asDouble();
        this->theta = root["theta"].asDouble();
        const Json::Value no_call_time = root["no_call_time"];
        for (int i=0;i<no_call_time.size();i++)
        {
            this->no_call_time.push_back(no_call_time[i].asInt());
        }

        const Json::Value no_put_time = root["no_put_time"];
        for (int i=0;i<no_put_time.size();i++)
        {
            this->no_put_time.push_back(no_put_time[i].asInt());
        }

        const Json::Value no_convert_time = root["no_convert_time"];
        for (int i=0;i<no_convert_time.size();i++)
        {
            this->no_convert_time.push_back(no_convert_time[i].asInt());
        }
        
        const Json::Value coupon_time_rate = root["coupon_time_rate"];
        for (int i=0;i<coupon_time_rate.size();i++)
        {
            this->coupon_time_rate.push_back(coupon_time_rate[i].asDouble());
        }

        const Json::Value risk_free_rate = root["risk_free_rate"];
        for (int i=0;i<risk_free_rate.size();i++)
        {
            this->risk_free_rate.push_back(risk_free_rate[i].asDouble());
        }

        const Json::Value rate_T = root["rate_T"];
        //cout<<typeid(rate_T).name()<<endl;
        for (int i=0;i<rate_T.size();i++)
        {
            this->rate_T.push_back(rate_T[i].asDouble());
        }

    }

    in.close();
}

void Parameter::calParam()
{
    this->Nt = ceil(500*this->T);
    this->Smax = 8*this->S0;
    this->k = this->F/this->conversion_price;
    this->ds = this->Smax/this->Ns;
    this->dt = this->T/this->Nt;
    this->rhopenltycall = 1000000/(this->dt*this->dt);
    this->rhopenltyput = 1000000/(this->dt*this->dt);
}

void Parameter::setBoundaryParam()
{
    //设置i
    for(int i=1;i<this->Ns+2;i++)
    {
        this->i.push_back(i);
    }
    
    //设置S
    for(int i=0;i<this->Ns+1;i++)
    {
        double tmp = i*this->ds;
        this->S.push_back(tmp);
    }
    
    this->Bc_T = this->Bc_star + this->final_coupon_rate*this->F;
    this->Bp_T = this->Bp_star + this->final_coupon_rate*this->F;

    //设置k_T
    for (int i=0;i<this->S.size();i++)
    {
        double tmp = this->F/(this->conversion_price-this->q*ceil(this->T-this->dtime));
        this->k_T.push_back(tmp);
    }
    
    //设置u
    for (int i =0;i<this->S.size();i++)
    {
        double tmp_1 = this->k_T[i]*this->S[i]+this->coupon_time_rate[0]*this->F;
        double tmp_2 = this->F+this->final_coupon_rate*this->F+this->coupon_time_rate[0]*this->F;
        double tmp_3 = max(tmp_1,tmp_2);
        u.push_back(tmp_3);
        //cout<<i+1<<" "<<tmp_3<<endl;
    }

    //设置B
    for(int i=0;i<this->Ns+1;i++)
    {
        double tmp = this->F + this->final_coupon_rate*this->F;
        B.push_back(tmp);
    }

    cout<<"***"<<this->u[0]<<"***"<<endl;
    cout<<"***"<<this->u[this->S.size()-1]<<"***";

}