#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <string>
#include <vector>
#include<Parameter.h>
using namespace std;

class PDESolver
{

public:
    PDESolver(int ver){this->version=ver;}
    void set_iter(int iter_time);
    void show_iter();
    Parameter* get_pde_param_ptr();
    double interp(vector<double> x,vector<double> y,double ind);
    void solve();

    ~PDESolver(){}

    Parameter* pde_param_ptr;
    
private:
    int version;
    int iteration;
    
    
};  

#endif