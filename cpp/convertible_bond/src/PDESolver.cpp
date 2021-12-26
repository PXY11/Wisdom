#include<PDESolver.h>
#include<iostream>
using namespace std;


void PDESolver::set_iter(int iter_time)
{
    this->iteration = iter_time;
}

void PDESolver::show_iter()
{
    cout<<"iter time = "<<this->iteration<<endl;
}

Parameter* PDESolver::get_pde_param_ptr()
{
    this->pde_param_ptr = new Parameter(11);
    return this->pde_param_ptr;
}