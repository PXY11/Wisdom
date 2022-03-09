#include<iostream>
#include<vector>
#include<Parameter.h>
#include<PDESolver.h>
//#include<dataStructure.h>
using namespace std;

int main(){

    PDESolver pdesolver(1);
    pdesolver.set_iter(9);
    pdesolver.show_iter();
    Parameter* param_ptr;
    param_ptr = pdesolver.get_pde_param_ptr();

    /* 测试convertible_code.m 功能*/
    pdesolver.pde_param_ptr->readParam();
    pdesolver.pde_param_ptr->calParam();
    pdesolver.pde_param_ptr->setBoundaryParam();
    double res = pdesolver.solve();
    cout<<"Final res = "<<res;
    // pdesolver.pde_param_ptr->show();  //用于输出pde parameters
    // cout<<"rhopenltycall="<<pdesolver.pde_param_ptr->rhopenltycall<<endl;
    //cout<<pdesolver.pde_param_ptr->dt<<endl;


    return 0;
}