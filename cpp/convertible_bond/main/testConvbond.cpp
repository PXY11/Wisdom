#include<iostream>
#include<vector>
#include<Parameter.h>
#include<PDESolver.h>
//#include<dataStructure.h>
using namespace std;

int main(){
    /*
    Parameter param(1);
    param.readParam();
    param.calParam();
    param.setBoundaryParam();
    param.show();
    */
    PDESolver pdesolver(1);
    pdesolver.set_iter(9);
    pdesolver.show_iter();
    Parameter* param_ptr;
    param_ptr = pdesolver.get_pde_param_ptr();

    pdesolver.pde_param_ptr->readParam();
    pdesolver.pde_param_ptr->calParam();
    pdesolver.pde_param_ptr->setBoundaryParam();
    //pdesolver.pde_param_ptr->show();
    cout<<"rhopenltycall="<<pdesolver.pde_param_ptr->rhopenltycall<<endl;
    //cout<<pdesolver.pde_param_ptr->dt<<endl;
/*
    pde_param_ptr->readParam();
    pde_param_ptr->calParam();
    pde_param_ptr->setBoundaryParam();
    pde_param_ptr->show();
    */
    return 0;
}