#include<iostream>
#include<vector>
#include<Parameter.h>
#include<PDESolver.h>
//#include<dataStructure.h>
using namespace std;

int main(){
    Parameter param(1);
    param.readParam();
    param.calParam();
    param.setBoundaryParam();
    param.show();
    PDESolver pdesolver(1);
    pdesolver.set_iter(9);
    pdesolver.show_iter();
    return 0;
}