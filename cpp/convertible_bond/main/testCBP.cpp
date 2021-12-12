#include<iostream>
#include<vector>
#include<Parameter.h>
//#include<dataStructure.h>
using namespace std;

int main(){
    Parameter param(1);
    param.readParam();
    param.calParam();
    param.show();
    return 0;
}