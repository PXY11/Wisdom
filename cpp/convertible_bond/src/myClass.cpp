#include <myClass.h>
#include <iostream>
using namespace std;
/*
Parameter::Parameter(double T)
{
    this-> T = T;
}
*/
void Parameter::show()
{
    cout<<"The value of T is:"<<this->T<<endl;
    cout<<"The value of Nt is:"<<this->Nt<<endl;
}

