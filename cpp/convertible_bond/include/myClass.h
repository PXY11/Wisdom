#ifndef PARAMETER_H
#define PARAMETER_H
#include<math.h>
class Parameter
{

public:
    Parameter(double T)
    {
        this->T = T;
        this->Nt = ceil(T*500);
        
    }
    void show();
    ~Parameter(){}

private:
    double T;
    double Nt;

};

#endif