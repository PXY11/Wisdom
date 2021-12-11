#ifndef PARAMETER_H
#define PARAMETER_H

class Parameter
{

public:
    Parameter(double T)
    {
        this->T = T;
    }
    void show();
    ~Parameter(){}


private:
    double T;

};

#endif