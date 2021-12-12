#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <string>
#include <vector>
using namespace std;
class PDESolver
{

public:
    PDESolver(int ver){this->version=ver;}
    ~PDESolver(){}

private:
    int version;

};  

#endif