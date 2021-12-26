#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <string>
#include <vector>
using namespace std;

class PDESolver
{

public:
    PDESolver(int ver){this->version=ver;}
    void set_iter(int iter_time);
    void show_iter();
    ~PDESolver(){}

private:
    int version;
    int iteration;
};  

#endif