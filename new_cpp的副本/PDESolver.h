#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <string>
#include <vector>
#include"Parameter.h"
// #define N 2
// #define M 2
using namespace std;

class PDESolver
{

public:
    PDESolver(Parameter* _paras, int _iters)
    {
        pde_param_ptr = _paras;
        iteration = _iters;
    }
    void set_iter(int iter_time);
    void show_iter();
    //Parameter* get_pde_param_ptr(); //返回一个Parameter类实例的指针，用于调用实例中的属性
    // double interpCB(vector<double> x,vector<double> y,double ind); //插值函数 
    double interpCB(MatrixXd x,MatrixXd y,double ind);
    double solve(); //解PDE，所有操作均是对get_pde_param_ptr()返回的实例去做操作，
                  //最终计算结果会保存在Parameter实例的属性u中，根据索引可查得对应的价格
    // double getA(double arcs[N][N],int n);
    // void  getAStart(double arcs[N][N],int n,double ans[N][N]);
    // bool GetMatrixInverse(double src[N][N],int n,double des[N][N]);
    // vector<vector<double>*>  MatMul(double left[N][M],double right[M][N]);

    double calBpstar();
    double calBcstar();
    ~PDESolver(){}

    Parameter* pde_param_ptr;
    
private:
    int iteration; //解PDE迭代次数
    
    
};  

#endif