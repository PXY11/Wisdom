#include <iostream>
#include <Eigen/Dense>

#include <vector>
using namespace Eigen;
using namespace std;
int main()
{
    MatrixXd m(2,2);
    m(0,0) = 1;
    m(1,0) = 2;
    m(0,1) = 1;
    m(1,1) = m(1,0) + m(0,1);
    MatrixXd n(2,2);
    n<<1,2,3,4;
    cout<<"Here is the matrix n \n"<<n<<endl;
    cout << "Here is the matrix m:\n" << m << endl;
    cout<<"m X n \n"<<n*m<<endl;
    cout<<"Inverse of m: \n"<<m.inverse()<<endl;
    MatrixXd minv;
    minv = m.inverse();
    cout << "Here is the matrix minv:\n" << minv << endl;
    cout<<"******************"<<endl;
    int ns = 24;
    MatrixXd mu(ns,ns);
    mu(0,0) = 998;
    // cout << "Here is the matrix mu:\n" << mu << endl;
    VectorXd v(ns);
    
    v(0) = 4;
    v(1) = v(0) - 1;
    v<<5;
    cout << "Here is the vector v:\n" << v << endl;
    // cout << "Here is the matrix m:\n" << m << endl;
    // cout<<"m X v \n"<<m*v<<endl;
    // // v<<9;
    // cout<<"v push 9 \n"<<v<<endl;

}
// https://blog.csdn.net/qq_41854911/article/details/119814660 例子参考
// https://www.cnblogs.com/goingupeveryday/p/5699053.html 下载参考