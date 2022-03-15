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
    cout<<"minv X n \n"<<minv*n<<endl;
    // cout<<"m左除 X n \n"<<m\n<<endl;
    cout<<"******************"<<endl;
    cout<<"3*m \n"<<3*m<<endl;
    double musize =13;
    MatrixXd mu(int(musize),1);
    mu(0,0) = 998;
    cout << "Here is the matrix mu:\n" << mu.transpose() << endl;
    cout<<mu.transpose()*mu<<endl;
    
    VectorXd v(3);
    v(0) = 4;
    v(1) = v(0) - 1;
    // v<<5;
    cout << "Here is the vector v:\n" << v << endl;
    // cout <<"max value of v = "<<max(v)<<endl;
    // cout << "Here is the matrix m:\n" << m << endl;
    // cout<<"m X v \n"<<m*v<<endl;
    // // v<<9;
    // cout<<"v push 9 \n"<<v<<endl;
    cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    MatrixXd P1(2,2);
    MatrixXd P2(2,2);
    P1<<1,1,1,1;
    P2<<1,1,1,1;
    bool sameP1=1;
    for(int i=0;i<2;i++)
    {
        for(int j =0;j<2;j++)
        {
            if (P2(i,j)!=P1(i,j))
            {   sameP1 = 0;
                break;
            }
        }
    }
    cout<<"samP1 = "<<sameP1;
}
// https://blog.csdn.net/qq_41854911/article/details/119814660 例子参考
// https://www.cnblogs.com/goingupeveryday/p/5699053.html 下载参考