#include<iostream>
#include<cmath>
#include <Eigen/Dense>  //https://www.cnblogs.com/rainbow70626/p/8819119.html 使用简介
#include <Eigen/LU>
#include <ctime>
using namespace Eigen;
using namespace std;


// void tridiag(int n, const double* a, const double* b, const double* c, double* u, const double* r)
// {
//     //n: Size of matrix.
//     //a: Below diagonal.
//     //b：Diagonal
//     //c: Above diagonal
//     //u: return result
//     //r: right
//     int j;
//     double bet;
//     MatrixXd gam(n,1);
//     if (b[0] == 0) throw " tridiag zero pivot";
//         u[0] = r[0] / (bet = b[0]);
    
//     for (j = 1; j < n; j++) 
//     {
//         gam(j,0) = c[j - 1] / bet;
//         bet = b[j] - a[j] * gam(j,0);
        
//         if (bet == 0) throw " tridiag zero pivot";
//         u[j] = (r[j] - a[j] * u[j - 1]) / bet;
//     }
//     for (j = n - 2; j >= 0; j--)
//         u[j] -= gam(j + 1,0) * u[j + 1];
// }

MatrixXd tridiag(int n, MatrixXd a, MatrixXd b, MatrixXd c, MatrixXd u, MatrixXd r)
{
    //n: Size of matrix.
    //a: Below diagonal.
    //b：Diagonal
    //c: Above diagonal
    //u: return result
    //r: right
    int j;
    double bet;
    MatrixXd gam(n,1);
    if(b(0,0)==0) throw "tridiag zero pivot";
        u(0,0) = r(0,0)/(bet = b(0,0));
    cout<<"if if if "<<endl;
    for (j=1;j<n;j++)
    {
        cout<<"*******************"<<endl;
        cout<<j<<" c(j-1,0) = "<<c(j-1,0)<<endl;
        gam(j,0) = c(j-1,0)/bet;
        bet = b(j,0) - a(j,0)*gam(j,0);
        cout<<j<<" "<<"b(j,0) = "<<b(j,0)<<" a(j,0) = "<<a(j,0)<<" gam(j,0) = "<<gam(j,0)<<endl;
        cout<<j<<" bet = "<<bet<<endl;
        if(bet==0) throw " tridiag zero pivot";
        u(j,0) = (r(j,0)-a(j,0)*u(j-1,0))/bet;
        cout<<"*******************"<<endl;
    }
    for (j=n-2;j>=0;j--)
        u(j,0) -= gam(j+1,0)*u(j+1,0);
    return u;
}

int main()
{
    int n = 5;
    MatrixXd a(n,1);
    MatrixXd b(n,1);
    MatrixXd c(n,1);
    MatrixXd u(n,1);
    MatrixXd r(n,1);

    a<<0,6,31,2,1; //a: Below diagonal. 前段补0
    b<<1,2,3,4,5;
    // b<<1,2,3,4,0;
    c<<13,7,9,7,0;  //c: Above diagonal 末尾补0
    r<<7,8,2,9,12;
    cout<<"Below diagonal a = \n"<<a.transpose()<<endl;
    cout<<"diagonal b = \n"<<b.transpose()<<endl;
    cout<<"Above diagonal c = \n"<<c.transpose()<<endl;
    cout<<"r = \n"<<r.transpose()<<endl;
    MatrixXd res = tridiag( n,  a,  b,  c,  u,  r);
    cout<<"tridiag res = \n"<<res.transpose()<<endl;
    MatrixXd mat(n,n);
    mat<<1,13,0,0,0,
    0,2,7,0,0,
    0,31,3,9,0,
    0,0,2,4,0,
    0,0,0,1,0;
    cout<<"mat = \n"<<mat<<endl;
    cout<<"mat^(-1)*r = \n"<<(mat.inverse()*r).transpose()<<endl;
    cout<<"LU res = \n"<< mat.partialPivLu().solve(r).transpose()<<endl;
}