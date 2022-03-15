#include <iostream>
#include <ctime>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;
/***********************
* solve equation: matrix_NN * x = v_Nd
************************/
const int MATRIX_SIZE = 99;
// https://blog.csdn.net/weixin_41074793/article/details/84241776
int main()
{
        Matrix< double, MATRIX_SIZE, MATRIX_SIZE > matrix_NN;
        matrix_NN = MatrixXd::Random( MATRIX_SIZE, MATRIX_SIZE );
        Matrix<double, MATRIX_SIZE, 1> v_Nd;
        v_Nd = MatrixXd::Random( MATRIX_SIZE,1);
        clock_t time_stt = clock();

        //直接求逆
        Matrix< double, MATRIX_SIZE, 1> x = matrix_NN.inverse()*v_Nd;
        cout << "time use in normal inverse is       " << 1000.0 * (clock() - time_stt) /
                                (double)CLOCKS_PER_SEC << "ms" << endl;
        //Qr分解
        time_stt = clock();
        x = matrix_NN.colPivHouseholderQr().solve(v_Nd);
        cout << "time use in Qr composition is       " << 1000 * (clock() - time_stt) / (double)
            CLOCKS_PER_SEC << "ms" << endl;

        //cholesky分解
        time_stt = clock();
        // 使得NN成为正定矩阵，才能分解
        matrix_NN = matrix_NN.transpose() * matrix_NN;
        x = matrix_NN.partialPivLu().solve(v_Nd);
        cout << "time use in cholesky composition is " << 1000 * (clock() - time_stt) / (double)
            CLOCKS_PER_SEC << "ms" << endl;

        //LU分解
        time_stt = clock();
        x = matrix_NN.partialPivLu().solve(v_Nd);
        cout << "time use in LU composition is       " << 1000 * (clock() - time_stt) / (double)
            CLOCKS_PER_SEC << "ms" << endl;

    return 0;
}



