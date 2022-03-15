#include<iostream>
#include "dataStructure.h"
using namespace std;

#define SWAP(a,b)   {temp=(a); (a)=(b); (b)=temp;}
#define EPS 1e-8
void GaussJordan(double** A, double** B, long n, long m)
{
	// 参数说明：
	// a. 线性方程组系数矩阵：A, (n)*(n)
	// n. 矩阵A的实际有效维度：n
	// b. 线性方程组rhs: (n)*m,  表示的是方程右端	// m. 表示rhs是由几列构成，
	//    Ax = B: 只求解一个线性方程组：m=1表示求解线性方程组，此时b中形态为b={{b1},......,{bn}}； 
	//    AX = B: 同时求解m个线性方程组：此时b中形态为B = {{b11,...,b1m},......,{bn1,...,bnm}}； 
	//    A求逆：  m=n表示求逆； 此时b是n*n的单位矩阵
	//--------------------------------------------------------------------------------------------
	//全选主元GJ算法介绍：
	//		  1. 选取矩阵中最大的元素，将其移到对角线上，比如aii, 将所在行同时除以aii，对角线值化为1；注：对[A,b]增广矩阵处理，所以左右两端做同样的操作
	//        2. 主元所在行*(-aji), j!=i, 加到第j行，将第i列所有除aii的位置化为0；
	//        3. 继续在划去第i行，第i列的其他位置寻主元，并按照上面的步骤操作，直到所有对角元素化为1，非对角元素化为0，此时右端的矩阵即为线性方程组的解或逆矩阵

	//举例：
	//	1. 求解线性方程组
	//  double A[3][3] = { {0,1,2},{1,1,4},{2,-1,0} };
	//  double b[3][1] = { {1},{1},{1} };
	//  long n = 3, m = 1;
	//  GaussJordan((double**)A, (double**)b, n, m);
	//-------------------------------------------------------
	// 2. 求逆	
	// double B[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	//  n = 3, m = 3;
	//  GaussJordan((double**)A, (double**)B, n, m);
	// ------------------------------------------------------
	// 3. 输入结果
	// 	for (int i = 0; i < n; i++)
	//  {
    //      for (int j = 0; j < m; j++)
	//		    cout << B[i][j] << ",";
	//	    cout << endl;
	//  }
	//
	
	//---保存系数矩阵从形参A到a中，最后a要释放空间---
	double** a;
	double** b;
	a = new double*[n];
	b = new double* [n];
	for (long i = 0; i < n; i++)
	{
		a[i] = new double[n];
		b[i] = new double[m];
		for (long j = 0; j < n; j++) a[i][j] = *((double*)A + i * n + j);
		for (long j = 0; j < m; j++) b[i][j] = *((double*)B + i * m + j);
	}

	//---定义中间变量---
	long *indxc, *indxr, *ipiv;
	indxc = new long[n];
	indxr = new long[n];
	ipiv = new long[n];
	long i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv, temp;//big表示主元值，dum, pivinv表示最终确定的主元值的倒数,
	for (j = 0; j < n; j++) ipiv[j] = 0;//对角线当前位置是否被选过主元，1表示选过主元，也就是所在列除了对角线为1，其他位置均为0.

	//------------step1: 全选主元------------------------
	for (i = 0; i < n; i++)
	{
		big = 0.0;
		//全选主元： 寻找矩阵a的最大值放入big, 所在行列数分别放入j,k
		for (j = 0; j < n; j++)//循环每一行j
			if (ipiv[j] != 1)//已经计算过的主元行列，在下一次计算中不考虑相应的行列
				for (k = 0; k < n; k++) //循环每一列k
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big) //寻找第j行的最大值
						{
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						throw("gaussj: Singular Matrix-1 !!");
				}
		++(ipiv[icol]);//big value对应的ipiv+1；表示已选过主元

		//------------step2: 主元不在对角线，将其换行至对角线------------------------
		if (irow != icol) //不在对角线上，交换相应的行，将big value移到对角线上
		{
			for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l])
		}

		//------------step3: 主元所在行标准化（*/aii）------------------------
		indxr[i] = irow;//主元原始位置
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)//目前a[icol][icol]为主元，[icol][icol]为主元所在位置
			a[icol][icol] = EPS;
		//          throw("gaussj: Singular Matrix-2 !!");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;//主元所在行除以主元化为1
		for (l = 0; l < n; l++) a[icol][l] *= pivinv;//主元所在行标准化，均除以a[icol][icol]，方程左/右两端均需要这样处理
		for (l = 0; l < m; l++) b[icol][l] *= pivinv;


		//------------step4: 消去非主元行在第icol列上的值，将其变为0---------------
		for (ll = 0; ll < n; ll++)
			if (ll != icol)
			{
				dum = a[ll][icol];//-dum*主元行+其他各行上面
				a[ll][icol] = 0.0;//非主元行在icol列的值均化为0；
				for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;//对于特定的非主元行，修正每一列，
				for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
			}
	}
	

	//存储结果到原来的B中
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < m; j++) *((double*)B + i * m + j) = b[i][j];
	}

	//释放中间变量指针
	delete[] indxc;
	delete[] indxr;
	delete[] ipiv;

	//释放二维数组指针
	for (long i = 0; i < n; i++)
	{
		delete[] a[i];
		delete[] b[i];
	}
	delete[] a;
	delete[] b;

}
void GaussJordan(Array2D& A, Array2D& B)
{
	// 同时求解多个线性方程组或矩阵求逆
	// AX = B, A(n*n)为方阵，B为n*n单位阵则求逆，B为n*m则同时求解多个同系数线性方程组
	// 返回结果保存在B中, 为保持A在求解过程中不被改变，建立中间变量a == A
	//--------------------------------------------------------------------------------------------
	//全选主元GJ算法介绍：
	//		  1. 选取矩阵中最大的元素，将其移到对角线上，比如aii, 将所在行同时除以aii，对角线值化为1；注：对[A,B]增广矩阵处理，所以左右两端做同样的操作
	//        2. 主元所在行*(-aji), j!=i, 加到第j行，将第i列所有除aii的位置化为0；
	//        3. 继续在划去第i行，第i列的其他位置寻主元，并按照上面的步骤操作，直到所有对角元素化为1，非对角元素化为0，此时右端的矩阵即为线性方程组的解或逆矩阵
	//
	// 示例：
	// 	Array2D AA = { {0,1,2},{1,1,4},{2,-1,0} };
	//	Array2D BB = { {1,0,0},{0,1,0},{0,0,1} };
	//	GaussJordan(AA, BB);
	//	cout << BB;
	//	cout << AA;
	// 
	//



	long Ar = A.size(), Ac = A[0].size(), Br = B.size(), Bc = B[0].size();

	if (Ar != Ac || Ar != Br)throw"matrix_size_error!!!";

	long n = Ar, m = Bc;
	Array2D a = A; //为了不改变A原始的值，否则A最终会变为单位阵


	//---定义中间变量---
	Array1L indxc(n), indxr(n), ipiv(n);
	long i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv, temp;//big表示主元值，dum, pivinv表示最终确定的主元值的倒数,
	for (j = 0; j < n; j++) ipiv[j] = 0;//对角线当前位置是否被选过主元，1表示选过主元，也就是所在列除了对角线为1，其他位置均为0.

	//------------step1: 全选主元------------------------
	for (i = 0; i < n; i++)
	{
		big = 0.0;
		//全选主元： 寻找矩阵a的最大值放入big, 所在行列数分别放入j,k
		for (j = 0; j < n; j++)//循环每一行j
			if (ipiv[j] != 1)//已经计算过的主元行列，在下一次计算中不考虑相应的行列
				for (k = 0; k < n; k++) //循环每一列k
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big) //寻找第j行的最大值
						{
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						throw("gaussj: Singular Matrix-1 !!");
				}
		++(ipiv[icol]);//big value对应的ipiv+1；表示已选过主元

		//------------step2: 主元不在对角线，将其换行至对角线------------------------
		if (irow != icol) //不在对角线上，交换相应的行，将big value移到对角线上
		{
			for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 0; l < m; l++) SWAP(B[irow][l], B[icol][l])
		}

		//------------step3: 主元所在行标准化（*/aii）------------------------
		indxr[i] = irow;//主元原始位置
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)//目前a[icol][icol]为主元，[icol][icol]为主元所在位置
			a[icol][icol] = EPS;
		//          throw("gaussj: Singular Matrix-2 !!");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;//主元所在行除以主元化为1
		for (l = 0; l < n; l++) a[icol][l] *= pivinv;//主元所在行标准化，均除以a[icol][icol]，方程左/右两端均需要这样处理
		for (l = 0; l < m; l++) B[icol][l] *= pivinv;


		//------------step4: 消去非主元行在第icol列上的值，将其变为0---------------
		for (ll = 0; ll < n; ll++)
			if (ll != icol)
			{
				dum = a[ll][icol];//-dum*主元行+其他各行上面
				a[ll][icol] = 0.0;//非主元行在icol列的值均化为0；
				for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;//对于特定的非主元行，修正每一列，
				for (l = 0; l < m; l++) B[ll][l] -= B[icol][l] * dum;
			}
	}

}


int main()
{
	//---------------------GJ方法：参数为C++原始变量类型-----------------------
	//1.线性方程组求解
	double A[3][3] = { {0,1,2},{1,1,4},{2,-1,0} };
	double b[3][1] = { {1},{1},{1} };
	long n = 3, m = 1;
	
	GaussJordan((double**)A, (double**)b, n, m);//变量A并不是严格为double**类型，需要强制转化
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << b[i][j] <<",";
		cout << endl;
	}


	//2.求解逆矩阵
	double B[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	n = 3, m = 3;
	GaussJordan((double**)A, (double**)B, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << B[i][j] << ",";
		cout << endl;
	}

	//---------------------GJ方法：参数为Array2D类型-----------------------
	Array2D AA = { {0,1,2},{1,1,4},{2,-1,0} };
	Array2D BB = { {1,0,0},{0,1,0},{0,0,1} };
	GaussJordan(AA, BB);
	cout << BB;
	cout << AA;

	return 0;
}