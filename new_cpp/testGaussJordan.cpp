#include<iostream>
#include "dataStructure.h"
using namespace std;

#define SWAP(a,b)   {temp=(a); (a)=(b); (b)=temp;}
#define EPS 1e-8
void GaussJordan(double** A, double** B, long n, long m)
{
	// ����˵����
	// a. ���Է�����ϵ������A, (n)*(n)
	// n. ����A��ʵ����Чά�ȣ�n
	// b. ���Է�����rhs: (n)*m,  ��ʾ���Ƿ����Ҷ�	// m. ��ʾrhs���ɼ��й��ɣ�
	//    Ax = B: ֻ���һ�����Է����飺m=1��ʾ������Է����飬��ʱb����̬Ϊb={{b1},......,{bn}}�� 
	//    AX = B: ͬʱ���m�����Է����飺��ʱb����̬ΪB = {{b11,...,b1m},......,{bn1,...,bnm}}�� 
	//    A���棺  m=n��ʾ���棻 ��ʱb��n*n�ĵ�λ����
	//--------------------------------------------------------------------------------------------
	//ȫѡ��ԪGJ�㷨���ܣ�
	//		  1. ѡȡ����������Ԫ�أ������Ƶ��Խ����ϣ�����aii, ��������ͬʱ����aii���Խ���ֵ��Ϊ1��ע����[A,b]�����������������������ͬ���Ĳ���
	//        2. ��Ԫ������*(-aji), j!=i, �ӵ���j�У�����i�����г�aii��λ�û�Ϊ0��
	//        3. �����ڻ�ȥ��i�У���i�е�����λ��Ѱ��Ԫ������������Ĳ��������ֱ�����жԽ�Ԫ�ػ�Ϊ1���ǶԽ�Ԫ�ػ�Ϊ0����ʱ�Ҷ˵ľ���Ϊ���Է�����Ľ�������

	//������
	//	1. ������Է�����
	//  double A[3][3] = { {0,1,2},{1,1,4},{2,-1,0} };
	//  double b[3][1] = { {1},{1},{1} };
	//  long n = 3, m = 1;
	//  GaussJordan((double**)A, (double**)b, n, m);
	//-------------------------------------------------------
	// 2. ����	
	// double B[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	//  n = 3, m = 3;
	//  GaussJordan((double**)A, (double**)B, n, m);
	// ------------------------------------------------------
	// 3. ������
	// 	for (int i = 0; i < n; i++)
	//  {
    //      for (int j = 0; j < m; j++)
	//		    cout << B[i][j] << ",";
	//	    cout << endl;
	//  }
	//
	
	//---����ϵ��������β�A��a�У����aҪ�ͷſռ�---
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

	//---�����м����---
	long *indxc, *indxr, *ipiv;
	indxc = new long[n];
	indxr = new long[n];
	ipiv = new long[n];
	long i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv, temp;//big��ʾ��Ԫֵ��dum, pivinv��ʾ����ȷ������Ԫֵ�ĵ���,
	for (j = 0; j < n; j++) ipiv[j] = 0;//�Խ��ߵ�ǰλ���Ƿ�ѡ����Ԫ��1��ʾѡ����Ԫ��Ҳ���������г��˶Խ���Ϊ1������λ�þ�Ϊ0.

	//------------step1: ȫѡ��Ԫ------------------------
	for (i = 0; i < n; i++)
	{
		big = 0.0;
		//ȫѡ��Ԫ�� Ѱ�Ҿ���a�����ֵ����big, �����������ֱ����j,k
		for (j = 0; j < n; j++)//ѭ��ÿһ��j
			if (ipiv[j] != 1)//�Ѿ����������Ԫ���У�����һ�μ����в�������Ӧ������
				for (k = 0; k < n; k++) //ѭ��ÿһ��k
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big) //Ѱ�ҵ�j�е����ֵ
						{
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						throw("gaussj: Singular Matrix-1 !!");
				}
		++(ipiv[icol]);//big value��Ӧ��ipiv+1����ʾ��ѡ����Ԫ

		//------------step2: ��Ԫ���ڶԽ��ߣ����任�����Խ���------------------------
		if (irow != icol) //���ڶԽ����ϣ�������Ӧ���У���big value�Ƶ��Խ�����
		{
			for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l])
		}

		//------------step3: ��Ԫ�����б�׼����*/aii��------------------------
		indxr[i] = irow;//��Ԫԭʼλ��
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)//Ŀǰa[icol][icol]Ϊ��Ԫ��[icol][icol]Ϊ��Ԫ����λ��
			a[icol][icol] = EPS;
		//          throw("gaussj: Singular Matrix-2 !!");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;//��Ԫ�����г�����Ԫ��Ϊ1
		for (l = 0; l < n; l++) a[icol][l] *= pivinv;//��Ԫ�����б�׼����������a[icol][icol]��������/�����˾���Ҫ��������
		for (l = 0; l < m; l++) b[icol][l] *= pivinv;


		//------------step4: ��ȥ����Ԫ���ڵ�icol���ϵ�ֵ�������Ϊ0---------------
		for (ll = 0; ll < n; ll++)
			if (ll != icol)
			{
				dum = a[ll][icol];//-dum*��Ԫ��+������������
				a[ll][icol] = 0.0;//����Ԫ����icol�е�ֵ����Ϊ0��
				for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;//�����ض��ķ���Ԫ�У�����ÿһ�У�
				for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
			}
	}
	

	//�洢�����ԭ����B��
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < m; j++) *((double*)B + i * m + j) = b[i][j];
	}

	//�ͷ��м����ָ��
	delete[] indxc;
	delete[] indxr;
	delete[] ipiv;

	//�ͷŶ�ά����ָ��
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
	// ͬʱ��������Է�������������
	// AX = B, A(n*n)Ϊ����BΪn*n��λ�������棬BΪn*m��ͬʱ�����ͬϵ�����Է�����
	// ���ؽ��������B��, Ϊ����A���������в����ı䣬�����м����a == A
	//--------------------------------------------------------------------------------------------
	//ȫѡ��ԪGJ�㷨���ܣ�
	//		  1. ѡȡ����������Ԫ�أ������Ƶ��Խ����ϣ�����aii, ��������ͬʱ����aii���Խ���ֵ��Ϊ1��ע����[A,B]�����������������������ͬ���Ĳ���
	//        2. ��Ԫ������*(-aji), j!=i, �ӵ���j�У�����i�����г�aii��λ�û�Ϊ0��
	//        3. �����ڻ�ȥ��i�У���i�е�����λ��Ѱ��Ԫ������������Ĳ��������ֱ�����жԽ�Ԫ�ػ�Ϊ1���ǶԽ�Ԫ�ػ�Ϊ0����ʱ�Ҷ˵ľ���Ϊ���Է�����Ľ�������
	//
	// ʾ����
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
	Array2D a = A; //Ϊ�˲��ı�Aԭʼ��ֵ������A���ջ��Ϊ��λ��


	//---�����м����---
	Array1L indxc(n), indxr(n), ipiv(n);
	long i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv, temp;//big��ʾ��Ԫֵ��dum, pivinv��ʾ����ȷ������Ԫֵ�ĵ���,
	for (j = 0; j < n; j++) ipiv[j] = 0;//�Խ��ߵ�ǰλ���Ƿ�ѡ����Ԫ��1��ʾѡ����Ԫ��Ҳ���������г��˶Խ���Ϊ1������λ�þ�Ϊ0.

	//------------step1: ȫѡ��Ԫ------------------------
	for (i = 0; i < n; i++)
	{
		big = 0.0;
		//ȫѡ��Ԫ�� Ѱ�Ҿ���a�����ֵ����big, �����������ֱ����j,k
		for (j = 0; j < n; j++)//ѭ��ÿһ��j
			if (ipiv[j] != 1)//�Ѿ����������Ԫ���У�����һ�μ����в�������Ӧ������
				for (k = 0; k < n; k++) //ѭ��ÿһ��k
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big) //Ѱ�ҵ�j�е����ֵ
						{
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						throw("gaussj: Singular Matrix-1 !!");
				}
		++(ipiv[icol]);//big value��Ӧ��ipiv+1����ʾ��ѡ����Ԫ

		//------------step2: ��Ԫ���ڶԽ��ߣ����任�����Խ���------------------------
		if (irow != icol) //���ڶԽ����ϣ�������Ӧ���У���big value�Ƶ��Խ�����
		{
			for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 0; l < m; l++) SWAP(B[irow][l], B[icol][l])
		}

		//------------step3: ��Ԫ�����б�׼����*/aii��------------------------
		indxr[i] = irow;//��Ԫԭʼλ��
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)//Ŀǰa[icol][icol]Ϊ��Ԫ��[icol][icol]Ϊ��Ԫ����λ��
			a[icol][icol] = EPS;
		//          throw("gaussj: Singular Matrix-2 !!");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;//��Ԫ�����г�����Ԫ��Ϊ1
		for (l = 0; l < n; l++) a[icol][l] *= pivinv;//��Ԫ�����б�׼����������a[icol][icol]��������/�����˾���Ҫ��������
		for (l = 0; l < m; l++) B[icol][l] *= pivinv;


		//------------step4: ��ȥ����Ԫ���ڵ�icol���ϵ�ֵ�������Ϊ0---------------
		for (ll = 0; ll < n; ll++)
			if (ll != icol)
			{
				dum = a[ll][icol];//-dum*��Ԫ��+������������
				a[ll][icol] = 0.0;//����Ԫ����icol�е�ֵ����Ϊ0��
				for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;//�����ض��ķ���Ԫ�У�����ÿһ�У�
				for (l = 0; l < m; l++) B[ll][l] -= B[icol][l] * dum;
			}
	}

}


int main()
{
	//---------------------GJ����������ΪC++ԭʼ��������-----------------------
	//1.���Է��������
	double A[3][3] = { {0,1,2},{1,1,4},{2,-1,0} };
	double b[3][1] = { {1},{1},{1} };
	long n = 3, m = 1;
	
	GaussJordan((double**)A, (double**)b, n, m);//����A�������ϸ�Ϊdouble**���ͣ���Ҫǿ��ת��
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << b[i][j] <<",";
		cout << endl;
	}


	//2.��������
	double B[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	n = 3, m = 3;
	GaussJordan((double**)A, (double**)B, n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << B[i][j] << ",";
		cout << endl;
	}

	//---------------------GJ����������ΪArray2D����-----------------------
	Array2D AA = { {0,1,2},{1,1,4},{2,-1,0} };
	Array2D BB = { {1,0,0},{0,1,0},{0,0,1} };
	GaussJordan(AA, BB);
	cout << BB;
	cout << AA;

	return 0;
}