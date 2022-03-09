/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * func1.cpp
 *
 * Code generation for function 'func1'
 *
 */

/* Include files */
#include "func1.h"
#include <cmath>
#include <cstring>
#include<iostream>
using namespace std;
/* Function Definitions */
void func1(double m, double b, double L, double n, double H, double Q[4], const
           double E[4], const double u[4], const double C0[4], double CC[16])
{
  double l;
  double h;
  double s;
  double a[4];
  int i;
  int b_i;
  double e[4];
  int jp1j;
  double A[16];
  int jy;
  double a_tmp[16];
  signed char ipiv[4];
  int j;
  signed char p[4];
  int mmj_tmp;
  int b_b;
  int jj;
  int iy;
  int ix;
  signed char i1;
  double D2[4];
  double D3_idx_1;
  double d;
  double D3_idx_2;
  double D3_idx_3;
  double D3_idx_0;

  /* ����Ĵ��벻�ù� */
  l = L / n;
  h = H / n;
  s = b * l;
  Q[0] /= s;
  a[0] = 0.0;
  Q[1] /= s;
  a[1] = 0.0;
  Q[2] /= s;
  a[2] = 0.0;
  Q[3] /= s;
  a[3] = 0.0;
  i = static_cast<int>(m);
  for (b_i = 0; b_i < i; b_i++) {
    a[b_i] = u[b_i] * h;
  }

  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 0.0;
  e[3] = 0.0;
  jp1j = static_cast<int>((m - 1.0));
  for (b_i = 0; b_i < jp1j; b_i++) {
    e[b_i] = E[b_i] * l / h;
  }

  std::memset(&A[0], 0, 16U * sizeof(double));
  A[0] = a[0] + e[0];
  jy = i - 1;
  A[(i + (jy << 2)) - 1] = a[jy] + e[jp1j - 1];
  A[4] = -e[0];
  jy = i - 2;
  A[(i + (jy << 2)) - 1] = -e[jy];
  for (b_i = 0; b_i <= i - 3; b_i++) {
    h = e[b_i + 1];
    A[(b_i + ((b_i + 1) << 2)) + 1] = (a[b_i + 1] + e[b_i]) + h;
    A[(b_i + (b_i << 2)) + 1] = -e[b_i];
    A[(b_i + ((b_i + 2) << 2)) + 1] = -h;
  }

  e[0] = a[0] * C0[0];
  e[1] = a[1] * C0[1];
  e[2] = a[2] * C0[2];
  e[3] = a[3] * C0[3];
  e[0] += Q[0] * l;
  std::memset(&a_tmp[0], 0, 16U * sizeof(double));
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    mmj_tmp = 2 - j;
    b_b = j * 5;
    jj = j * 5;
    jp1j = b_b + 2;
    jy = 4 - j;
    iy = 0;
    ix = b_b;
    h = std::abs(A[jj]);
    for (b_i = 2; b_i <= jy; b_i++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > h) {
        iy = b_i - 1;
        h = s;
      }
    }

    if (A[jj + iy] != 0.0) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = static_cast<signed char>((iy + 1));
        h = A[j];
        A[j] = A[iy];
        A[iy] = h;
        ix = j + 4;
        iy += 4;
        h = A[ix];
        A[ix] = A[iy];
        A[iy] = h;
        ix += 4;
        iy += 4;
        h = A[ix];
        A[ix] = A[iy];
        A[iy] = h;
        ix += 4;
        iy += 4;
        h = A[ix];
        A[ix] = A[iy];
        A[iy] = h;
      }

      i = (jj - j) + 4;
      for (b_i = jp1j; b_i <= i; b_i++) {
        A[b_i - 1] /= A[jj];
      }
    }

    jy = b_b + 4;
    iy = jj;
    for (b_i = 0; b_i <= mmj_tmp; b_i++) {
      h = A[jy];
      if (A[jy] != 0.0) {
        ix = jj + 1;
        i = iy + 6;
        jp1j = (iy - j) + 8;
        for (b_b = i; b_b <= jp1j; b_b++) {
          A[b_b - 1] += A[ix] * -h;
          ix++;
        }
      }

      jy += 4;
      iy += 4;
    }
  }

  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 4;
  if (ipiv[0] > 1) {
    jy = ipiv[0] - 1;
    iy = p[jy];
    p[jy] = 1;
    p[0] = static_cast<signed char>(iy);
  }

  if (ipiv[1] > 2) {
    jy = ipiv[1] - 1;
    iy = p[jy];
    p[jy] = p[1];
    p[1] = static_cast<signed char>(iy);
  }

  if (ipiv[2] > 3) {
    jy = ipiv[2] - 1;
    iy = p[jy];
    p[jy] = p[2];
    p[2] = static_cast<signed char>(iy);
  }

  i1 = p[0];
  a_tmp[(p[0] - 1) << 2] = 1.0;
  for (j = 1; j < 5; j++) {
    if (a_tmp[(j + ((i1 - 1) << 2)) - 1] != 0.0) {
      i = j + 1;
      for (b_i = i; b_i < 5; b_i++) {
        jy = (b_i + ((i1 - 1) << 2)) - 1;
        a_tmp[jy] -= a_tmp[(j + ((i1 - 1) << 2)) - 1] * A[(b_i + ((j - 1) << 2))
          - 1];
      }
    }
  }

  i1 = p[1];
  a_tmp[((p[1] - 1) << 2) + 1] = 1.0;
  for (j = 2; j < 5; j++) {
    if (a_tmp[(j + ((i1 - 1) << 2)) - 1] != 0.0) {
      i = j + 1;
      for (b_i = i; b_i < 5; b_i++) {
        jy = (b_i + ((i1 - 1) << 2)) - 1;
        a_tmp[jy] -= a_tmp[(j + ((i1 - 1) << 2)) - 1] * A[(b_i + ((j - 1) << 2))
          - 1];
      }
    }
  }

  i1 = p[2];
  a_tmp[((p[2] - 1) << 2) + 2] = 1.0;
  for (j = 3; j < 5; j++) {
    jp1j = (j + ((i1 - 1) << 2)) - 1;
    if (a_tmp[jp1j] != 0.0) {
      i = j + 1;
      for (b_i = i; b_i < 5; b_i++) {
        jy = ((i1 - 1) << 2) + 3;
        a_tmp[jy] -= a_tmp[jp1j] * A[((j - 1) << 2) + 3];
      }
    }
  }

  a_tmp[((p[3] - 1) << 2) + 3] = 1.0;
  for (j = 0; j < 4; j++) {
    iy = j << 2;
    d = a_tmp[iy + 3];
    if (d != 0.0) {
      a_tmp[iy + 3] = d / A[15];
      for (b_i = 0; b_i < 3; b_i++) {
        jy = b_i + iy;
        a_tmp[jy] -= a_tmp[iy + 3] * A[b_i + 12];
      }
    }

    d = a_tmp[iy + 2];
    if (d != 0.0) {
      a_tmp[iy + 2] = d / A[10];
      for (b_i = 0; b_i < 2; b_i++) {
        jy = b_i + iy;
        a_tmp[jy] -= a_tmp[iy + 2] * A[b_i + 8];
      }
    }

    d = a_tmp[iy + 1];
    if (d != 0.0) {
      a_tmp[iy + 1] = d / A[5];
      for (b_i = 0; b_i < 1; b_i++) {
        a_tmp[iy] -= a_tmp[iy + 1] * A[4];
      }
    }

    if (a_tmp[iy] != 0.0) {
      a_tmp[iy] /= A[0];
    }

    D2[j] = a[j] * e[j];
  }

  D2[0] += Q[1] * l;
  D3_idx_1 = a[1] * D2[1];
  D3_idx_2 = a[2] * D2[2];
  D3_idx_3 = a[3] * D2[3];
  D3_idx_0 = a[0] * D2[0] + Q[2] * l;
  a[0] *= D3_idx_0;
  a[1] *= D3_idx_1;
  a[2] *= D3_idx_2;
  a[3] *= D3_idx_3;
  a[0] += Q[3] * l;
  for (i = 0; i < 4; i++) {
    d = a_tmp[i + 4];
    h = a_tmp[i + 8];
    s = a_tmp[i + 12];
    CC[i] = ((a_tmp[i] * e[0] + d * e[1]) + h * e[2]) + s * e[3];
    CC[i + 4] = ((a_tmp[i] * D2[0] + d * D2[1]) + h * D2[2]) + s * D2[3];
    CC[i + 8] = ((a_tmp[i] * D3_idx_0 + d * D3_idx_1) + h * D3_idx_2) + s *
      D3_idx_3;
    CC[i + 12] = ((a_tmp[i] * a[0] + d * a[1]) + h * a[2]) + s * a[3];
  }
  // for(int i=0;i<16;i++)
  // {
  //   cout<<CC[i]<<endl;
  // }
  // cout<<"e0: "<<e[0];
  cout<<"in func1"<<endl;
  cout<<"p0: "<<p[0]<<"p1: "<<p[1]<<endl;
  cout<<"A"<<A[0]<<"e0"<<e[0]<<endl;
  cout<<"func1 end"<<endl;
  /* ����Ĵ��벻�ù� */
}

/* End of code generation (func1.cpp) */
