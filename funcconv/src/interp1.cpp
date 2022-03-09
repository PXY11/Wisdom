/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * interp1.cpp
 *
 * Code generation for function 'interp1'
 *
 */

/* Include files */
#include "interp1.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "rt_nonfinite.h"
#include <cstring>

/* Function Declarations */
static void interp1Linear(const emxArray_real_T *y, const double xi[2401],
  double yi[2401], const double varargin_1[2401]);

/* Function Definitions */
static void interp1Linear(const emxArray_real_T *y, const double xi[2401],
  double yi[2401], const double varargin_1[2401])
{
  double minx;
  double maxx;
  int k;
  int low_i;
  int low_ip1;
  int high_i;
  double r;
  int mid_i;
  double d;
  minx = varargin_1[0];
  maxx = varargin_1[2400];

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(low_i,low_ip1,high_i,r,mid_i,d)

  for (k = 0; k < 2401; k++) {
    if (rtIsNaN(xi[k])) {
      yi[k] = rtNaN;
    } else {
      if ((!(xi[k] > maxx)) && (!(xi[k] < minx))) {
        low_i = 1;
        low_ip1 = 2;
        high_i = 2401;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (xi[k] >= varargin_1[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        r = varargin_1[low_i - 1];
        r = (xi[k] - r) / (varargin_1[low_i] - r);
        if (r == 0.0) {
          yi[k] = y->data[low_i - 1];
        } else if (r == 1.0) {
          yi[k] = y->data[low_i];
        } else {
          d = y->data[low_i - 1];
          if (d == y->data[low_i]) {
            yi[k] = d;
          } else {
            yi[k] = (1.0 - r) * d + r * y->data[low_i];
          }
        }
      }
    }
  }
}

void b_interp1(const double varargin_1[2401], const emxArray_real_T *varargin_2,
               const double varargin_3[2401], double Vq[2401])
{
  emxArray_real_T *y;
  int i;
  int n;
  double x[2401];
  int exitg1;
  double xtmp;
  int nd2;
  int y_tmp;
  emxInit_real_T(&y, 1);
  i = y->size[0];
  y->size[0] = varargin_2->size[0];
  emxEnsureCapacity_real_T(y, i);
  n = varargin_2->size[0];
  for (i = 0; i < n; i++) {
    y->data[i] = varargin_2->data[i];
  }

  std::memcpy(&x[0], &varargin_1[0], 2401U * sizeof(double));
  i = 0;
  do {
    exitg1 = 0;
    if (i < 2401) {
      if (rtIsNaN(varargin_1[i])) {
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      if (varargin_1[1] < varargin_1[0]) {
        for (i = 0; i < 1200; i++) {
          xtmp = x[i];
          x[i] = x[2400 - i];
          x[2400 - i] = xtmp;
        }

        if ((varargin_2->size[0] != 0) && (varargin_2->size[0] > 1)) {
          n = varargin_2->size[0] - 1;
          nd2 = varargin_2->size[0] >> 1;
          for (i = 0; i < nd2; i++) {
            xtmp = y->data[i];
            y_tmp = n - i;
            y->data[i] = y->data[y_tmp];
            y->data[y_tmp] = xtmp;
          }
        }
      }

      for (i = 0; i < 2401; i++) {
        Vq[i] = rtNaN;
      }

      interp1Linear(y, varargin_3, Vq, x);
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_real_T(&y);
}

double interp1(const double varargin_1[6], const double varargin_2[6], double
               varargin_3)
{
  double Vq;
  int low_i;
  double y[6];
  double x[6];
  int exitg1;
  double xtmp;
  int low_ip1;
  int high_i;
  int mid_i;
  for (low_i = 0; low_i < 6; low_i++) {
    y[low_i] = varargin_2[low_i];
    x[low_i] = varargin_1[low_i];
  }

  low_i = 0;
  do {
    exitg1 = 0;
    if (low_i < 6) {
      if (rtIsNaN(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      if (varargin_1[1] < varargin_1[0]) {
        xtmp = x[0];
        x[0] = x[5];
        x[5] = xtmp;
        xtmp = y[0];
        y[0] = y[5];
        y[5] = xtmp;
        xtmp = x[1];
        x[1] = x[4];
        x[4] = xtmp;
        xtmp = y[1];
        y[1] = y[4];
        y[4] = xtmp;
        xtmp = x[2];
        x[2] = x[3];
        x[3] = xtmp;
        xtmp = y[2];
        y[2] = y[3];
        y[3] = xtmp;
      }

      Vq = rtNaN;
      if ((!rtIsNaN(varargin_3)) && (!(varargin_3 > x[5])) && (!(varargin_3 < x
            [0]))) {
        low_i = 1;
        low_ip1 = 2;
        high_i = 6;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (varargin_3 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        xtmp = x[low_i - 1];
        xtmp = (varargin_3 - xtmp) / (x[low_i] - xtmp);
        if (xtmp == 0.0) {
          Vq = y[low_i - 1];
        } else if (xtmp == 1.0) {
          Vq = y[low_i];
        } else if (y[low_i - 1] == y[low_i]) {
          Vq = y[low_i - 1];
        } else {
          Vq = (1.0 - xtmp) * y[low_i - 1] + xtmp * y[low_i];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return Vq;
}

/* End of code generation (interp1.cpp) */
