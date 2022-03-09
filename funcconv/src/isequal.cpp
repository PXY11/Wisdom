/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * isequal.cpp
 *
 * Code generation for function 'isequal'
 *
 */

/* Include files */
#include "isequal.h"
#include "funcconv.h"
#include "rt_nonfinite.h"

/* Function Definitions */
bool isequal(const emxArray_real_T *varargin_1_d, const emxArray_int32_T
             *varargin_1_colidx, const emxArray_int32_T *varargin_1_rowidx, int
             varargin_1_m, int varargin_1_n, const emxArray_real_T *varargin_2_d,
             const emxArray_int32_T *varargin_2_colidx, const emxArray_int32_T
             *varargin_2_rowidx, int varargin_2_m, int varargin_2_n)
{
  bool p;
  bool b_p;
  int nzx1;
  int k;
  int exitg3;
  int exitg2;
  bool exitg1;
  p = false;
  b_p = false;
  if ((varargin_1_m == varargin_2_m) && (varargin_1_n == varargin_2_n)) {
    b_p = true;
  }

  if (b_p) {
    nzx1 = varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 2;
    b_p = (varargin_1_colidx->data[varargin_1_colidx->size[0] - 1] - 1 ==
           varargin_2_colidx->data[varargin_2_colidx->size[0] - 1] - 1);
    if (b_p) {
      k = 0;
      do {
        exitg3 = 0;
        if (k <= varargin_1_n) {
          if (varargin_1_colidx->data[k] != varargin_2_colidx->data[k]) {
            b_p = false;
            exitg3 = 1;
          } else {
            k++;
          }
        } else {
          k = 0;
          exitg3 = 2;
        }
      } while (exitg3 == 0);

      if (exitg3 != 1) {
        do {
          exitg2 = 0;
          if (k <= nzx1) {
            if (varargin_1_rowidx->data[k] != varargin_2_rowidx->data[k]) {
              b_p = false;
              exitg2 = 1;
            } else {
              k++;
            }
          } else {
            k = 0;
            exitg2 = 2;
          }
        } while (exitg2 == 0);

        if (exitg2 != 1) {
          exitg1 = false;
          while ((!exitg1) && (k <= nzx1)) {
            if (!(varargin_1_d->data[k] == varargin_2_d->data[k])) {
              b_p = false;
              exitg1 = true;
            } else {
              k++;
            }
          }
        }
      }
    }
  }

  return b_p || p;
}

/* End of code generation (isequal.cpp) */
