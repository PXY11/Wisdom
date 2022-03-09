/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.cpp
 *
 * Code generation for function 'sparse'
 *
 */

/* Include files */
#include "sparse.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void sparse(const emxArray_real_T *varargin_1, emxArray_real_T *y_d,
            emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int *y_m,
            int *y_n)
{
  int mInt;
  int nInt;
  int numalloc;
  int row;
  int ctr;
  double xrc;
  mInt = varargin_1->size[0];
  nInt = varargin_1->size[1];
  numalloc = 0;
  row = varargin_1->size[0] * varargin_1->size[1];
  for (ctr = 0; ctr < row; ctr++) {
    if (varargin_1->data[ctr] != 0.0) {
      numalloc++;
    }
  }

  *y_m = varargin_1->size[0];
  *y_n = varargin_1->size[1];
  if (numalloc < 1) {
    numalloc = 1;
  }

  row = y_d->size[0];
  y_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(y_d, row);
  for (row = 0; row < numalloc; row++) {
    y_d->data[row] = 0.0;
  }

  row = y_colidx->size[0];
  y_colidx->size[0] = varargin_1->size[1] + 1;
  emxEnsureCapacity_int32_T(y_colidx, row);
  ctr = varargin_1->size[1];
  for (row = 0; row <= ctr; row++) {
    y_colidx->data[row] = 0;
  }

  y_colidx->data[0] = 1;
  row = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(y_rowidx, row);
  for (row = 0; row < numalloc; row++) {
    y_rowidx->data[row] = 0;
  }

  y_rowidx->data[0] = 1;
  ctr = 0;
  for (numalloc = 0; numalloc < nInt; numalloc++) {
    for (row = 0; row < mInt; row++) {
      xrc = varargin_1->data[row + varargin_1->size[0] * numalloc];
      if (xrc != 0.0) {
        y_rowidx->data[ctr] = row + 1;
        y_d->data[ctr] = xrc;
        ctr++;
      }
    }

    y_colidx->data[numalloc + 1] = ctr + 1;
  }
}

/* End of code generation (sparse.cpp) */
