/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eye.cpp
 *
 * Code generation for function 'eye'
 *
 */

/* Include files */
#include "eye.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void eye(double varargin_1, emxArray_real_T *b_I)
{
  double t;
  int m;
  int k;
  int loop_ub;
  if (varargin_1 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_1;
  }

  m = static_cast<int>(t);
  k = b_I->size[0] * b_I->size[1];
  b_I->size[0] = static_cast<int>(t);
  b_I->size[1] = static_cast<int>(t);
  emxEnsureCapacity_real_T(b_I, k);
  loop_ub = static_cast<int>(t) * static_cast<int>(t);
  for (k = 0; k < loop_ub; k++) {
    b_I->data[k] = 0.0;
  }

  if (static_cast<int>(t) > 0) {
    for (k = 0; k < m; k++) {
      b_I->data[k + b_I->size[0] * k] = 1.0;
    }
  }
}

/* End of code generation (eye.cpp) */
