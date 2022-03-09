/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes1.cpp
 *
 * Code generation for function 'mtimes1'
 *
 */

/* Include files */
#include "mtimes1.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mtimes(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            C)
{
  int m;
  int inner;
  int i;
  int k;
  int aoffset;
  m = A->size[0];
  inner = A->size[1];
  i = C->size[0];
  C->size[0] = A->size[0];
  emxEnsureCapacity_real_T(C, i);
  for (i = 0; i < m; i++) {
    C->data[i] = 0.0;
  }

  for (k = 0; k < inner; k++) {
    aoffset = k * m;
    for (i = 0; i < m; i++) {
      C->data[i] += B->data[k] * A->data[aoffset + i];
    }
  }
}

/* End of code generation (mtimes1.cpp) */
