/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse1.h
 *
 * Code generation for function 'sparse1'
 *
 */

#ifndef SPARSE1_H
#define SPARSE1_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void b_sparse_minus(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d, const
  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m, int b_n,
  emxArray_real_T *s_d, emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx,
  int *s_m, int *s_n);
extern void b_sparse_plus(const emxArray_real_T *a_d, const emxArray_int32_T
  *a_colidx, const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d, const
  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m, int b_n,
  emxArray_real_T *s_d, emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx,
  int *s_m, int *s_n);
extern void b_sparse_times(double a, const emxArray_boolean_T *b_d, const
  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m, int b_n,
  coder_internal_sparse *s);
extern void sparse_minus(const emxArray_real_T *a, const emxArray_real_T *b_d,
  const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m,
  int b_n, emxArray_real_T *s);
extern void sparse_mldivide(const emxArray_real_T *A_d, const emxArray_int32_T
  *A_colidx, const emxArray_int32_T *A_rowidx, int A_m, int A_n, const double
  b_data[], const int b_size[1], emxArray_real_T *y);
extern void sparse_plus(const emxArray_real_T *a, const emxArray_real_T *b_d,
  const emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m,
  int b_n, emxArray_real_T *s);
extern void sparse_times(double a, const emxArray_real_T *b_d, const
  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx, int b_m, int b_n,
  coder_internal_sparse *s);

#endif

/* End of code generation (sparse1.h) */
