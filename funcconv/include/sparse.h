/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.h
 *
 * Code generation for function 'sparse'
 *
 */

#ifndef SPARSE_H
#define SPARSE_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void sparse(const emxArray_real_T *varargin_1, emxArray_real_T *y_d,
                   emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int
                   *y_m, int *y_n);

#endif

/* End of code generation (sparse.h) */
