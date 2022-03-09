/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spdiags.h
 *
 * Code generation for function 'spdiags'
 *
 */

#ifndef SPDIAGS_H
#define SPDIAGS_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void b_spdiags(const bool arg1[2401], double arg3, double arg4,
                      emxArray_boolean_T *res1_d, emxArray_int32_T *res1_colidx,
                      emxArray_int32_T *res1_rowidx, int *res1_m, int *res1_n);
extern void spdiags(const double arg1_data[], const int arg1_size[2], double
                    arg3, double arg4, coder_internal_sparse *res1);

#endif

/* End of code generation (spdiags.h) */
