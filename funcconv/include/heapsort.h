/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * heapsort.h
 *
 * Code generation for function 'heapsort'
 *
 */

#ifndef HEAPSORT_H
#define HEAPSORT_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void b_heapsort(emxArray_int32_T *x, int xstart, int xend, const
  cell_wrap_2 cmp_tunableEnvironment[2]);

#endif

/* End of code generation (heapsort.h) */
