/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * introsort.h
 *
 * Code generation for function 'introsort'
 *
 */

#ifndef INTROSORT_H
#define INTROSORT_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void introsort(emxArray_int32_T *x, int xend, const cell_wrap_2
                      cmp_tunableEnvironment[2]);

#endif

/* End of code generation (introsort.h) */
