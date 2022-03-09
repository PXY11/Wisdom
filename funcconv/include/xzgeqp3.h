/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeqp3.h
 *
 * Code generation for function 'xzgeqp3'
 *
 */

#ifndef XZGEQP3_H
#define XZGEQP3_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void qrpf(emxArray_real_T *A, int m, int n, emxArray_real_T *tau,
                 emxArray_int32_T *jpvt);

#endif

/* End of code generation (xzgeqp3.h) */
