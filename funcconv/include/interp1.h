/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * interp1.h
 *
 * Code generation for function 'interp1'
 *
 */

#ifndef INTERP1_H
#define INTERP1_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern void b_interp1(const double varargin_1[2401], const emxArray_real_T
                      *varargin_2, const double varargin_3[2401], double Vq[2401]);
extern double interp1(const double varargin_1[6], const double varargin_2[6],
                      double varargin_3);

#endif

/* End of code generation (interp1.h) */
