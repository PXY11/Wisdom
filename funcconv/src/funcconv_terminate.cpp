/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * funcconv_terminate.cpp
 *
 * Code generation for function 'funcconv_terminate'
 *
 */

/* Include files */
#include "funcconv_terminate.h"
#include "funcconv.h"
#include "funcconv_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void funcconv_terminate()
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
  isInitialized_funcconv = false;
}

/* End of code generation (funcconv_terminate.cpp) */
