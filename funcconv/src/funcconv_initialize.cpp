/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * funcconv_initialize.cpp
 *
 * Code generation for function 'funcconv_initialize'
 *
 */

/* Include files */
#include "funcconv_initialize.h"
#include "funcconv.h"
#include "funcconv_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void funcconv_initialize()
{
  rt_InitInfAndNaN();
  omp_init_nest_lock(&emlrtNestLockGlobal);
  isInitialized_funcconv = true;
}

/* End of code generation (funcconv_initialize.cpp) */
