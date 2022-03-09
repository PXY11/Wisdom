/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "power.h"
#include "funcconv.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void power(const double a_data[], const int a_size[1], double y_data[], int
           y_size[1])
{
  int nx;
  int k;
  y_size[0] = a_size[0];
  nx = a_size[0];
  for (k = 0; k < nx; k++) {
    y_data[k] = a_data[k] * a_data[k];
  }
}

/* End of code generation (power.cpp) */
