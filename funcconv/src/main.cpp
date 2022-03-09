/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.cpp
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "funcconv.h"
#include "funcconv_emxAPI.h"
#include "funcconv_terminate.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_1x12_real_T(double result[12]);
static void argInit_1x2_real_T(double result[2]);
static void argInit_1x6_real_T(double result[6]);
static void argInit_2401x1_real_T(double result[2401]);
static emxArray_real_T *argInit_Unboundedx1_real_T();
static double argInit_real_T();
static void main_funcconv();

/* Function Definitions */
static void argInit_1x12_real_T(double result[12])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 12; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

static void argInit_1x2_real_T(double result[2])
{
  double result_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result[0] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp;
}

static void argInit_1x6_real_T(double result[6])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 6; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

static void argInit_2401x1_real_T(double result[2401])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2401; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

static emxArray_real_T *argInit_Unboundedx1_real_T()
{
  emxArray_real_T *result;
  static const int iv[1] = { 2 };

  int idx0;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(1, iv);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[idx0] = argInit_real_T();
  }

  return result;
}

static double argInit_real_T()
{
  return 0.0;
}

static void main_funcconv()
{
  double T_tmp_tmp_tmp_tmp;
  double S0;
  double Smax;
  double sigma;
  double final_coupon_rate;
  double q;
  double p;
  double R;
  double eta;
  double no_call_time_tmp_tmp[2];
  double Bc_star;
  double Bp_star;
  double theta;
  double coupon_time_rate[12];
  double rate_T_tmp[6];
  double dtime;
  double risk_free_rate[6];
  double k;
  double ds;
  double dt;
  double rhopenltycall;
  double rhopenltyput;
  double i_tmp[2401];
  double Bc_T;
  double Bp_T;
  double k_T[2401];
  emxArray_real_T *u;
  emxArray_real_T *B;

  /* Initialize function 'funcconv' input arguments. */
  T_tmp_tmp_tmp_tmp = argInit_real_T();
  S0 = argInit_real_T();
  Smax = argInit_real_T();
  sigma = argInit_real_T();
  final_coupon_rate = argInit_real_T();
  q = argInit_real_T();
  p = argInit_real_T();
  R = argInit_real_T();
  eta = argInit_real_T();

  /* Initialize function input argument 'no_call_time'. */
  argInit_1x2_real_T(no_call_time_tmp_tmp);

  /* Initialize function input argument 'no_put_time'. */
  /* Initialize function input argument 'no_convert_time'. */
  Bc_star = argInit_real_T();
  Bp_star = argInit_real_T();
  theta = argInit_real_T();

  /* Initialize function input argument 'coupon_time_rate'. */
  argInit_1x12_real_T(coupon_time_rate);

  /* Initialize function input argument 'rate_T'. */
  argInit_1x6_real_T(rate_T_tmp);
  dtime = argInit_real_T();

  /* Initialize function input argument 'risk_free_rate'. */
  argInit_1x6_real_T(risk_free_rate);
  k = argInit_real_T();
  ds = argInit_real_T();
  dt = argInit_real_T();
  rhopenltycall = argInit_real_T();
  rhopenltyput = argInit_real_T();

  /* Initialize function input argument 'i'. */
  argInit_2401x1_real_T(i_tmp);

  /* Initialize function input argument 'S'. */
  Bc_T = argInit_real_T();
  Bp_T = argInit_real_T();

  /* Initialize function input argument 'k_T'. */
  argInit_2401x1_real_T(k_T);

  /* Initialize function input argument 'u'. */
  u = argInit_Unboundedx1_real_T();

  /* Initialize function input argument 'B'. */
  B = argInit_Unboundedx1_real_T();

  /* Call the entry-point 'funcconv'. */
  T_tmp_tmp_tmp_tmp = funcconv(T_tmp_tmp_tmp_tmp, T_tmp_tmp_tmp_tmp,
    T_tmp_tmp_tmp_tmp, T_tmp_tmp_tmp_tmp, T_tmp_tmp_tmp_tmp, S0, Smax, sigma,
    final_coupon_rate, q, p, R, eta, no_call_time_tmp_tmp, no_call_time_tmp_tmp,
    no_call_time_tmp_tmp, Bc_star, Bp_star, theta, coupon_time_rate, rate_T_tmp,
    dtime, risk_free_rate, k, ds, dt, rhopenltycall, rhopenltyput, i_tmp, i_tmp,
    Bc_T, Bp_T, k_T, u, B);
  emxDestroyArray_real_T(B);
  emxDestroyArray_real_T(u);
}

int main(int, const char * const [])
{
  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_funcconv();

  /* Terminate the application.
     You do not need to do this more than one time. */
  funcconv_terminate();
  return 0;
}

/* End of code generation (main.cpp) */
