/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * funcconv.h
 *
 * Code generation for function 'funcconv'
 *
 */

#ifndef FUNCCONV_H
#define FUNCCONV_H

/* Include files */
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "funcconv_types.h"

/* Function Declarations */
extern double funcconv(double T, double Nt, double Ns, double F, double
  conversion_price, double S0, double Smax, double sigma, double
  final_coupon_rate, double q, double p, double R, double eta, const double
  no_call_time[2], const double no_put_time[2], const double no_convert_time[2],
  double Bc_star, double Bp_star, double theta, const double coupon_time_rate[12],
  const double rate_T[6], double dtime, const double risk_free_rate[6], double k,
  double ds, double dt, double rhopenltycall, double rhopenltyput, const double
  i[2401], const double S[2401], double Bc_T, double Bp_T, const double k_T[2401],
  emxArray_real_T *u, emxArray_real_T *B);

#endif

/* End of code generation (funcconv.h) */
