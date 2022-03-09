/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * funcconv.cpp
 *
 * Code generation for function 'funcconv'
 *
 */

/* Include files */
#include "funcconv.h"
#include "eye.h"
#include "funcconv_data.h"
#include "funcconv_emxutil.h"
#include "funcconv_initialize.h"
#include "interp1.h"
#include "isequal.h"
#include "mldivide.h"
#include "mtimes.h"
#include "mtimes1.h"
#include "power.h"
#include "rt_nonfinite.h"
#include "sparse1.h"
#include "spdiags.h"
#include "speye.h"
#include <cmath>
#include <cstring>

/* Function Definitions */
double funcconv(double T, double Nt, double Ns, double F, double
                conversion_price, double S0, double, double sigma, double,
                double q, double p, double R, double eta, const double
                no_call_time[2], const double no_put_time[2], const double
                no_convert_time[2], double Bc_star, double Bp_star, double theta,
                const double coupon_time_rate[12], const double rate_T[6],
                double dtime, const double risk_free_rate[6], double, double ds,
                double dt, double rhopenltycall, double rhopenltyput, const
                double [2401], const double S[2401], double, double, const
                double [2401], emxArray_real_T *u, emxArray_real_T *B)
{
  double U;
  double count;
  int flag;
  int n;
  emxArray_real_T *Mu_d;
  emxArray_int32_T *Mu_colidx;
  emxArray_int32_T *Mu_rowidx;
  emxArray_real_T *MB_d;
  emxArray_int32_T *MB_colidx;
  emxArray_int32_T *MB_rowidx;
  emxArray_real_T *u_n;
  emxArray_real_T *P1_old_d;
  emxArray_int32_T *P1_old_colidx;
  emxArray_int32_T *P1_old_rowidx;
  emxArray_real_T *P2_old_d;
  emxArray_int32_T *P2_old_colidx;
  emxArray_int32_T *P2_old_rowidx;
  emxArray_real_T *u_old;
  emxArray_real_T *P1_d;
  emxArray_int32_T *P1_colidx;
  emxArray_int32_T *P1_rowidx;
  emxArray_real_T *P2_d;
  emxArray_int32_T *P2_colidx;
  emxArray_int32_T *P2_rowidx;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *c_d;
  emxArray_int32_T *c_colidx;
  emxArray_int32_T *c_rowidx;
  emxArray_real_T *b_c_d;
  emxArray_int32_T *b_c_colidx;
  emxArray_int32_T *b_c_rowidx;
  emxArray_real_T *y;
  emxArray_real_T *a;
  emxArray_boolean_T *b_d;
  emxArray_real_T *maxval;
  coder_internal_sparse expl_temp;
  emxArray_real_T *r2;
  double t_tmp;
  double b_r;
  double Bp_n;
  int b_i;
  double k_n[2401];
  double beta_data[2401];
  int c_i;
  int b_n;
  double a_tmp;
  double AccI;
  int loop_ub;
  int S_size[1];
  double alpha_data[2401];
  int alpha_size[1];
  double d;
  double d1;
  int beta_size[1];
  double alpha_MuMB_data[2403];
  double beta_MuMB_data[2403];
  int alpha_MuMB_size[2];
  static double b_alpha_MuMB_data[7209];
  int Mu_m;
  int Mu_n;
  int b_k;
  double Bc_n;
  double b_Bp_n;
  int c_m;
  int c_n;
  int tmp_size[1];
  double tmp_data[2402];
  int x_tmp;
  bool z1[2401];
  int P1_old_m;
  int P1_old_n;
  int P2_old_m;
  int P2_old_n;
  int exitg1;
  int P1_m;
  int P1_n;
  int P2_m;
  int P2_n;
  bool exitg2;
  if (isInitialized_funcconv == false) {
    funcconv_initialize();
  }

  count = 0.0;

  /* total iteration time */
  flag = 0;

  /* if not convergenc, break will stop the whole loop */
  n = 0;
  emxInit_real_T(&Mu_d, 1);
  emxInit_int32_T(&Mu_colidx, 1);
  emxInit_int32_T(&Mu_rowidx, 1);
  emxInit_real_T(&MB_d, 1);
  emxInit_int32_T(&MB_colidx, 1);
  emxInit_int32_T(&MB_rowidx, 1);
  emxInit_real_T(&u_n, 1);
  emxInit_real_T(&P1_old_d, 1);
  emxInit_int32_T(&P1_old_colidx, 1);
  emxInit_int32_T(&P1_old_rowidx, 1);
  emxInit_real_T(&P2_old_d, 1);
  emxInit_int32_T(&P2_old_colidx, 1);
  emxInit_int32_T(&P2_old_rowidx, 1);
  emxInit_real_T(&u_old, 1);
  emxInit_real_T(&P1_d, 1);
  emxInit_int32_T(&P1_colidx, 1);
  emxInit_int32_T(&P1_rowidx, 1);
  emxInit_real_T(&P2_d, 1);
  emxInit_int32_T(&P2_colidx, 1);
  emxInit_int32_T(&P2_rowidx, 1);
  emxInit_real_T(&r, 1);
  emxInit_real_T(&r1, 1);
  emxInit_real_T(&c_d, 1);
  emxInit_int32_T(&c_colidx, 1);
  emxInit_int32_T(&c_rowidx, 1);
  emxInit_real_T(&b_c_d, 1);
  emxInit_int32_T(&b_c_colidx, 1);
  emxInit_int32_T(&b_c_rowidx, 1);
  emxInit_real_T(&y, 1);
  emxInit_real_T(&a, 2);
  emxInit_boolean_T(&b_d, 1);
  emxInit_real_T(&maxval, 1);
  c_emxInitStruct_coder_internal_(&expl_temp);
  emxInit_real_T(&r2, 2);
  while ((n <= static_cast<int>(Nt) - 1) && (flag != 1)) {
    /* t=T-n*dt,n=1 means t=T-dt,n=N tmeans t=0 %want to u^n→u^(n+1) */
    /* computation using r(t) */
    t_tmp = (static_cast<double>(n) + 1.0) * dt;

    /* time to maturaty */
    /*  discount rate for time step n */
    /*  discount rate for time step n+1 */
    b_r = -(std::log(interp1(rate_T, risk_free_rate, ((static_cast<double>(n) +
                1.0) - 1.0) * dt - 2.2204460492503131E-16) / interp1(rate_T,
              risk_free_rate, t_tmp - 2.2204460492503131E-16)) / dt);

    /* risk free rate in time step n */
    /* r=0.03; */
    Bp_n = std::ceil(t_tmp - dtime);
    if (Bp_n - std::ceil((t_tmp - dt) - dtime) == 1.0) {
      /* u dividend adjustment */
      for (b_i = 0; b_i < 2401; b_i++) {
        k_n[b_i] = (1.0 - q) * S[b_i];
      }

      b_interp1(S, B, k_n, beta_data);
      c_i = B->size[0];
      B->size[0] = 2401;
      emxEnsureCapacity_real_T(B, c_i);
      for (c_i = 0; c_i < 2401; c_i++) {
        B->data[c_i] = beta_data[c_i];
      }

      /* B数组在C++中要循环赋值 2401 x 1 */
      b_interp1(S, u, k_n, beta_data);
      c_i = u->size[0];
      u->size[0] = 2401;
      emxEnsureCapacity_real_T(u, c_i);
      for (c_i = 0; c_i < 2401; c_i++) {
        u->data[c_i] = beta_data[c_i];
      }

      /* u数组在C++中要循环赋值 2401 x 1 */
    }

    /* S=(0:Ns)'*ds+q.*S.*ceil(t-dtime);%dividend influence to stock price */
    /* alpha 2399 x 1     */
    if (2.0 > Ns) {
      c_i = 0;
      b_i = 0;
      b_n = 0;
    } else {
      c_i = 1;
      b_i = static_cast<int>(Ns);
      b_n = 1;
    }

    a_tmp = sigma * sigma;
    AccI = (b_r + p * eta) - q;
    loop_ub = b_i - c_i;
    S_size[0] = loop_ub;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      k_n[b_i] = S[c_i + b_i];
    }

    power(k_n, S_size, alpha_data, alpha_size);
    d = 2.0 * (ds * ds);
    d1 = 2.0 * ds;
    loop_ub = alpha_size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      alpha_data[c_i] = (a_tmp * alpha_data[c_i] / d - AccI * S[b_n + c_i] / d1)
        * dt;
    }

    /* Ns-1, namely alpha_1~alpha_(Ns-1) means ds~(Ns-1)*ds */
    /* beta 2399 x 1 */
    if (2.0 > Ns) {
      c_i = 0;
      b_i = 0;
      b_n = 0;
    } else {
      c_i = 1;
      b_i = static_cast<int>(Ns);
      b_n = 1;
    }

    loop_ub = b_i - c_i;
    S_size[0] = loop_ub;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      k_n[b_i] = S[c_i + b_i];
    }

    power(k_n, S_size, beta_data, beta_size);
    d1 = 2.0 * ds;
    loop_ub = beta_size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      beta_data[c_i] = (a_tmp * beta_data[c_i] / d + AccI * S[b_n + c_i] / d1) *
        dt;
    }

    /* alpha_MuMB 2401 x 1 */
    b_n = alpha_size[0] + 2;
    loop_ub = alpha_size[0];
    if (0 <= loop_ub - 1) {
      std::memcpy(&alpha_MuMB_data[0], &alpha_data[0], loop_ub * sizeof(double));
    }

    alpha_MuMB_data[alpha_size[0]] = 0.0;
    alpha_MuMB_data[alpha_size[0] + 1] = 0.0;

    /* take Mu and MB down diagonals as alpha_MuMB */
    /* beta_MuMB 2401 x 1 */
    b_i = beta_size[0] + 2;
    beta_MuMB_data[0] = 0.0;
    beta_MuMB_data[1] = 0.0;
    loop_ub = beta_size[0];
    if (0 <= loop_ub - 1) {
      std::memcpy(&beta_MuMB_data[2], &beta_data[0], loop_ub * sizeof(double));
    }

    /* take Mu and MB up diagonals as beta_MuMB */
    /* gamma_Mu 2401 x 1 */
    /* take Mu diagonals as gamma_Mu */
    /* Mu 2401 x 2401 */
    AccI = b_r + p;
    a_tmp = AccI * dt;
    alpha_MuMB_size[0] = b_n;
    alpha_MuMB_size[1] = 3;
    if (0 <= b_n - 1) {
      std::memcpy(&b_alpha_MuMB_data[0], &alpha_MuMB_data[0], b_n * sizeof
                  (double));
    }

    b_alpha_MuMB_data[b_n] = -AccI * dt;
    loop_ub = alpha_size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      b_alpha_MuMB_data[(c_i + b_n) + 1] = -((alpha_data[c_i] + beta_data[c_i])
        + a_tmp);
    }

    b_alpha_MuMB_data[(alpha_size[0] + b_n) + 1] = 0.0;
    for (c_i = 0; c_i < b_i; c_i++) {
      b_alpha_MuMB_data[c_i + b_n * 2] = beta_MuMB_data[c_i];
    }

    spdiags(b_alpha_MuMB_data, alpha_MuMB_size, Ns + 1.0, Ns + 1.0, &expl_temp);
    c_i = Mu_d->size[0];
    Mu_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(Mu_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      Mu_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = Mu_colidx->size[0];
    Mu_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(Mu_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      Mu_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = Mu_rowidx->size[0];
    Mu_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(Mu_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      Mu_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    Mu_m = expl_temp.m;
    Mu_n = expl_temp.n;

    /* this is Mu Mu的定义在这里！！！ */
    /* gamma_MB 2401 x 1 */
    /* take MB diagonals as gamma_MB */
    /* MB 2401 x 2401 */
    AccI = b_r + p * (1.0 - R);
    b_r = AccI * dt;
    alpha_MuMB_size[0] = b_n;
    alpha_MuMB_size[1] = 3;
    if (0 <= b_n - 1) {
      std::memcpy(&b_alpha_MuMB_data[0], &alpha_MuMB_data[0], b_n * sizeof
                  (double));
    }

    b_alpha_MuMB_data[b_n] = -AccI * dt;
    loop_ub = alpha_size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      b_alpha_MuMB_data[(c_i + b_n) + 1] = -((alpha_data[c_i] + beta_data[c_i])
        + b_r);
    }

    b_alpha_MuMB_data[(alpha_size[0] + b_n) + 1] = 0.0;
    for (c_i = 0; c_i < b_i; c_i++) {
      b_alpha_MuMB_data[c_i + b_n * 2] = beta_MuMB_data[c_i];
    }

    spdiags(b_alpha_MuMB_data, alpha_MuMB_size, Ns + 1.0, Ns + 1.0, &expl_temp);
    c_i = MB_d->size[0];
    MB_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(MB_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      MB_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = MB_colidx->size[0];
    MB_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(MB_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      MB_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = MB_rowidx->size[0];
    MB_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(MB_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      MB_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    b_n = expl_temp.m;
    b_k = expl_temp.n;

    /* This is MB */
    /* spdiags用法 知乎 https://zhuanlan.zhihu.com/p/365797391 */
    b_r = std::ceil(t_tmp);
    AccI = (b_r - t_tmp) * coupon_time_rate[static_cast<int>(b_r) - 1] * F;

    /* counpond rate, ACCI */
    if ((no_call_time[0] <= t_tmp) && (t_tmp <= no_call_time[1])) {
      Bc_n = 1000.0;

      /* when n value of Bc,take a large number means no call */
    } else {
      Bc_n = Bc_star + AccI;

      /* dirty price  */
    }

    if ((no_put_time[0] <= t_tmp) && (t_tmp <= no_put_time[1])) {
      b_Bp_n = 0.01;

      /* small number means no put */
    } else {
      b_Bp_n = Bp_star + AccI;
    }

    if ((no_convert_time[0] <= t_tmp) && (t_tmp <= no_convert_time[1])) {
      for (b_i = 0; b_i < 2401; b_i++) {
        k_n[b_i] = 0.01;
      }

      /* small means no converte */
    } else {
      a_tmp = std::ceil(T - dtime);
      for (b_i = 0; b_i < 2401; b_i++) {
        k_n[b_i] = F / (conversion_price - q * S[b_i] * (a_tmp - Bp_n));
      }

      /* when T=X.5, stock pay dividends */
    }

    /* 测试减少缩进 */
    c_i = u_n->size[0];
    u_n->size[0] = u->size[0];
    emxEnsureCapacity_real_T(u_n, c_i);
    loop_ub = u->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      u_n->data[c_i] = u->data[c_i];
    }

    /* now value of u is n-1 step, will converte to n step following  */
    /*      coder.varsize('B'); */
    sparse_times(theta, MB_d, MB_colidx, MB_rowidx, expl_temp.m, expl_temp.n,
                 &expl_temp);
    c_i = c_d->size[0];
    c_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(c_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      c_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = c_colidx->size[0];
    c_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(c_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = c_rowidx->size[0];
    c_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(c_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    c_m = expl_temp.m;
    c_n = expl_temp.n;
    sparse_times(1.0 - theta, MB_d, MB_colidx, MB_rowidx, b_n, b_k, &expl_temp);
    c_i = b_c_d->size[0];
    b_c_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(b_c_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      b_c_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = b_c_colidx->size[0];
    b_c_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(b_c_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      b_c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = b_c_rowidx->size[0];
    b_c_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(b_c_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      b_c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    eye(Ns + 1.0, r2);
    sparse_plus(r2, b_c_d, b_c_colidx, b_c_rowidx, expl_temp.m, expl_temp.n, a);
    if ((a->size[1] == 1) || (B->size[0] == 1)) {
      c_i = y->size[0];
      y->size[0] = a->size[0];
      emxEnsureCapacity_real_T(y, c_i);
      loop_ub = a->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        y->data[c_i] = 0.0;
        b_n = a->size[1];
        for (b_i = 0; b_i < b_n; b_i++) {
          y->data[c_i] += a->data[c_i + a->size[0] * b_i] * B->data[b_i];
        }
      }
    } else {
      mtimes(a, B, y);
    }

    c_i = B->size[0];
    B->size[0] = y->size[0];
    emxEnsureCapacity_real_T(B, c_i);
    loop_ub = y->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      B->data[c_i] = y->data[c_i];
    }

    eye(Ns + 1.0, r2);
    sparse_minus(r2, c_d, c_colidx, c_rowidx, c_m, c_n, a);
    mldiv(a, B);

    /* similar to eqation (3.45)，without constratint, B's eqation, */
    /*  solution B_n is n step iteration initial value */
    b_n = B->size[0];
    for (b_i = 0; b_i < b_n; b_i++) {
      if (B->data[b_i] > Bc_n) {
        B->data[b_i] = Bc_n;
      }
    }

    /* →explict B<=Bc */
    /* 20220115改写到这 */
    c_i = y->size[0];
    y->size[0] = B->size[0];
    emxEnsureCapacity_real_T(y, c_i);
    loop_ub = B->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      y->data[c_i] = R * B->data[c_i];
    }

    a_tmp = p * dt;
    for (b_k = 0; b_k < 2401; b_k++) {
      Bp_n = k_n[b_k] * (1.0 - eta) * S[b_k];
      if ((!(Bp_n > y->data[b_k])) && (!rtIsNaN(y->data[b_k]))) {
        Bp_n = y->data[b_k];
      }

      Bp_n *= a_tmp;
      alpha_data[b_k] = Bp_n;
    }

    /* similar to (3.44) the right term */
    if (1.0 > Ns) {
      loop_ub = -1;
    } else {
      loop_ub = static_cast<int>(Ns) - 1;
    }

    sparse_times(1.0 - theta, Mu_d, Mu_colidx, Mu_rowidx, Mu_m, Mu_n, &expl_temp);
    c_i = c_d->size[0];
    c_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(c_d, c_i);
    b_n = expl_temp.d->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = c_colidx->size[0];
    c_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(c_colidx, c_i);
    b_n = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = c_rowidx->size[0];
    c_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(c_rowidx, c_i);
    b_n = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    c_m = expl_temp.m;
    c_n = expl_temp.n;
    speye(Ns + 1.0, &expl_temp);
    c_i = MB_d->size[0];
    MB_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(MB_d, c_i);
    b_n = expl_temp.d->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = MB_colidx->size[0];
    MB_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(MB_colidx, c_i);
    b_n = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = MB_rowidx->size[0];
    MB_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(MB_rowidx, c_i);
    b_n = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    b_sparse_plus(MB_d, MB_colidx, MB_rowidx, c_d, c_colidx, c_rowidx, c_m, c_n,
                  b_c_d, b_c_colidx, b_c_rowidx, &b_i, &b_n);
    sparse_mtimes(b_c_d, b_c_colidx, b_c_rowidx, b_i, b_n, u, r);
    sparse_times(theta, Mu_d, Mu_colidx, Mu_rowidx, Mu_m, Mu_n, &expl_temp);
    c_i = c_d->size[0];
    c_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(c_d, c_i);
    b_n = expl_temp.d->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = c_colidx->size[0];
    c_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(c_colidx, c_i);
    b_n = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = c_rowidx->size[0];
    c_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(c_rowidx, c_i);
    b_n = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    c_m = expl_temp.m;
    c_n = expl_temp.n;
    speye(Ns + 1.0, &expl_temp);
    c_i = MB_d->size[0];
    MB_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(MB_d, c_i);
    b_n = expl_temp.d->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = MB_colidx->size[0];
    MB_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(MB_colidx, c_i);
    b_n = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = MB_rowidx->size[0];
    MB_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(MB_rowidx, c_i);
    b_n = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < b_n; c_i++) {
      MB_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    b_sparse_minus(MB_d, MB_colidx, MB_rowidx, c_d, c_colidx, c_rowidx, c_m, c_n,
                   b_c_d, b_c_colidx, b_c_rowidx, &b_i, &b_n);
    tmp_size[0] = loop_ub + 2;
    for (c_i = 0; c_i <= loop_ub; c_i++) {
      tmp_data[c_i] = r->data[c_i] + alpha_data[c_i];
    }

    tmp_data[loop_ub + 1] = r->data[loop_ub + 1];
    sparse_mldivide(b_c_d, b_c_colidx, b_c_rowidx, b_i, b_n, tmp_data, tmp_size,
                    u);

    /* similar to eqation (3.44)，without constratint, u's eqation, solution u_n is n step iteration initial value */
    a_tmp = t_tmp + dt;
    b_r = std::ceil(a_tmp) - b_r;
    x_tmp = static_cast<int>(std::ceil(a_tmp + 1.0E-6)) - 1;
    a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
    AccI = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
    for (b_k = 0; b_k < 2401; b_k++) {
      Bp_n = k_n[b_k] * S[b_k];
      k_n[b_k] = Bp_n;
      d = Bp_n + a_tmp;
      if ((Bc_n > d) || rtIsNaN(d)) {
        d = Bc_n;
      }

      if ((u->data[b_k] < d) || rtIsNaN(d)) {
        d1 = u->data[b_k];
      } else {
        d1 = d;
      }

      d = Bp_n + AccI;
      if ((b_Bp_n > d) || rtIsNaN(d)) {
        d = b_Bp_n;
      }

      if ((d1 > d) || rtIsNaN(d)) {
        beta_data[b_k] = d1;
      } else {
        beta_data[b_k] = d;
      }
    }

    c_i = u->size[0];
    u->size[0] = 2401;
    emxEnsureCapacity_real_T(u, c_i);

    /* →explict constraint u, initial value of u_n */
    a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
    for (b_k = 0; b_k < 2401; b_k++) {
      u->data[b_k] = beta_data[b_k];
      Bp_n = k_n[b_k] + a_tmp;
      if ((b_Bp_n > Bp_n) || rtIsNaN(Bp_n)) {
        Bp_n = b_Bp_n;
      }

      z1[b_k] = (beta_data[b_k] <= Bp_n);
    }

    b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
    b_sparse_times(rhopenltyput, b_d, MB_colidx, MB_rowidx, b_i, b_n, &expl_temp);
    c_i = P1_old_d->size[0];
    P1_old_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(P1_old_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P1_old_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = P1_old_colidx->size[0];
    P1_old_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(P1_old_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P1_old_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = P1_old_rowidx->size[0];
    P1_old_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(P1_old_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P1_old_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    P1_old_m = expl_temp.m;
    P1_old_n = expl_temp.n;
    a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
    for (b_k = 0; b_k < 2401; b_k++) {
      Bp_n = k_n[b_k] + a_tmp;
      if ((Bc_n > Bp_n) || rtIsNaN(Bp_n)) {
        Bp_n = Bc_n;
      }

      z1[b_k] = (beta_data[b_k] >= Bp_n);
    }

    b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
    b_sparse_times(-rhopenltycall, b_d, MB_colidx, MB_rowidx, b_i, b_n,
                   &expl_temp);
    c_i = P2_old_d->size[0];
    P2_old_d->size[0] = expl_temp.d->size[0];
    emxEnsureCapacity_real_T(P2_old_d, c_i);
    loop_ub = expl_temp.d->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P2_old_d->data[c_i] = expl_temp.d->data[c_i];
    }

    c_i = P2_old_colidx->size[0];
    P2_old_colidx->size[0] = expl_temp.colidx->size[0];
    emxEnsureCapacity_int32_T(P2_old_colidx, c_i);
    loop_ub = expl_temp.colidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P2_old_colidx->data[c_i] = expl_temp.colidx->data[c_i];
    }

    c_i = P2_old_rowidx->size[0];
    P2_old_rowidx->size[0] = expl_temp.rowidx->size[0];
    emxEnsureCapacity_int32_T(P2_old_rowidx, c_i);
    loop_ub = expl_temp.rowidx->size[0];
    for (c_i = 0; c_i < loop_ub; c_i++) {
      P2_old_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
    }

    P2_old_m = expl_temp.m;
    P2_old_n = expl_temp.n;

    /*  P1,P2 penalty parameter initial value, when this item =0, penalty parameter is 0,  */
    /*  otherwise is rho, P1_n^0 is max(Bp_n,kS),P2_n^0 is max(Bc_n,kS) */
    /* 测试减少缩进 */
    do {
      exitg1 = 0;
      count++;
      c_i = u_old->size[0];
      u_old->size[0] = u->size[0];
      emxEnsureCapacity_real_T(u_old, c_i);
      loop_ub = u->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        u_old->data[c_i] = u->data[c_i];
      }

      /* u_old is u^0_n,i.e.,n step initial value */
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((b_Bp_n > Bp_n) || rtIsNaN(Bp_n)) {
          Bp_n = b_Bp_n;
        }

        z1[b_k] = (u->data[b_k] <= Bp_n);
      }

      b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
      b_sparse_times(rhopenltyput, b_d, MB_colidx, MB_rowidx, b_i, b_n,
                     &expl_temp);
      c_i = P1_d->size[0];
      P1_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(P1_d, c_i);
      loop_ub = expl_temp.d->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = P1_colidx->size[0];
      P1_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(P1_colidx, c_i);
      loop_ub = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = P1_rowidx->size[0];
      P1_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(P1_rowidx, c_i);
      loop_ub = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      P1_m = expl_temp.m;
      P1_n = expl_temp.n;
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((Bc_n > Bp_n) || rtIsNaN(Bp_n)) {
          Bp_n = Bc_n;
        }

        z1[b_k] = (u->data[b_k] >= Bp_n);
      }

      b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
      b_sparse_times(-rhopenltycall, b_d, MB_colidx, MB_rowidx, b_i, b_n,
                     &expl_temp);
      c_i = P2_d->size[0];
      P2_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(P2_d, c_i);
      loop_ub = expl_temp.d->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = P2_colidx->size[0];
      P2_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(P2_colidx, c_i);
      loop_ub = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = P2_rowidx->size[0];
      P2_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(P2_rowidx, c_i);
      loop_ub = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      P2_m = expl_temp.m;
      P2_n = expl_temp.n;
      if (1.0 > Ns) {
        loop_ub = -1;
      } else {
        loop_ub = static_cast<int>(Ns) - 1;
      }

      sparse_times(1.0 - theta, Mu_d, Mu_colidx, Mu_rowidx, Mu_m, Mu_n,
                   &expl_temp);
      c_i = c_d->size[0];
      c_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(c_d, c_i);
      b_n = expl_temp.d->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = c_colidx->size[0];
      c_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(c_colidx, c_i);
      b_n = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = c_rowidx->size[0];
      c_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(c_rowidx, c_i);
      b_n = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      c_m = expl_temp.m;
      c_n = expl_temp.n;
      speye(Ns + 1.0, &expl_temp);
      c_i = MB_d->size[0];
      MB_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(MB_d, c_i);
      b_n = expl_temp.d->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = MB_colidx->size[0];
      MB_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(MB_colidx, c_i);
      b_n = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = MB_rowidx->size[0];
      MB_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(MB_rowidx, c_i);
      b_n = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      b_sparse_plus(MB_d, MB_colidx, MB_rowidx, c_d, c_colidx, c_rowidx, c_m,
                    c_n, b_c_d, b_c_colidx, b_c_rowidx, &b_i, &b_n);
      sparse_mtimes(b_c_d, b_c_colidx, b_c_rowidx, b_i, b_n, u_n, r);
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((b_Bp_n > Bp_n) || rtIsNaN(Bp_n)) {
          beta_data[b_k] = b_Bp_n;
        } else {
          beta_data[b_k] = Bp_n;
        }
      }

      b_sparse_mtimes(P1_d, P1_colidx, P1_rowidx, P1_m, P1_n, beta_data, y);
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((Bc_n > Bp_n) || rtIsNaN(Bp_n)) {
          beta_data[b_k] = Bc_n;
        } else {
          beta_data[b_k] = Bp_n;
        }
      }

      b_sparse_mtimes(P2_d, P2_colidx, P2_rowidx, P2_m, P2_n, beta_data, r1);
      sparse_times(theta, Mu_d, Mu_colidx, Mu_rowidx, Mu_m, Mu_n, &expl_temp);
      c_i = c_d->size[0];
      c_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(c_d, c_i);
      b_n = expl_temp.d->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = c_colidx->size[0];
      c_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(c_colidx, c_i);
      b_n = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = c_rowidx->size[0];
      c_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(c_rowidx, c_i);
      b_n = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        c_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      c_m = expl_temp.m;
      c_n = expl_temp.n;
      speye(Ns + 1.0, &expl_temp);
      c_i = MB_d->size[0];
      MB_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(MB_d, c_i);
      b_n = expl_temp.d->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = MB_colidx->size[0];
      MB_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(MB_colidx, c_i);
      b_n = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = MB_rowidx->size[0];
      MB_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(MB_rowidx, c_i);
      b_n = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < b_n; c_i++) {
        MB_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      b_sparse_minus(MB_d, MB_colidx, MB_rowidx, c_d, c_colidx, c_rowidx, c_m,
                     c_n, b_c_d, b_c_colidx, b_c_rowidx, &b_i, &b_n);
      b_sparse_plus(b_c_d, b_c_colidx, b_c_rowidx, P1_d, P1_colidx, P1_rowidx,
                    P1_m, P1_n, MB_d, MB_colidx, MB_rowidx, &b_n, &b_k);
      b_sparse_minus(MB_d, MB_colidx, MB_rowidx, P2_d, P2_colidx, P2_rowidx,
                     P2_m, P2_n, b_c_d, b_c_colidx, b_c_rowidx, &b_i, &b_n);
      tmp_size[0] = loop_ub + 2;
      for (c_i = 0; c_i <= loop_ub; c_i++) {
        tmp_data[c_i] = ((r->data[c_i] + alpha_data[c_i]) + y->data[c_i]) -
          r1->data[c_i];
      }

      tmp_data[loop_ub + 1] = (r->data[loop_ub + 1] + y->data[loop_ub + 1]) -
        r1->data[loop_ub + 1];
      sparse_mldivide(b_c_d, b_c_colidx, b_c_rowidx, b_i, b_n, tmp_data,
                      tmp_size, u);

      /* similar to (4.28)，solution is u_n^(k+1) */
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((b_Bp_n > Bp_n) || rtIsNaN(Bp_n)) {
          Bp_n = b_Bp_n;
        }

        z1[b_k] = (u->data[b_k] <= Bp_n);
      }

      b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
      b_sparse_times(rhopenltyput, b_d, MB_colidx, MB_rowidx, b_i, b_n,
                     &expl_temp);
      c_i = P1_d->size[0];
      P1_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(P1_d, c_i);
      loop_ub = expl_temp.d->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = P1_colidx->size[0];
      P1_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(P1_colidx, c_i);
      loop_ub = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = P1_rowidx->size[0];
      P1_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(P1_rowidx, c_i);
      loop_ub = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P1_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      P1_m = expl_temp.m;
      P1_n = expl_temp.n;
      a_tmp = static_cast<double>((b_r == 1.0)) * coupon_time_rate[x_tmp] * F;
      for (b_k = 0; b_k < 2401; b_k++) {
        Bp_n = k_n[b_k] + a_tmp;
        if ((Bc_n > Bp_n) || rtIsNaN(Bp_n)) {
          Bp_n = Bc_n;
        }

        z1[b_k] = (u->data[b_k] >= Bp_n);
      }

      b_spdiags(z1, Ns + 1.0, Ns + 1.0, b_d, MB_colidx, MB_rowidx, &b_i, &b_n);
      b_sparse_times(-rhopenltycall, b_d, MB_colidx, MB_rowidx, b_i, b_n,
                     &expl_temp);
      c_i = P2_d->size[0];
      P2_d->size[0] = expl_temp.d->size[0];
      emxEnsureCapacity_real_T(P2_d, c_i);
      loop_ub = expl_temp.d->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_d->data[c_i] = expl_temp.d->data[c_i];
      }

      c_i = P2_colidx->size[0];
      P2_colidx->size[0] = expl_temp.colidx->size[0];
      emxEnsureCapacity_int32_T(P2_colidx, c_i);
      loop_ub = expl_temp.colidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_colidx->data[c_i] = expl_temp.colidx->data[c_i];
      }

      c_i = P2_rowidx->size[0];
      P2_rowidx->size[0] = expl_temp.rowidx->size[0];
      emxEnsureCapacity_int32_T(P2_rowidx, c_i);
      loop_ub = expl_temp.rowidx->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        P2_rowidx->data[c_i] = expl_temp.rowidx->data[c_i];
      }

      if (isequal(P1_d, P1_colidx, P1_rowidx, P1_m, P1_n, P1_old_d,
                  P1_old_colidx, P1_old_rowidx, P1_old_m, P1_old_n) && isequal
          (P2_d, P2_colidx, P2_rowidx, expl_temp.m, expl_temp.n, P2_old_d,
           P2_old_colidx, P2_old_rowidx, P2_old_m, P2_old_n)) {
        exitg1 = 1;
      } else {
        b_i = u->size[0];
        c_i = y->size[0];
        y->size[0] = u->size[0];
        emxEnsureCapacity_real_T(y, c_i);
        for (b_k = 0; b_k < b_i; b_k++) {
          y->data[b_k] = std::abs(u->data[b_k]);
        }

        c_i = maxval->size[0];
        maxval->size[0] = y->size[0];
        emxEnsureCapacity_real_T(maxval, c_i);
        b_i = y->size[0];
        for (b_k = 0; b_k < b_i; b_k++) {
          if ((1.0 > y->data[b_k]) || rtIsNaN(y->data[b_k])) {
            maxval->data[b_k] = 1.0;
          } else {
            maxval->data[b_k] = y->data[b_k];
          }
        }

        c_i = u_old->size[0];
        u_old->size[0] = u->size[0];
        emxEnsureCapacity_real_T(u_old, c_i);
        loop_ub = u->size[0];
        for (c_i = 0; c_i < loop_ub; c_i++) {
          u_old->data[c_i] = u->data[c_i] - u_old->data[c_i];
        }

        b_i = u_old->size[0];
        c_i = y->size[0];
        y->size[0] = u_old->size[0];
        emxEnsureCapacity_real_T(y, c_i);
        for (b_k = 0; b_k < b_i; b_k++) {
          y->data[b_k] = std::abs(u_old->data[b_k]);
        }

        loop_ub = y->size[0];
        for (c_i = 0; c_i < loop_ub; c_i++) {
          y->data[c_i] /= maxval->data[c_i];
        }

        b_n = y->size[0];
        if (y->size[0] <= 2) {
          if (y->size[0] == 1) {
            a_tmp = y->data[0];
          } else if ((y->data[0] < y->data[1]) || (rtIsNaN(y->data[0]) &&
                      (!rtIsNaN(y->data[1])))) {
            a_tmp = y->data[1];
          } else {
            a_tmp = y->data[0];
          }
        } else {
          if (!rtIsNaN(y->data[0])) {
            b_i = 1;
          } else {
            b_i = 0;
            b_k = 2;
            exitg2 = false;
            while ((!exitg2) && (b_k <= y->size[0])) {
              if (!rtIsNaN(y->data[b_k - 1])) {
                b_i = b_k;
                exitg2 = true;
              } else {
                b_k++;
              }
            }
          }

          if (b_i == 0) {
            a_tmp = y->data[0];
          } else {
            a_tmp = y->data[b_i - 1];
            c_i = b_i + 1;
            for (b_k = c_i; b_k <= b_n; b_k++) {
              Bp_n = y->data[b_k - 1];
              if (a_tmp < Bp_n) {
                a_tmp = Bp_n;
              }
            }
          }
        }

        if (a_tmp < 10000.0) {
          exitg1 = 1;
        } else if (count > Nt * 40.0) {
          /* disp(['Probably does not converge, current n is ' num2str(n)]); */
          flag = 1;
          exitg1 = 1;
        } else {
          c_i = P1_old_d->size[0];
          P1_old_d->size[0] = P1_d->size[0];
          emxEnsureCapacity_real_T(P1_old_d, c_i);
          loop_ub = P1_d->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P1_old_d->data[c_i] = P1_d->data[c_i];
          }

          c_i = P1_old_colidx->size[0];
          P1_old_colidx->size[0] = P1_colidx->size[0];
          emxEnsureCapacity_int32_T(P1_old_colidx, c_i);
          loop_ub = P1_colidx->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P1_old_colidx->data[c_i] = P1_colidx->data[c_i];
          }

          c_i = P1_old_rowidx->size[0];
          P1_old_rowidx->size[0] = P1_rowidx->size[0];
          emxEnsureCapacity_int32_T(P1_old_rowidx, c_i);
          loop_ub = P1_rowidx->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P1_old_rowidx->data[c_i] = P1_rowidx->data[c_i];
          }

          P1_old_m = P1_m;
          P1_old_n = P1_n;
          c_i = P2_old_d->size[0];
          P2_old_d->size[0] = P2_d->size[0];
          emxEnsureCapacity_real_T(P2_old_d, c_i);
          loop_ub = P2_d->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P2_old_d->data[c_i] = P2_d->data[c_i];
          }

          c_i = P2_old_colidx->size[0];
          P2_old_colidx->size[0] = P2_colidx->size[0];
          emxEnsureCapacity_int32_T(P2_old_colidx, c_i);
          loop_ub = P2_colidx->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P2_old_colidx->data[c_i] = P2_colidx->data[c_i];
          }

          c_i = P2_old_rowidx->size[0];
          P2_old_rowidx->size[0] = P2_rowidx->size[0];
          emxEnsureCapacity_int32_T(P2_old_rowidx, c_i);
          loop_ub = P2_rowidx->size[0];
          for (c_i = 0; c_i < loop_ub; c_i++) {
            P2_old_rowidx->data[c_i] = P2_rowidx->data[c_i];
          }

          P2_old_m = expl_temp.m;
          P2_old_n = expl_temp.n;

          /* fix n,every interation of k, u,P1,P2 will change value, so add 'old' */
        }
      }
    } while (exitg1 == 0);

    if (b_r == 1.0) {
      /* coupon rate payment */
      AccI = coupon_time_rate[x_tmp] * F;
      loop_ub = u->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        u->data[c_i] += AccI;
      }

      loop_ub = B->size[0];
      for (c_i = 0; c_i < loop_ub; c_i++) {
        B->data[c_i] += AccI;
      }
    }

    n++;
  }

  emxFree_real_T(&r2);
  c_emxFreeStruct_coder_internal_(&expl_temp);
  emxFree_real_T(&maxval);
  emxFree_boolean_T(&b_d);
  emxFree_real_T(&a);
  emxFree_real_T(&y);
  emxFree_int32_T(&b_c_rowidx);
  emxFree_int32_T(&b_c_colidx);
  emxFree_real_T(&b_c_d);
  emxFree_int32_T(&c_rowidx);
  emxFree_int32_T(&c_colidx);
  emxFree_real_T(&c_d);
  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_int32_T(&P2_rowidx);
  emxFree_int32_T(&P2_colidx);
  emxFree_real_T(&P2_d);
  emxFree_int32_T(&P1_rowidx);
  emxFree_int32_T(&P1_colidx);
  emxFree_real_T(&P1_d);
  emxFree_real_T(&u_old);
  emxFree_int32_T(&P2_old_rowidx);
  emxFree_int32_T(&P2_old_colidx);
  emxFree_real_T(&P2_old_d);
  emxFree_int32_T(&P1_old_rowidx);
  emxFree_int32_T(&P1_old_colidx);
  emxFree_real_T(&P1_old_d);
  emxFree_real_T(&u_n);
  emxFree_int32_T(&MB_rowidx);
  emxFree_int32_T(&MB_colidx);
  emxFree_real_T(&MB_d);
  emxFree_int32_T(&Mu_rowidx);
  emxFree_int32_T(&Mu_colidx);
  emxFree_real_T(&Mu_d);

  /*  */
  AccI = Ns * S0 / S[2400] + 1.0;
  if (AccI > 1.0) {
    U = u->data[static_cast<int>(std::ceil(AccI)) - 1] * (AccI - std::floor(AccI
      + 0.001)) + u->data[static_cast<int>(std::floor(AccI)) - 1] * (std::ceil
      (AccI + 0.001) - AccI);
  } else {
    U = u->data[0];
  }

  /*  disp(U) */
  return U;
}

/* End of code generation (funcconv.cpp) */
