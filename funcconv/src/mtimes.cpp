/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.cpp
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T
                     *a_colidx, const emxArray_int32_T *a_rowidx, int a_m, int
                     a_n, const double b[2401], emxArray_real_T *c)
{
  int i;
  int acol;
  double bc;
  int apend;
  int nap;
  int ap;
  int c_tmp;
  i = c->size[0];
  c->size[0] = a_m;
  emxEnsureCapacity_real_T(c, i);
  for (i = 0; i < a_m; i++) {
    c->data[i] = 0.0;
  }

  if ((a_n != 0) && (a_m != 0) && (a_colidx->data[a_colidx->size[0] - 1] - 1 !=
       0)) {
    for (acol = 0; acol < a_n; acol++) {
      bc = b[acol];
      i = a_colidx->data[acol];
      apend = a_colidx->data[acol + 1];
      nap = apend - a_colidx->data[acol];
      if (nap >= 4) {
        apend = (apend - nap) + ((nap / 4) << 2);
        for (ap = i; ap <= apend - 1; ap += 4) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
          c->data[a_rowidx->data[ap] - 1] += a_d->data[ap] * bc;
          c_tmp = a_rowidx->data[ap + 1] - 1;
          c->data[c_tmp] += a_d->data[ap + 1] * bc;
          c_tmp = a_rowidx->data[ap + 2] - 1;
          c->data[c_tmp] += a_d->data[ap + 2] * bc;
        }

        nap = a_colidx->data[acol + 1] - 1;
        for (ap = apend; ap <= nap; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
        }
      } else {
        apend--;
        for (ap = i; ap <= apend; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
        }
      }
    }
  }
}

void sparse_mtimes(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, int a_m, int a_n, const
                   emxArray_real_T *b, emxArray_real_T *c)
{
  int i;
  int acol;
  double bc;
  int apend;
  int nap;
  int ap;
  int c_tmp;
  i = c->size[0];
  c->size[0] = a_m;
  emxEnsureCapacity_real_T(c, i);
  for (i = 0; i < a_m; i++) {
    c->data[i] = 0.0;
  }

  if ((a_n != 0) && (a_m != 0) && (a_colidx->data[a_colidx->size[0] - 1] - 1 !=
       0)) {
    for (acol = 0; acol < a_n; acol++) {
      bc = b->data[acol];
      i = a_colidx->data[acol];
      apend = a_colidx->data[acol + 1];
      nap = apend - a_colidx->data[acol];
      if (nap >= 4) {
        apend = (apend - nap) + ((nap / 4) << 2);
        for (ap = i; ap <= apend - 1; ap += 4) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
          c->data[a_rowidx->data[ap] - 1] += a_d->data[ap] * bc;
          c_tmp = a_rowidx->data[ap + 1] - 1;
          c->data[c_tmp] += a_d->data[ap + 1] * bc;
          c_tmp = a_rowidx->data[ap + 2] - 1;
          c->data[c_tmp] += a_d->data[ap + 2] * bc;
        }

        nap = a_colidx->data[acol + 1] - 1;
        for (ap = apend; ap <= nap; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
        }
      } else {
        apend--;
        for (ap = i; ap <= apend; ap++) {
          c_tmp = a_rowidx->data[ap - 1] - 1;
          c->data[c_tmp] += a_d->data[ap - 1] * bc;
        }
      }
    }
  }
}

/* End of code generation (mtimes.cpp) */
