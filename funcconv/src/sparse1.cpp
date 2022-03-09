/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse1.cpp
 *
 * Code generation for function 'sparse1'
 *
 */

/* Include files */
#include "sparse1.h"
#include "CXSparseAPI.h"
#include "binOp.h"
#include "fillIn.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "makeCXSparseMatrix.h"
#include "rt_nonfinite.h"
#include "solve_from_lu.h"
#include "sparse.h"

/* Function Definitions */
void b_sparse_minus(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                    const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d,
                    const emxArray_int32_T *b_colidx, const emxArray_int32_T
                    *b_rowidx, int b_m, int b_n, emxArray_real_T *s_d,
                    emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx, int *
                    s_m, int *s_n)
{
  coder_internal_sparse expl_temp;
  int c;
  int didx;
  int b_c;
  int aidx_tmp;
  int aidx;
  int bidx_tmp;
  int bidx;
  bool moreAToDo;
  bool moreBToDo;
  double val;
  c_emxInitStruct_coder_internal_(&expl_temp);
  allocEqsizeBinop(a_colidx, b_colidx, b_n, b_m, &expl_temp);
  c = s_d->size[0];
  s_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(s_d, c);
  didx = expl_temp.d->size[0];
  for (c = 0; c < didx; c++) {
    s_d->data[c] = expl_temp.d->data[c];
  }

  c = s_colidx->size[0];
  s_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(s_colidx, c);
  didx = expl_temp.colidx->size[0];
  for (c = 0; c < didx; c++) {
    s_colidx->data[c] = expl_temp.colidx->data[c];
  }

  c = s_rowidx->size[0];
  s_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(s_rowidx, c);
  didx = expl_temp.rowidx->size[0];
  for (c = 0; c < didx; c++) {
    s_rowidx->data[c] = expl_temp.rowidx->data[c];
  }

  *s_m = expl_temp.m;
  *s_n = expl_temp.n;
  didx = 1;
  s_colidx->data[0] = 1;
  c_emxFreeStruct_coder_internal_(&expl_temp);
  for (c = 0; c < *s_n; c++) {
    b_c = c + 1;
    aidx_tmp = a_colidx->data[b_c - 1];
    aidx = aidx_tmp - 1;
    bidx_tmp = b_colidx->data[b_c - 1];
    bidx = bidx_tmp - 1;
    moreAToDo = (aidx_tmp < a_colidx->data[b_c]);
    moreBToDo = (bidx_tmp < b_colidx->data[b_c]);
    while (moreAToDo || moreBToDo) {
      while ((aidx + 1 < a_colidx->data[b_c]) && ((!moreBToDo) ||
              (a_rowidx->data[aidx] < b_rowidx->data[bidx]))) {
        if (a_d->data[aidx] != 0.0) {
          s_d->data[didx - 1] = a_d->data[aidx];
          s_rowidx->data[didx - 1] = a_rowidx->data[aidx];
          didx++;
        }

        aidx++;
      }

      moreAToDo = (aidx + 1 < a_colidx->data[b_c]);
      while ((bidx + 1 < b_colidx->data[b_c]) && ((!moreAToDo) ||
              (b_rowidx->data[bidx] < a_rowidx->data[aidx]))) {
        if (0.0 - b_d->data[bidx] != 0.0) {
          s_d->data[didx - 1] = 0.0 - b_d->data[bidx];
          s_rowidx->data[didx - 1] = b_rowidx->data[bidx];
          didx++;
        }

        bidx++;
      }

      while ((aidx + 1 < a_colidx->data[b_c]) && (bidx + 1 < b_colidx->data[b_c])
             && (a_rowidx->data[aidx] == b_rowidx->data[bidx])) {
        val = a_d->data[aidx] - b_d->data[bidx];
        if (val != 0.0) {
          s_d->data[didx - 1] = val;
          s_rowidx->data[didx - 1] = b_rowidx->data[bidx];
          didx++;
        }

        bidx++;
        aidx++;
      }

      moreAToDo = (aidx + 1 < a_colidx->data[b_c]);
      moreBToDo = (bidx + 1 < b_colidx->data[b_c]);
    }

    s_colidx->data[b_c] = didx;
  }
}

void b_sparse_plus(const emxArray_real_T *a_d, const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, const emxArray_real_T *b_d,
                   const emxArray_int32_T *b_colidx, const emxArray_int32_T
                   *b_rowidx, int b_m, int b_n, emxArray_real_T *s_d,
                   emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx, int
                   *s_m, int *s_n)
{
  coder_internal_sparse expl_temp;
  int c;
  int didx;
  int b_c;
  int aidx_tmp;
  int aidx;
  int bidx_tmp;
  int bidx;
  bool moreAToDo;
  bool moreBToDo;
  double val;
  c_emxInitStruct_coder_internal_(&expl_temp);
  allocEqsizeBinop(a_colidx, b_colidx, b_n, b_m, &expl_temp);
  c = s_d->size[0];
  s_d->size[0] = expl_temp.d->size[0];
  emxEnsureCapacity_real_T(s_d, c);
  didx = expl_temp.d->size[0];
  for (c = 0; c < didx; c++) {
    s_d->data[c] = expl_temp.d->data[c];
  }

  c = s_colidx->size[0];
  s_colidx->size[0] = expl_temp.colidx->size[0];
  emxEnsureCapacity_int32_T(s_colidx, c);
  didx = expl_temp.colidx->size[0];
  for (c = 0; c < didx; c++) {
    s_colidx->data[c] = expl_temp.colidx->data[c];
  }

  c = s_rowidx->size[0];
  s_rowidx->size[0] = expl_temp.rowidx->size[0];
  emxEnsureCapacity_int32_T(s_rowidx, c);
  didx = expl_temp.rowidx->size[0];
  for (c = 0; c < didx; c++) {
    s_rowidx->data[c] = expl_temp.rowidx->data[c];
  }

  *s_m = expl_temp.m;
  *s_n = expl_temp.n;
  didx = 1;
  s_colidx->data[0] = 1;
  c_emxFreeStruct_coder_internal_(&expl_temp);
  for (c = 0; c < *s_n; c++) {
    b_c = c + 1;
    aidx_tmp = a_colidx->data[b_c - 1];
    aidx = aidx_tmp - 1;
    bidx_tmp = b_colidx->data[b_c - 1];
    bidx = bidx_tmp - 1;
    moreAToDo = (aidx_tmp < a_colidx->data[b_c]);
    moreBToDo = (bidx_tmp < b_colidx->data[b_c]);
    while (moreAToDo || moreBToDo) {
      while ((aidx + 1 < a_colidx->data[b_c]) && ((!moreBToDo) ||
              (a_rowidx->data[aidx] < b_rowidx->data[bidx]))) {
        if (a_d->data[aidx] != 0.0) {
          s_d->data[didx - 1] = a_d->data[aidx];
          s_rowidx->data[didx - 1] = a_rowidx->data[aidx];
          didx++;
        }

        aidx++;
      }

      moreAToDo = (aidx + 1 < a_colidx->data[b_c]);
      while ((bidx + 1 < b_colidx->data[b_c]) && ((!moreAToDo) ||
              (b_rowidx->data[bidx] < a_rowidx->data[aidx]))) {
        if (b_d->data[bidx] != 0.0) {
          s_d->data[didx - 1] = b_d->data[bidx];
          s_rowidx->data[didx - 1] = b_rowidx->data[bidx];
          didx++;
        }

        bidx++;
      }

      while ((aidx + 1 < a_colidx->data[b_c]) && (bidx + 1 < b_colidx->data[b_c])
             && (a_rowidx->data[aidx] == b_rowidx->data[bidx])) {
        val = a_d->data[aidx] + b_d->data[bidx];
        if (val != 0.0) {
          s_d->data[didx - 1] = val;
          s_rowidx->data[didx - 1] = b_rowidx->data[bidx];
          didx++;
        }

        bidx++;
        aidx++;
      }

      moreAToDo = (aidx + 1 < a_colidx->data[b_c]);
      moreBToDo = (bidx + 1 < b_colidx->data[b_c]);
    }

    s_colidx->data[b_c] = didx;
  }
}

void b_sparse_times(double a, const emxArray_boolean_T *b_d, const
                    emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx,
                    int b_m, int b_n, coder_internal_sparse *s)
{
  emxArray_real_T *S;
  int nzs_tmp;
  int i;
  int numalloc;
  emxArray_real_T *tmpd;
  int idx;
  if (a * 0.0 == 0.0) {
    nzs_tmp = b_colidx->data[b_colidx->size[0] - 1];
    if (1 > b_colidx->data[b_colidx->size[0] - 1] - 1) {
      numalloc = 0;
    } else {
      numalloc = b_colidx->data[b_colidx->size[0] - 1] - 1;
    }

    emxInit_real_T(&tmpd, 1);
    i = tmpd->size[0];
    tmpd->size[0] = numalloc;
    emxEnsureCapacity_real_T(tmpd, i);
    for (i = 0; i < numalloc; i++) {
      tmpd->data[i] = a * static_cast<double>(b_d->data[i]);
    }

    s->m = b_m;
    s->n = b_n;
    if (b_colidx->data[b_colidx->size[0] - 1] - 1 >= 1) {
      numalloc = b_colidx->data[b_colidx->size[0] - 1] - 2;
    } else {
      numalloc = 0;
    }

    i = s->d->size[0];
    s->d->size[0] = numalloc + 1;
    emxEnsureCapacity_real_T(s->d, i);
    for (i = 0; i <= numalloc; i++) {
      s->d->data[i] = 0.0;
    }

    i = s->colidx->size[0];
    s->colidx->size[0] = b_n + 1;
    emxEnsureCapacity_int32_T(s->colidx, i);
    s->colidx->data[0] = 1;
    i = s->rowidx->size[0];
    s->rowidx->size[0] = numalloc + 1;
    emxEnsureCapacity_int32_T(s->rowidx, i);
    for (i = 0; i <= numalloc; i++) {
      s->rowidx->data[i] = 0;
    }

    for (numalloc = 0; numalloc < b_n; numalloc++) {
      s->colidx->data[numalloc + 1] = 1;
    }

    sparse_fillIn(s);
    if (1 > b_colidx->data[b_colidx->size[0] - 1] - 1) {
      numalloc = 1;
    } else {
      numalloc = b_colidx->data[b_colidx->size[0] - 1];
    }

    for (i = 0; i <= numalloc - 2; i++) {
      s->rowidx->data[i] = b_rowidx->data[i];
    }

    i = s->colidx->size[0];
    s->colidx->size[0] = b_colidx->size[0];
    emxEnsureCapacity_int32_T(s->colidx, i);
    numalloc = b_colidx->size[0];
    for (i = 0; i < numalloc; i++) {
      s->colidx->data[i] = b_colidx->data[i];
    }

    for (numalloc = 0; numalloc <= nzs_tmp - 2; numalloc++) {
      s->d->data[numalloc] = tmpd->data[numalloc];
    }

    emxFree_real_T(&tmpd);
    b_sparse_fillIn(s);
  } else {
    emxInit_real_T(&S, 2);
    i = S->size[0] * S->size[1];
    S->size[0] = b_m;
    S->size[1] = b_n;
    emxEnsureCapacity_real_T(S, i);
    numalloc = b_m * b_n;
    for (i = 0; i < numalloc; i++) {
      S->data[i] = rtNaN;
    }

    for (numalloc = 0; numalloc < b_n; numalloc++) {
      i = b_colidx->data[numalloc];
      nzs_tmp = b_colidx->data[numalloc + 1] - 1;
      for (idx = i; idx <= nzs_tmp; idx++) {
        S->data[(b_rowidx->data[idx - 1] + S->size[0] * numalloc) - 1] = a *
          static_cast<double>(b_d->data[idx - 1]);
      }
    }

    sparse(S, s->d, s->colidx, s->rowidx, &s->m, &s->n);
    emxFree_real_T(&S);
  }
}

void sparse_minus(const emxArray_real_T *a, const emxArray_real_T *b_d, const
                  emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx,
                  int b_m, int b_n, emxArray_real_T *s)
{
  int col;
  int idx;
  int row;
  col = s->size[0] * s->size[1];
  s->size[0] = b_m;
  s->size[1] = b_n;
  emxEnsureCapacity_real_T(s, col);
  idx = b_m * b_n;
  for (col = 0; col < idx; col++) {
    s->data[col] = 0.0;
  }

  for (col = 0; col < b_n; col++) {
    idx = b_colidx->data[col];
    for (row = 0; row < b_m; row++) {
      if ((idx < b_colidx->data[col + 1]) && (row + 1 == b_rowidx->data[idx - 1]))
      {
        s->data[row + s->size[0] * col] = a->data[row + a->size[0] * col] -
          b_d->data[idx - 1];
        idx++;
      } else {
        s->data[row + s->size[0] * col] = a->data[row + a->size[0] * col];
      }
    }
  }
}

void sparse_mldivide(const emxArray_real_T *A_d, const emxArray_int32_T
                     *A_colidx, const emxArray_int32_T *A_rowidx, int A_m, int
                     A_n, const double b_data[], const int b_size[1],
                     emxArray_real_T *y)
{
  int idx;
  coder_internal_sparse in;
  cs_di* cxA;
  cs_dis * S;
  cs_din * N;
  int numalloc;
  emxArray_int32_T *counts;
  int outridx_tmp;
  int outridx;
  if ((A_m == 0) || (A_n == 0)) {
    idx = y->size[0];
    y->size[0] = A_n;
    emxEnsureCapacity_real_T(y, idx);
    for (idx = 0; idx < A_n; idx++) {
      y->data[idx] = 0.0;
    }
  } else if (b_size[0] == A_n) {
    if (A_m < A_n) {
      c_emxInitStruct_coder_internal_(&in);
      in.m = A_n;
      in.n = A_m;
      if (A_colidx->data[A_colidx->size[0] - 1] - 1 >= 1) {
        numalloc = A_colidx->data[A_colidx->size[0] - 1] - 2;
      } else {
        numalloc = 0;
      }

      idx = in.d->size[0];
      in.d->size[0] = numalloc + 1;
      emxEnsureCapacity_real_T(in.d, idx);
      for (idx = 0; idx <= numalloc; idx++) {
        in.d->data[idx] = 0.0;
      }

      idx = in.colidx->size[0];
      in.colidx->size[0] = A_m + 1;
      emxEnsureCapacity_int32_T(in.colidx, idx);
      in.colidx->data[0] = 1;
      idx = in.rowidx->size[0];
      in.rowidx->size[0] = numalloc + 1;
      emxEnsureCapacity_int32_T(in.rowidx, idx);
      for (idx = 0; idx <= numalloc; idx++) {
        in.rowidx->data[idx] = 0;
      }

      for (numalloc = 0; numalloc < A_m; numalloc++) {
        in.colidx->data[numalloc + 1] = 1;
      }

      sparse_fillIn(&in);
      if ((A_m != 0) && (A_n != 0)) {
        numalloc = in.colidx->size[0];
        for (idx = 0; idx < numalloc; idx++) {
          in.colidx->data[idx] = 0;
        }

        idx = A_colidx->data[A_colidx->size[0] - 1];
        for (numalloc = 0; numalloc <= idx - 2; numalloc++) {
          in.colidx->data[A_rowidx->data[numalloc]]++;
        }

        in.colidx->data[0] = 1;
        idx = A_m + 1;
        for (numalloc = 2; numalloc <= idx; numalloc++) {
          in.colidx->data[numalloc - 1] += in.colidx->data[numalloc - 2];
        }

        emxInit_int32_T(&counts, 1);
        idx = counts->size[0];
        counts->size[0] = A_m;
        emxEnsureCapacity_int32_T(counts, idx);
        for (idx = 0; idx < A_m; idx++) {
          counts->data[idx] = 0;
        }

        for (numalloc = 0; numalloc < A_n; numalloc++) {
          for (idx = A_colidx->data[numalloc] - 1; idx + 1 < A_colidx->
               data[numalloc + 1]; idx++) {
            outridx_tmp = counts->data[A_rowidx->data[idx] - 1];
            outridx = (outridx_tmp + in.colidx->data[A_rowidx->data[idx] - 1]) -
              1;
            in.d->data[outridx] = A_d->data[idx];
            in.rowidx->data[outridx] = numalloc + 1;
            counts->data[A_rowidx->data[idx] - 1] = outridx_tmp + 1;
          }
        }

        emxFree_int32_T(&counts);
      }

      cxA = makeCXSparseMatrix(in.colidx->data[in.colidx->size[0] - 1] - 1, in.n,
        in.m, &in.colidx->data[0], &in.rowidx->data[0], &in.d->data[0]);
      c_emxFreeStruct_coder_internal_(&in);
    } else {
      cxA = makeCXSparseMatrix(A_colidx->data[A_colidx->size[0] - 1] - 1, A_n,
        A_m, &A_colidx->data[0], &A_rowidx->data[0], &A_d->data[0]);
    }

    S = cs_di_sqr(2, cxA, 0);
    N = cs_di_lu(cxA, S, 1);
    cs_di_spfree(cxA);
    if (N == NULL) {
      cs_di_sfree(S);
      cs_di_nfree(N);
      CXSparseAPI_iteratedQR(A_d, A_colidx, A_rowidx, A_m, A_n, b_data, b_size,
        A_n, y);
    } else {
      idx = y->size[0];
      y->size[0] = b_size[0];
      emxEnsureCapacity_real_T(y, idx);
      numalloc = b_size[0];
      for (idx = 0; idx < numalloc; idx++) {
        y->data[idx] = b_data[idx];
      }

      solve_from_lu_di(N, S, (double *)&y->data[0], b_size[0]);
      cs_di_sfree(S);
      cs_di_nfree(N);
    }
  } else {
    CXSparseAPI_iteratedQR(A_d, A_colidx, A_rowidx, A_m, A_n, b_data, b_size,
      A_n, y);
  }
}

void sparse_plus(const emxArray_real_T *a, const emxArray_real_T *b_d, const
                 emxArray_int32_T *b_colidx, const emxArray_int32_T *b_rowidx,
                 int b_m, int b_n, emxArray_real_T *s)
{
  int col;
  int idx;
  int row;
  col = s->size[0] * s->size[1];
  s->size[0] = b_m;
  s->size[1] = b_n;
  emxEnsureCapacity_real_T(s, col);
  idx = b_m * b_n;
  for (col = 0; col < idx; col++) {
    s->data[col] = 0.0;
  }

  for (col = 0; col < b_n; col++) {
    idx = b_colidx->data[col];
    for (row = 0; row < b_m; row++) {
      if ((idx < b_colidx->data[col + 1]) && (row + 1 == b_rowidx->data[idx - 1]))
      {
        s->data[row + s->size[0] * col] = a->data[row + a->size[0] * col] +
          b_d->data[idx - 1];
        idx++;
      } else {
        s->data[row + s->size[0] * col] = a->data[row + a->size[0] * col];
      }
    }
  }
}

void sparse_times(double a, const emxArray_real_T *b_d, const emxArray_int32_T
                  *b_colidx, const emxArray_int32_T *b_rowidx, int b_m, int b_n,
                  coder_internal_sparse *s)
{
  emxArray_real_T *S;
  int nzs_tmp;
  int i;
  int numalloc;
  emxArray_real_T *tmpd;
  int idx;
  if (a * 0.0 == 0.0) {
    nzs_tmp = b_colidx->data[b_colidx->size[0] - 1];
    if (1 > b_colidx->data[b_colidx->size[0] - 1] - 1) {
      numalloc = 0;
    } else {
      numalloc = b_colidx->data[b_colidx->size[0] - 1] - 1;
    }

    emxInit_real_T(&tmpd, 1);
    i = tmpd->size[0];
    tmpd->size[0] = numalloc;
    emxEnsureCapacity_real_T(tmpd, i);
    for (i = 0; i < numalloc; i++) {
      tmpd->data[i] = a * b_d->data[i];
    }

    s->m = b_m;
    s->n = b_n;
    if (b_colidx->data[b_colidx->size[0] - 1] - 1 >= 1) {
      numalloc = b_colidx->data[b_colidx->size[0] - 1] - 2;
    } else {
      numalloc = 0;
    }

    i = s->d->size[0];
    s->d->size[0] = numalloc + 1;
    emxEnsureCapacity_real_T(s->d, i);
    for (i = 0; i <= numalloc; i++) {
      s->d->data[i] = 0.0;
    }

    i = s->colidx->size[0];
    s->colidx->size[0] = b_n + 1;
    emxEnsureCapacity_int32_T(s->colidx, i);
    s->colidx->data[0] = 1;
    i = s->rowidx->size[0];
    s->rowidx->size[0] = numalloc + 1;
    emxEnsureCapacity_int32_T(s->rowidx, i);
    for (i = 0; i <= numalloc; i++) {
      s->rowidx->data[i] = 0;
    }

    for (numalloc = 0; numalloc < b_n; numalloc++) {
      s->colidx->data[numalloc + 1] = 1;
    }

    sparse_fillIn(s);
    if (1 > b_colidx->data[b_colidx->size[0] - 1] - 1) {
      numalloc = 1;
    } else {
      numalloc = b_colidx->data[b_colidx->size[0] - 1];
    }

    for (i = 0; i <= numalloc - 2; i++) {
      s->rowidx->data[i] = b_rowidx->data[i];
    }

    i = s->colidx->size[0];
    s->colidx->size[0] = b_colidx->size[0];
    emxEnsureCapacity_int32_T(s->colidx, i);
    numalloc = b_colidx->size[0];
    for (i = 0; i < numalloc; i++) {
      s->colidx->data[i] = b_colidx->data[i];
    }

    for (numalloc = 0; numalloc <= nzs_tmp - 2; numalloc++) {
      s->d->data[numalloc] = tmpd->data[numalloc];
    }

    emxFree_real_T(&tmpd);
    b_sparse_fillIn(s);
  } else {
    emxInit_real_T(&S, 2);
    i = S->size[0] * S->size[1];
    S->size[0] = b_m;
    S->size[1] = b_n;
    emxEnsureCapacity_real_T(S, i);
    numalloc = b_m * b_n;
    for (i = 0; i < numalloc; i++) {
      S->data[i] = rtNaN;
    }

    for (numalloc = 0; numalloc < b_n; numalloc++) {
      i = b_colidx->data[numalloc];
      nzs_tmp = b_colidx->data[numalloc + 1] - 1;
      for (idx = i; idx <= nzs_tmp; idx++) {
        S->data[(b_rowidx->data[idx - 1] + S->size[0] * numalloc) - 1] = a *
          b_d->data[idx - 1];
      }
    }

    sparse(S, s->d, s->colidx, s->rowidx, &s->m, &s->n);
    emxFree_real_T(&S);
  }
}

/* End of code generation (sparse1.cpp) */
