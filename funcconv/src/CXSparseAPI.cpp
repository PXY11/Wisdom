/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CXSparseAPI.cpp
 *
 * Code generation for function 'CXSparseAPI'
 *
 */

/* Include files */
#include "CXSparseAPI.h"
#include "fillIn.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "makeCXSparseMatrix.h"
#include "rt_nonfinite.h"
#include "solve_from_qr.h"

/* Function Definitions */
void CXSparseAPI_iteratedQR(const emxArray_real_T *A_d, const emxArray_int32_T
  *A_colidx, const emxArray_int32_T *A_rowidx, int A_m, int A_n, const double
  b_data[], const int b_size[1], int n, emxArray_real_T *out)
{
  coder_internal_sparse in;
  cs_di* cxA;
  cs_dis * S;
  cs_din * N;
  int numalloc;
  double tol;
  int idx;
  emxArray_real_T *outBuff;
  emxArray_int32_T *counts;
  int outridx_tmp;
  int outridx;
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
          outridx = (outridx_tmp + in.colidx->data[A_rowidx->data[idx] - 1]) - 1;
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
    cxA = makeCXSparseMatrix(A_colidx->data[A_colidx->size[0] - 1] - 1, A_n, A_m,
      &A_colidx->data[0], &A_rowidx->data[0], &A_d->data[0]);
  }

  S = cs_di_sqr(2, cxA, 1);
  N = cs_di_qr(cxA, S);
  cs_di_spfree(cxA);
  qr_rank_di(N, &tol);
  idx = out->size[0];
  out->size[0] = n;
  emxEnsureCapacity_real_T(out, idx);
  emxInit_real_T(&outBuff, 1);
  if (b_size[0] < n) {
    idx = outBuff->size[0];
    outBuff->size[0] = n;
    emxEnsureCapacity_real_T(outBuff, idx);
  } else {
    idx = outBuff->size[0];
    outBuff->size[0] = b_size[0];
    emxEnsureCapacity_real_T(outBuff, idx);
  }

  numalloc = b_size[0];
  for (idx = 0; idx < numalloc; idx++) {
    outBuff->data[idx] = b_data[idx];
  }

  solve_from_qr_di(N, S, (double *)&outBuff->data[0], b_size[0], n);
  if (1 > n) {
    numalloc = 0;
  } else {
    numalloc = n;
  }

  for (idx = 0; idx < numalloc; idx++) {
    out->data[idx] = outBuff->data[idx];
  }

  emxFree_real_T(&outBuff);
  cs_di_sfree(S);
  cs_di_nfree(N);
}

/* End of code generation (CXSparseAPI.cpp) */
