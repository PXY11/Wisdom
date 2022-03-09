/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spdiags.cpp
 *
 * Code generation for function 'spdiags'
 *
 */

/* Include files */
#include "spdiags.h"
#include "fillIn.h"
#include "funcconv.h"
#include "funcconv_emxutil.h"
#include "introsort.h"
#include "rt_nonfinite.h"
#include <cmath>

/* Function Definitions */
void b_spdiags(const bool arg1[2401], double arg3, double arg4,
               emxArray_boolean_T *res1_d, emxArray_int32_T *res1_colidx,
               emxArray_int32_T *res1_rowidx, int *res1_m, int *res1_n)
{
  int ns;
  double len[2];
  emxArray_int32_T *aRows;
  double minAdjustedDim_data_idx_0;
  int ridx;
  int i;
  emxArray_int32_T *aCols;
  emxArray_boolean_T *aDat;
  emxArray_real_T *idx;
  emxArray_real_T *b_i;
  emxArray_int32_T *r;
  int currRowIdx;
  emxArray_int32_T *ridxInt;
  int nc;
  emxArray_int32_T *sortedIndices;
  cell_wrap_2 this_tunableEnvironment[2];
  emxArray_int32_T *t;
  int res1_n_tmp;
  int c;
  bool val;
  ns = 0;
  if ((0.0 >= -arg3 + 1.0) && (0.0 <= arg4 - 1.0)) {
    ns = 1;
  }

  len[0] = 0.0;
  if (0 <= ns - 1) {
    if ((arg3 < arg4) || rtIsNaN(arg4)) {
      minAdjustedDim_data_idx_0 = arg3;
    } else {
      minAdjustedDim_data_idx_0 = arg4;
    }

    len[1] = (minAdjustedDim_data_idx_0 - 1.0) + 1.0;
  }

  emxInit_int32_T(&aRows, 1);
  ridx = static_cast<int>(len[ns]);
  i = aRows->size[0];
  aRows->size[0] = ridx;
  emxEnsureCapacity_int32_T(aRows, i);
  for (i = 0; i < ridx; i++) {
    aRows->data[i] = 0;
  }

  emxInit_int32_T(&aCols, 1);
  i = aCols->size[0];
  aCols->size[0] = ridx;
  emxEnsureCapacity_int32_T(aCols, i);
  for (i = 0; i < ridx; i++) {
    aCols->data[i] = 0;
  }

  emxInit_boolean_T(&aDat, 1);
  i = aDat->size[0];
  aDat->size[0] = ridx;
  emxEnsureCapacity_boolean_T(aDat, i);
  for (i = 0; i < ridx; i++) {
    aDat->data[i] = false;
  }

  emxInit_real_T(&idx, 2);
  emxInit_real_T(&b_i, 1);
  emxInit_int32_T(&r, 2);
  for (currRowIdx = 0; currRowIdx < ns; currRowIdx++) {
    if (rtIsNaN(minAdjustedDim_data_idx_0)) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else if (minAdjustedDim_data_idx_0 < 1.0) {
      idx->size[0] = 1;
      idx->size[1] = 0;
    } else if (rtIsInf(minAdjustedDim_data_idx_0) && (1.0 ==
                minAdjustedDim_data_idx_0)) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      ridx = static_cast<int>(std::floor(minAdjustedDim_data_idx_0 - 1.0));
      idx->size[1] = ridx + 1;
      emxEnsureCapacity_real_T(idx, i);
      for (i = 0; i <= ridx; i++) {
        idx->data[i] = static_cast<double>(i) + 1.0;
      }
    }

    ridx = idx->size[1];
    i = b_i->size[0];
    b_i->size[0] = idx->size[1];
    emxEnsureCapacity_real_T(b_i, i);
    for (i = 0; i < ridx; i++) {
      b_i->data[i] = idx->data[i];
    }

    if (rtIsNaN(len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else if (len[1] < 1.0) {
      idx->size[0] = 1;
      idx->size[1] = 0;
    } else if (rtIsInf(len[1]) && (1.0 == len[1])) {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, i);
      idx->data[0] = rtNaN;
    } else {
      i = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      ridx = static_cast<int>(std::floor(len[1] - 1.0));
      idx->size[1] = ridx + 1;
      emxEnsureCapacity_real_T(idx, i);
      for (i = 0; i <= ridx; i++) {
        idx->data[i] = static_cast<double>(i) + 1.0;
      }
    }

    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, i);
    ridx = idx->size[0] * idx->size[1];
    for (i = 0; i < ridx; i++) {
      r->data[i] = static_cast<int>(idx->data[i]);
    }

    ridx = r->size[0] * r->size[1];
    for (i = 0; i < ridx; i++) {
      aRows->data[r->data[i] - 1] = static_cast<int>(b_i->data[i]);
    }

    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, i);
    ridx = idx->size[0] * idx->size[1];
    for (i = 0; i < ridx; i++) {
      r->data[i] = static_cast<int>(idx->data[i]);
    }

    ridx = r->size[0] * r->size[1];
    for (i = 0; i < ridx; i++) {
      aCols->data[r->data[i] - 1] = static_cast<int>(b_i->data[i]);
    }

    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, i);
    ridx = idx->size[0] * idx->size[1];
    for (i = 0; i < ridx; i++) {
      r->data[i] = static_cast<int>(idx->data[i]);
    }

    ridx = r->size[1];
    for (i = 0; i < ridx; i++) {
      aDat->data[r->data[i] - 1] = arg1[static_cast<int>(b_i->data[i]) - 1];
    }
  }

  emxFree_int32_T(&r);
  emxFree_real_T(&b_i);
  emxFree_real_T(&idx);
  emxInit_int32_T(&ridxInt, 1);
  nc = aCols->size[0];
  ns = aRows->size[0];
  i = ridxInt->size[0];
  ridxInt->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(ridxInt, i);
  for (currRowIdx = 0; currRowIdx < ns; currRowIdx++) {
    ridxInt->data[currRowIdx] = aRows->data[currRowIdx];
  }

  ns = aCols->size[0];
  i = aRows->size[0];
  aRows->size[0] = aCols->size[0];
  emxEnsureCapacity_int32_T(aRows, i);
  for (currRowIdx = 0; currRowIdx < ns; currRowIdx++) {
    aRows->data[currRowIdx] = aCols->data[currRowIdx];
  }

  emxInit_int32_T(&sortedIndices, 1);
  i = sortedIndices->size[0];
  sortedIndices->size[0] = aCols->size[0];
  emxEnsureCapacity_int32_T(sortedIndices, i);
  for (currRowIdx = 0; currRowIdx < nc; currRowIdx++) {
    sortedIndices->data[currRowIdx] = currRowIdx + 1;
  }

  emxInitMatrix_cell_wrap_2(this_tunableEnvironment);
  i = this_tunableEnvironment[0].f1->size[0];
  this_tunableEnvironment[0].f1->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[0].f1, i);
  ridx = aRows->size[0];
  for (i = 0; i < ridx; i++) {
    this_tunableEnvironment[0].f1->data[i] = aRows->data[i];
  }

  i = this_tunableEnvironment[1].f1->size[0];
  this_tunableEnvironment[1].f1->size[0] = ridxInt->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[1].f1, i);
  ridx = ridxInt->size[0];
  for (i = 0; i < ridx; i++) {
    this_tunableEnvironment[1].f1->data[i] = ridxInt->data[i];
  }

  emxInit_int32_T(&t, 1);
  introsort(sortedIndices, aRows->size[0], this_tunableEnvironment);
  ns = aRows->size[0];
  i = t->size[0];
  t->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(t, i);
  ridx = aRows->size[0];
  emxFreeMatrix_cell_wrap_2(this_tunableEnvironment);
  for (i = 0; i < ridx; i++) {
    t->data[i] = aRows->data[i];
  }

  for (currRowIdx = 0; currRowIdx < ns; currRowIdx++) {
    aRows->data[currRowIdx] = t->data[sortedIndices->data[currRowIdx] - 1];
  }

  ns = ridxInt->size[0];
  i = t->size[0];
  t->size[0] = ridxInt->size[0];
  emxEnsureCapacity_int32_T(t, i);
  ridx = ridxInt->size[0];
  for (i = 0; i < ridx; i++) {
    t->data[i] = ridxInt->data[i];
  }

  for (currRowIdx = 0; currRowIdx < ns; currRowIdx++) {
    ridxInt->data[currRowIdx] = t->data[sortedIndices->data[currRowIdx] - 1];
  }

  emxFree_int32_T(&t);
  res1_n_tmp = static_cast<int>(arg4);
  if (aCols->size[0] >= 1) {
    ns = aCols->size[0];
  } else {
    ns = 1;
  }

  i = res1_d->size[0];
  res1_d->size[0] = ns;
  emxEnsureCapacity_boolean_T(res1_d, i);
  emxFree_int32_T(&aCols);
  for (i = 0; i < ns; i++) {
    res1_d->data[i] = false;
  }

  i = res1_colidx->size[0];
  res1_colidx->size[0] = res1_n_tmp + 1;
  emxEnsureCapacity_int32_T(res1_colidx, i);
  res1_colidx->data[0] = 1;
  i = res1_rowidx->size[0];
  res1_rowidx->size[0] = ns;
  emxEnsureCapacity_int32_T(res1_rowidx, i);
  for (i = 0; i < ns; i++) {
    res1_rowidx->data[i] = 0;
  }

  ns = 0;
  for (c = 0; c < res1_n_tmp; c++) {
    ridx = c + 1;
    while ((ns + 1 <= nc) && (aRows->data[ns] == ridx)) {
      res1_rowidx->data[ns] = ridxInt->data[ns];
      ns++;
    }

    res1_colidx->data[ridx] = ns + 1;
  }

  emxFree_int32_T(&ridxInt);
  emxFree_int32_T(&aRows);
  for (currRowIdx = 0; currRowIdx < nc; currRowIdx++) {
    res1_d->data[currRowIdx] = aDat->data[sortedIndices->data[currRowIdx] - 1];
  }

  emxFree_int32_T(&sortedIndices);
  emxFree_boolean_T(&aDat);
  ns = 1;
  i = res1_colidx->size[0];
  for (c = 0; c <= i - 2; c++) {
    ridx = res1_colidx->data[c];
    res1_colidx->data[c] = ns;
    while (ridx < res1_colidx->data[c + 1]) {
      currRowIdx = res1_rowidx->data[ridx - 1];
      val = res1_d->data[ridx - 1];
      ridx++;
      if (val) {
        res1_d->data[ns - 1] = true;
        res1_rowidx->data[ns - 1] = currRowIdx;
        ns++;
      }
    }
  }

  res1_colidx->data[res1_colidx->size[0] - 1] = ns;
  *res1_m = static_cast<int>(arg3);
  *res1_n = res1_n_tmp;
}

void spdiags(const double arg1_data[], const int arg1_size[2], double arg3,
             double arg4, coder_internal_sparse *res1)
{
  int trueCount;
  int nm1d2;
  signed char d_data[3];
  int mGEn;
  double len_data[4];
  emxArray_int32_T *aRows;
  int idx_tmp;
  int loop_ub;
  signed char i;
  signed char maxNegD_data[3];
  double u1;
  emxArray_int32_T *aCols;
  double d;
  double minAdjustedDim_data[3];
  emxArray_real_T *aDat;
  emxArray_real_T *idx;
  emxArray_real_T *b_i;
  emxArray_int32_T *r;
  int k;
  emxArray_int32_T *ridxInt;
  bool guard1 = false;
  emxArray_int32_T *sortedIndices;
  double ndbl;
  double apnd;
  double cdiff;
  double u0;
  cell_wrap_2 this_tunableEnvironment[2];
  int n;
  emxArray_int32_T *t;
  trueCount = 0;
  if ((-1.0 >= -arg3 + 1.0) && (-1.0 <= arg4 - 1.0)) {
    trueCount = 1;
  }

  if ((0.0 >= -arg3 + 1.0) && (0.0 <= arg4 - 1.0)) {
    trueCount++;
  }

  if ((1.0 >= -arg3 + 1.0) && (1.0 <= arg4 - 1.0)) {
    trueCount++;
  }

  nm1d2 = 0;
  if ((-1.0 >= -arg3 + 1.0) && (-1.0 <= arg4 - 1.0)) {
    d_data[0] = -1;
    nm1d2 = 1;
  }

  if ((0.0 >= -arg3 + 1.0) && (0.0 <= arg4 - 1.0)) {
    d_data[nm1d2] = 0;
    nm1d2++;
  }

  if ((1.0 >= -arg3 + 1.0) && (1.0 <= arg4 - 1.0)) {
    d_data[nm1d2] = 1;
  }

  mGEn = (arg3 >= arg4);
  len_data[0] = 0.0;
  for (nm1d2 = 0; nm1d2 < trueCount; nm1d2++) {
    idx_tmp = 1 - d_data[nm1d2];
    if (1 > idx_tmp) {
      i = 1;
      maxNegD_data[nm1d2] = 1;
    } else {
      i = static_cast<signed char>(idx_tmp);
      maxNegD_data[nm1d2] = static_cast<signed char>(idx_tmp);
    }

    u1 = arg4 - static_cast<double>(d_data[nm1d2]);
    if ((arg3 < u1) || rtIsNaN(u1)) {
      d = arg3;
    } else {
      d = u1;
    }

    minAdjustedDim_data[nm1d2] = d;
    len_data[nm1d2 + 1] = ((len_data[nm1d2] + d) - static_cast<double>(i)) + 1.0;
  }

  emxInit_int32_T(&aRows, 1);
  idx_tmp = aRows->size[0];
  loop_ub = static_cast<int>(len_data[trueCount]);
  aRows->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aRows, idx_tmp);
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    aRows->data[idx_tmp] = 0;
  }

  emxInit_int32_T(&aCols, 1);
  idx_tmp = aCols->size[0];
  aCols->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(aCols, idx_tmp);
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    aCols->data[idx_tmp] = 0;
  }

  emxInit_real_T(&aDat, 1);
  idx_tmp = aDat->size[0];
  aDat->size[0] = loop_ub;
  emxEnsureCapacity_real_T(aDat, idx_tmp);
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    aDat->data[idx_tmp] = 0.0;
  }

  emxInit_real_T(&idx, 2);
  emxInit_real_T(&b_i, 1);
  emxInit_int32_T(&r, 2);
  for (k = 0; k < trueCount; k++) {
    if (rtIsNaN(minAdjustedDim_data[k])) {
      idx_tmp = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, idx_tmp);
      idx->data[0] = rtNaN;
    } else if (minAdjustedDim_data[k] < maxNegD_data[k]) {
      idx->size[0] = 1;
      idx->size[1] = 0;
    } else if (rtIsInf(minAdjustedDim_data[k]) && (maxNegD_data[k] ==
                minAdjustedDim_data[k])) {
      idx_tmp = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, idx_tmp);
      idx->data[0] = rtNaN;
    } else {
      idx_tmp = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      loop_ub = static_cast<int>(std::floor(minAdjustedDim_data[k] -
        static_cast<double>(maxNegD_data[k])));
      idx->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(idx, idx_tmp);
      for (idx_tmp = 0; idx_tmp <= loop_ub; idx_tmp++) {
        idx->data[idx_tmp] = static_cast<double>(maxNegD_data[k]) + static_cast<
          double>(idx_tmp);
      }
    }

    idx_tmp = b_i->size[0];
    b_i->size[0] = idx->size[1];
    emxEnsureCapacity_real_T(b_i, idx_tmp);
    loop_ub = idx->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      b_i->data[idx_tmp] = idx->data[idx_tmp];
    }

    guard1 = false;
    if (rtIsNaN(len_data[k] + 1.0)) {
      guard1 = true;
    } else {
      d = len_data[k + 1];
      if (rtIsNaN(d)) {
        guard1 = true;
      } else if (d < len_data[k] + 1.0) {
        idx->size[0] = 1;
        idx->size[1] = 0;
      } else if ((rtIsInf(len_data[k] + 1.0) || rtIsInf(d)) && (len_data[k] +
                  1.0 == d)) {
        idx_tmp = idx->size[0] * idx->size[1];
        idx->size[0] = 1;
        idx->size[1] = 1;
        emxEnsureCapacity_real_T(idx, idx_tmp);
        idx->data[0] = rtNaN;
      } else if (std::floor(len_data[k] + 1.0) == len_data[k] + 1.0) {
        idx_tmp = idx->size[0] * idx->size[1];
        idx->size[0] = 1;
        loop_ub = static_cast<int>(std::floor(d - (len_data[k] + 1.0)));
        idx->size[1] = loop_ub + 1;
        emxEnsureCapacity_real_T(idx, idx_tmp);
        for (idx_tmp = 0; idx_tmp <= loop_ub; idx_tmp++) {
          idx->data[idx_tmp] = (len_data[k] + 1.0) + static_cast<double>(idx_tmp);
        }
      } else {
        ndbl = std::floor((d - (len_data[k] + 1.0)) + 0.5);
        apnd = (len_data[k] + 1.0) + ndbl;
        cdiff = apnd - d;
        u0 = std::abs(len_data[k] + 1.0);
        u1 = std::abs(d);
        if ((u0 > u1) || rtIsNaN(u1)) {
          u1 = u0;
        }

        if (std::abs(cdiff) < 4.4408920985006262E-16 * u1) {
          ndbl++;
          apnd = d;
        } else if (cdiff > 0.0) {
          apnd = (len_data[k] + 1.0) + (ndbl - 1.0);
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = static_cast<int>(ndbl);
        } else {
          n = 0;
        }

        idx_tmp = idx->size[0] * idx->size[1];
        idx->size[0] = 1;
        idx->size[1] = n;
        emxEnsureCapacity_real_T(idx, idx_tmp);
        if (n > 0) {
          idx->data[0] = len_data[k] + 1.0;
          if (n > 1) {
            idx->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (loop_ub = 0; loop_ub <= nm1d2 - 2; loop_ub++) {
              idx_tmp = loop_ub + 1;
              idx->data[loop_ub + 1] = (len_data[k] + 1.0) + static_cast<double>
                (idx_tmp);
              idx->data[(n - loop_ub) - 2] = apnd - static_cast<double>(idx_tmp);
            }

            if (nm1d2 << 1 == n - 1) {
              idx->data[nm1d2] = ((len_data[k] + 1.0) + apnd) / 2.0;
            } else {
              idx->data[nm1d2] = (len_data[k] + 1.0) + static_cast<double>(nm1d2);
              idx->data[nm1d2 + 1] = apnd - static_cast<double>(nm1d2);
            }
          }
        }
      }
    }

    if (guard1) {
      idx_tmp = idx->size[0] * idx->size[1];
      idx->size[0] = 1;
      idx->size[1] = 1;
      emxEnsureCapacity_real_T(idx, idx_tmp);
      idx->data[0] = rtNaN;
    }

    idx_tmp = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, idx_tmp);
    loop_ub = idx->size[0] * idx->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      r->data[idx_tmp] = static_cast<int>(idx->data[idx_tmp]);
    }

    loop_ub = r->size[0] * r->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      aRows->data[r->data[idx_tmp] - 1] = static_cast<int>(b_i->data[idx_tmp]);
    }

    idx_tmp = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, idx_tmp);
    loop_ub = idx->size[0] * idx->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      r->data[idx_tmp] = static_cast<int>(idx->data[idx_tmp]);
    }

    loop_ub = r->size[0] * r->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      aCols->data[r->data[idx_tmp] - 1] = static_cast<int>((b_i->data[idx_tmp] +
        static_cast<double>(d_data[k])));
    }

    idx_tmp = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = idx->size[1];
    emxEnsureCapacity_int32_T(r, idx_tmp);
    loop_ub = idx->size[0] * idx->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      r->data[idx_tmp] = static_cast<int>(idx->data[idx_tmp]);
    }

    nm1d2 = mGEn * d_data[k];
    loop_ub = r->size[1];
    for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
      aDat->data[r->data[idx_tmp] - 1] = arg1_data[(static_cast<int>((b_i->
        data[idx_tmp] + static_cast<double>(nm1d2))) + arg1_size[0] * k) - 1];
    }
  }

  emxFree_int32_T(&r);
  emxFree_real_T(&b_i);
  emxFree_real_T(&idx);
  emxInit_int32_T(&ridxInt, 1);
  mGEn = aCols->size[0];
  nm1d2 = aRows->size[0];
  idx_tmp = ridxInt->size[0];
  ridxInt->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(ridxInt, idx_tmp);
  for (k = 0; k < nm1d2; k++) {
    ridxInt->data[k] = aRows->data[k];
  }

  nm1d2 = aCols->size[0];
  idx_tmp = aRows->size[0];
  aRows->size[0] = aCols->size[0];
  emxEnsureCapacity_int32_T(aRows, idx_tmp);
  for (k = 0; k < nm1d2; k++) {
    aRows->data[k] = aCols->data[k];
  }

  emxInit_int32_T(&sortedIndices, 1);
  idx_tmp = sortedIndices->size[0];
  sortedIndices->size[0] = aCols->size[0];
  emxEnsureCapacity_int32_T(sortedIndices, idx_tmp);
  for (k = 0; k < mGEn; k++) {
    sortedIndices->data[k] = k + 1;
  }

  emxInitMatrix_cell_wrap_2(this_tunableEnvironment);
  idx_tmp = this_tunableEnvironment[0].f1->size[0];
  this_tunableEnvironment[0].f1->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[0].f1, idx_tmp);
  loop_ub = aRows->size[0];
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    this_tunableEnvironment[0].f1->data[idx_tmp] = aRows->data[idx_tmp];
  }

  idx_tmp = this_tunableEnvironment[1].f1->size[0];
  this_tunableEnvironment[1].f1->size[0] = ridxInt->size[0];
  emxEnsureCapacity_int32_T(this_tunableEnvironment[1].f1, idx_tmp);
  loop_ub = ridxInt->size[0];
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    this_tunableEnvironment[1].f1->data[idx_tmp] = ridxInt->data[idx_tmp];
  }

  emxInit_int32_T(&t, 1);
  introsort(sortedIndices, aRows->size[0], this_tunableEnvironment);
  nm1d2 = aRows->size[0];
  idx_tmp = t->size[0];
  t->size[0] = aRows->size[0];
  emxEnsureCapacity_int32_T(t, idx_tmp);
  loop_ub = aRows->size[0];
  emxFreeMatrix_cell_wrap_2(this_tunableEnvironment);
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    t->data[idx_tmp] = aRows->data[idx_tmp];
  }

  for (k = 0; k < nm1d2; k++) {
    aRows->data[k] = t->data[sortedIndices->data[k] - 1];
  }

  nm1d2 = ridxInt->size[0];
  idx_tmp = t->size[0];
  t->size[0] = ridxInt->size[0];
  emxEnsureCapacity_int32_T(t, idx_tmp);
  loop_ub = ridxInt->size[0];
  for (idx_tmp = 0; idx_tmp < loop_ub; idx_tmp++) {
    t->data[idx_tmp] = ridxInt->data[idx_tmp];
  }

  for (k = 0; k < nm1d2; k++) {
    ridxInt->data[k] = t->data[sortedIndices->data[k] - 1];
  }

  emxFree_int32_T(&t);
  res1->m = static_cast<int>(arg3);
  n = static_cast<int>(arg4);
  res1->n = n;
  if (aCols->size[0] >= 1) {
    nm1d2 = aCols->size[0];
  } else {
    nm1d2 = 1;
  }

  idx_tmp = res1->d->size[0];
  res1->d->size[0] = nm1d2;
  emxEnsureCapacity_real_T(res1->d, idx_tmp);
  emxFree_int32_T(&aCols);
  for (idx_tmp = 0; idx_tmp < nm1d2; idx_tmp++) {
    res1->d->data[idx_tmp] = 0.0;
  }

  idx_tmp = res1->colidx->size[0];
  res1->colidx->size[0] = n + 1;
  emxEnsureCapacity_int32_T(res1->colidx, idx_tmp);
  res1->colidx->data[0] = 1;
  idx_tmp = res1->rowidx->size[0];
  res1->rowidx->size[0] = nm1d2;
  emxEnsureCapacity_int32_T(res1->rowidx, idx_tmp);
  for (idx_tmp = 0; idx_tmp < nm1d2; idx_tmp++) {
    res1->rowidx->data[idx_tmp] = 0;
  }

  nm1d2 = 0;
  for (loop_ub = 0; loop_ub < n; loop_ub++) {
    idx_tmp = loop_ub + 1;
    while ((nm1d2 + 1 <= mGEn) && (aRows->data[nm1d2] == idx_tmp)) {
      res1->rowidx->data[nm1d2] = ridxInt->data[nm1d2];
      nm1d2++;
    }

    res1->colidx->data[idx_tmp] = nm1d2 + 1;
  }

  emxFree_int32_T(&ridxInt);
  emxFree_int32_T(&aRows);
  for (k = 0; k < mGEn; k++) {
    res1->d->data[k] = aDat->data[sortedIndices->data[k] - 1];
  }

  emxFree_int32_T(&sortedIndices);
  emxFree_real_T(&aDat);
  sparse_fillIn(res1);
}

/* End of code generation (spdiags.cpp) */
