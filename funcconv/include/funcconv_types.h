/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * funcconv_types.h
 *
 * Code generation for function 'funcconv_types'
 *
 */

#ifndef FUNCCONV_TYPES_H
#define FUNCCONV_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#include "cs.h"

/* Type Definitions */
struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

struct cell_wrap_2
{
  emxArray_int32_T *f1;
};

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

struct coder_internal_sparse
{
  emxArray_real_T *d;
  emxArray_int32_T *colidx;
  emxArray_int32_T *rowidx;
  int m;
  int n;
};

struct emxArray_boolean_T
{
  bool *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};

#endif

/* End of code generation (funcconv_types.h) */
