/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVMATH_H
#define CRVMATH_H

#include "crv.h"
#include <mthQR.h>

/** \file crvMath.h
    \brief main file for math functions used in crv */

namespace crv {

/** \brief binomial function n!/(i!(n-i)!) */
int binomial(int n, int i);
/** \brief trinomial function n!/(i!j!(n-i-j)!) */
int trinomial(int n, int i, int j);
/** \brief "quadnomial" function n!/(i!j!k!(n-i-j-k)!) */
int quadnomial(int n, int i, int j, int k);

/** \brief faster power for integers */
inline double intpow(const double b, const int e)
{
  switch (e) {
  case 0: return 1.0;
  case 1: return b;
  case 2: return b*b;
  case 3: return b*b*b;
  case 4: return b*b*b*b;
  case 5: return b*b*b*b*b;
  case 6: return b*b*b*b*b*b;
  default:
    return intpow(b, e-6) * intpow(b, 6);
  }
}

/** \brief invert a matrix using QR factorization */
void invertMatrixWithQR(int n, mth::Matrix<double>& A,
    mth::Matrix<double>& Ai);
/** \brief invert a matrix using Pivoting and LU decomposition */
void invertMatrixWithPLU(int n, mth::Matrix<double>& A,
    mth::Matrix<double>& Ai);
}

#endif

