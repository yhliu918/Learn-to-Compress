#ifndef _POLYNOMIAL_REGRESSION_ABC_INTERNALS_H
#define _POLYNOMIAL_REGRESSION_ABC_INTERNALS_H  __POLYNOMIAL_REGRESSION_ABC_INTERNALS_H

/**
 * PURPOSE:
 *
 *  Polynomial Regression aims to fit a non-linear relationship to a set of
 *  points. It approximates this by solving a series of linear equations using
 *  a least-squares approach.
 *
 *  We can model the expected value y as an nth degree polynomial, yielding
 *  the general polynomial regression model:
 *
 *  y = a0 + a1 * x + a2 * x^2 + ... + an * x^n
 *
 * LICENSE:
 *
 * MIT License
 *
 * Copyright (c) 2020 Chris Engelsma, Audrius Meskauskas
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * @author Chris Engelsma (initial version, all algorithm)
 * @author Audrius Meskauskas (later changes starting from August 19, 2020)
 */

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <array>
#include <string.h>

namespace andviane {

// Build the matrix of X values raised in degree.
  template<int till_degree, typename PRECISION, typename ITERATOR_X>
  static void build_x_matrix(ITERATOR_X x_iter, size_t N, std::array<std::vector<PRECISION>, till_degree> &x_raised) {
    for (int degree = 0; degree < till_degree; degree++) {
      switch (degree) {
        case 0:
          x_raised[degree] = std::vector<PRECISION>(N, 1); // x^0
          break;
        case 1:
          x_raised[degree] = std::vector<PRECISION>(N);   // x^1
          for (int ix = 0; ix < N; ix++) {
            x_raised[degree][ix] = (PRECISION) *x_iter;
            ++x_iter;
          }
          break;
        default:
          x_raised[degree] = std::vector<PRECISION>(N);   // x^degree
          for (int ix = 0; ix < N; ix++)
            x_raised[degree][ix] = x_raised[degree - 1][ix] * x_raised[1][ix];
          break;
      }
    }
  }

// Build the fixed (enumerating) matrix of X values raised in degree.
  template<int till_degree, typename PRECISION>
  static void build_x_matrix(size_t N, std::array<std::vector<PRECISION>, till_degree> &x_raised) {
    for (int degree = 0; degree < till_degree; degree++) {
      switch (degree) {
        case 0:
          x_raised[degree] = std::vector<PRECISION>(N, 1); // x^0
          break;
        case 1:
          x_raised[degree] = std::vector<PRECISION>(N);   // x^1
          for (int ix = 0; ix < N; ix++) {
            x_raised[degree][ix] = (PRECISION) ix;
          }
          break;
        default:
          x_raised[degree] = std::vector<PRECISION>(N);   // x^degree
          for (int ix = 0; ix < N; ix++)
            x_raised[degree][ix] = x_raised[degree - 1][ix] * x_raised[1][ix];
          break;
      }
    }
  }

// Main algorithm of polynomial regression
  template<int n, typename TYPE=double, typename PRECISION=TYPE, typename ITERATOR_Y>
  Polynomial <n, TYPE, PRECISION> polynomial_regression_iter(const std::array<std::vector<PRECISION>, 2 * n + 1> &x_raised,
                                                        ITERATOR_Y y_iter,
                                                        bool compute_residual, size_t N) {
    constexpr int np1 = n + 1;
    constexpr int np2 = n + 2;
    constexpr int tnp1 = 2 * n + 1;

    // X = vector that stores values of sigma(xi^2n)
    PRECISION X[tnp1];
    for (int i = 0; i < tnp1; ++i) {
      X[i] = 0;
      for (int j = 0; j < N; ++j)
        X[i] += x_raised[i][j];
    }

    // a = vector to store final coefficients.
    Polynomial<n, TYPE, PRECISION> a(N);

    // B = normal augmented matrix that stores the equations.
    PRECISION B[np1][np2];
    memset(&B, 0, sizeof(B));

    for (int i = 0; i <= n; ++i)
      for (int j = 0; j <= n; ++j)
        B[i][j] = X[i + j];

    // Y = vector to store values of sigma(xi^n * yi)
    PRECISION Y[np1];
    for (int i = 0; i < np1; ++i) {
      Y[i] = 0;
      ITERATOR_Y y_iter_iter = y_iter;
      for (int j = 0; j < N; ++j) {
        Y[i] += x_raised[i][j] * (*y_iter_iter++);
      }
    }

    // Load values of Y as last column of B
    for (int i = 0; i <= n; ++i)
      B[i][np1] = Y[i];

    // Pivotisation of the B matrix.
    for (int i = 0; i < np1; ++i)
      for (int k = i + 1; k < np1; ++k)
        if (B[i][i] < B[k][i])
          for (int j = 0; j <= np1; ++j) {
            std::swap(B[i][j], B[k][j]);
          }

    // Performs the Gaussian elimination.
    // (1) Make all elements below the pivot equals to zero
    //     or eliminate the variable.
    for (int i = 0; i < n; ++i)
      for (int k = i + 1; k < np1; ++k) {
        PRECISION t = B[k][i] / B[i][i];
        for (int j = 0; j <= np1; ++j)
          B[k][j] -= t * B[i][j];         // (1)
      }

    // Back substitution.
    // (1) Set the variable as the rhs of last equation
    // (2) Subtract all lhs values except the target coefficient.
    // (3) Divide rhs by coefficient of variable being calculated.
    for (int i = n; i >= 0; --i) {
      a[i] = B[i][np1];                 // (1)
      for (int j = 0; j < np1; ++j)
        if (j != i)
          a[i] -= B[i][j] * a[j];       // (2)
      a[i] /= B[i][i];                  // (3)
    }

    if (compute_residual) {
      PRECISION r = 0;
      ITERATOR_Y y_iter_iter = y_iter;
      for (int i = 0; i < N; i++) {
        PRECISION x = x_raised[1][i];
        PRECISION diff = a(x) - (*y_iter_iter++);
        r = r + diff * diff;
      }
      a.residual(r);
    }
    return a;
  }
}
#endif
