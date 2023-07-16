#ifndef _POLYNOMIAL_REGRESSION_ABC_H
#define _POLYNOMIAL_REGRESSION_ABC_H  __POLYNOMIAL_REGRESSION_ABC_H

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
#include <array>
#include <vector>

#include "Polynomial.hpp"
#include "polynomial_regression_internals.hpp"

namespace andviane {

// Perform polynomial regression over two collections that may have different type but expecting the same size
// This function only works with containers that provide the size operator.
  template<int order, typename TYPE=double, typename PRECISION=TYPE,
      typename COLLECTION_X=std::vector<TYPE>, typename COLLECTION_Y=std::vector<TYPE>>
  Polynomial<order, TYPE, PRECISION> polynomial_regression(const COLLECTION_X &x,
                                                           const COLLECTION_Y &y, bool compute_residual = false);

// Perform polynomial regression using X and Y iterators. This function also works with containers that do not provide
// the size operator (like std::forward_list)
  template<int order, typename TYPE=double, typename PRECISION=TYPE,
      typename COLLECTION_X=std::vector<TYPE>, typename COLLECTION_Y=std::vector<TYPE>>
  Polynomial<order, TYPE, PRECISION> polynomial_regression(COLLECTION_X &x,
                                                           COLLECTION_Y &y,
                                                           bool compute_residual, size_t size);

// Perform polynomial regression over single collection (x simply changes 0 to N)
  template<int order, typename TYPE=double, typename PRECISION=TYPE, typename COLLECTION_Y=std::vector<TYPE>>
  Polynomial<order, TYPE, PRECISION> polynomial_regression(const COLLECTION_Y &y, bool compute_residual = false);

// Perform polynomial regression over single collection assuming the fixed sample size (x simply changes 0 to N)
// Assuming fixed size allows to compute the x_raised matrix only once.
  template<int order, int fixed_size, typename TYPE=double, typename PRECISION=TYPE, typename COLLECTION_Y=std::vector<TYPE>>
  Polynomial<order, TYPE, PRECISION> polynomial_regression_fixed(const COLLECTION_Y &y, bool compute_residual = false);

// Perform polynomial regression using X and Y iterators.
  template<int n, typename TYPE=double, typename PRECISION=TYPE, typename ITERATOR_X, typename ITERATOR_Y>
  Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_X x_iter,
                                                       ITERATOR_Y y_iter,
                                                       size_t N, bool compute_residual = false);

  // Perform polynomial regression Y iterator only (X enumerates 0 to N)
  template<int n, typename TYPE=double, typename PRECISION=TYPE, typename ITERATOR_Y>
  Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_Y y_iter,
                                                       size_t N, bool compute_residual = false);

// Perform polynomial regression using Y iterator only (X enumerates 0 to N assuming the fixed sample size)
  template<int n, int fixed_size, typename TYPE=double, typename PRECISION=TYPE, typename ITERATOR_Y>
  Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_Y y_iter, bool compute_residual = false);

#include "polynomial_regression_internals.tpp"
}
#endif
