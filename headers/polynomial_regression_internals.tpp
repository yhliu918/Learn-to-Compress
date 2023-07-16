#include <string.h>
#include <cassert>

// Perform polynomial regression using X and Y iterators
template<int n, typename TYPE, typename PRECISION, typename ITERATOR_X, typename ITERATOR_Y>
Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_X x_iter,
                                                     ITERATOR_Y y_iter,
                                                     bool compute_residual,
                                                     size_t N) {
  static_assert(n >= 0);

  // X values raised in degree.
  std::array<std::vector<PRECISION>, 2 * n + 1> x_raised;
  build_x_matrix<2 * n + 1, PRECISION, ITERATOR_X>(x_iter, N, x_raised);
  return polynomial_regression_iter<n, TYPE, PRECISION, ITERATOR_Y>(x_raised, y_iter, compute_residual, N);
}

// Perform polynomial regression Y iterator only (X enumerates 0 to N)
template<int n, typename TYPE, typename PRECISION, typename ITERATOR_Y>
Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_Y y_iter,
                                                     bool compute_residual, size_t N) {
  static_assert(n >= 0);

  // X values raised in degree.
  std::array<std::vector<PRECISION>, 2 * n + 1> x_raised;
  build_x_matrix<2 * n + 1, PRECISION>(N, x_raised);
  return polynomial_regression_iter<n, TYPE, PRECISION, ITERATOR_Y>(x_raised, y_iter, compute_residual, N);
}

// Perform polynomial regression using Y iterator only (X enumerates 0 to N assuming the fixed sample size)
template<int n, int fixed_size, typename TYPE, typename PRECISION, typename ITERATOR_Y>
Polynomial<n, TYPE, PRECISION> polynomial_regression_iter(ITERATOR_Y y_iter, bool compute_residual) {
  static_assert(n >= 0);

  // X values raised in degree, static to compute oly once
  static std::array<std::vector<PRECISION>, 2 * n + 1> x_raised;
  if (x_raised[0].empty()) {
    // Build once on the first call.
    build_x_matrix<2 * n + 1, PRECISION>(fixed_size, x_raised);
  }

  return polynomial_regression_iter<n, TYPE, PRECISION, ITERATOR_Y>(x_raised, y_iter, compute_residual, fixed_size);
}

// Perform polynomial regression over two collections that may have different type but expecting the same size
template<int order, typename TYPE, typename PRECISION,
    typename COLLECTION_X, typename COLLECTION_Y>
Polynomial<order, TYPE, PRECISION> polynomial_regression(const COLLECTION_X &x,
                                                         const COLLECTION_Y &y, bool compute_residual) {
  assert(x.size() == y.size());
  assert(x.size() > 0);
  return polynomial_regression_iter<order, TYPE, PRECISION>(x.cbegin(), y.cbegin(), compute_residual, x.size());
}

template<int order, typename TYPE, typename PRECISION,
    typename COLLECTION_X, typename COLLECTION_Y>
Polynomial<order, TYPE, PRECISION> polynomial_regression(COLLECTION_X &x,
                                                         COLLECTION_Y &y, bool compute_residual, size_t size) {
  return polynomial_regression_iter<order, TYPE, PRECISION>(x.cbegin(), y.cbegin(), compute_residual, size);
}


// Perform polynomial regression over single collection (x simply changes 0 to N)
template<int order, typename TYPE, typename PRECISION, typename COLLECTION_Y>
Polynomial<order, TYPE, PRECISION> polynomial_regression(const COLLECTION_Y &y, bool compute_residual) {
  assert(y.size() > 0);
  return polynomial_regression_iter<order, TYPE, PRECISION>(y.cbegin(), compute_residual, y.size());
}

// Perform polynomial regression over single collection assuming the fixed sample size (x simply changes 0 to N)
// Assuming fixed size allows to compute the x_raised matrix only once.
template<int order, int fixed_size, typename TYPE, typename PRECISION, typename COLLECTION_Y>
Polynomial<order, TYPE, PRECISION> polynomial_regression_fixed(const COLLECTION_Y &y, bool compute_residual) {
  assert(y.size() <= fixed_size);
  return polynomial_regression_iter<order, fixed_size, TYPE, PRECISION>(y.cbegin(), compute_residual);
}

