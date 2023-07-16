#include <cassert>

template<int p_order, typename TYPE, typename PRECISION>
Polynomial<p_order, TYPE, PRECISION>::Polynomial(int n, bool valid) : valid_(valid), data_size_(n) {
  static_assert(p_order >= 0);
  coefficients_.fill(0);
}

template<int p_order, typename TYPE, typename PRECISION>
Polynomial<p_order, TYPE, PRECISION>::Polynomial(std::array<PRECISION, p_order + 1> coefficients, bool valid, int n) :
    valid_(valid), data_size_(n) {
  static_assert(p_order >= 0);
  coefficients_ = coefficients;
}

// Override (), allowing to use polynomial as function that interpolates
template<int p_order, typename TYPE, typename PRECISION>
TYPE Polynomial<p_order, TYPE, PRECISION>::operator()(TYPE x) {
  PRECISION xx = 1;
  PRECISION s = 0;
  for (PRECISION a: coefficients_) {
    s = s + xx * a;
    xx = xx * x;
  }

  // If the "official type" happens to be integer or the like, we need a proper rounding.
  return std::is_integral<TYPE>::value ? (TYPE) std::round((double) s) : (TYPE) s;
}

template<int p_order, typename TYPE, typename PRECISION>
Polynomial<p_order - 1, TYPE, PRECISION> Polynomial<p_order, TYPE, PRECISION>::differentiate() {
  static_assert(p_order > 0);
  Polynomial<p_order - 1, TYPE, PRECISION> diff;
  for (int n = 1; n <= p_order; n++) {
    diff[n - 1] = n * coefficients_[n];
  }
  return diff;
}

template<int p_order, typename TYPE, typename PRECISION>
Polynomial<p_order + 1, TYPE, PRECISION> Polynomial<p_order, TYPE, PRECISION>::integrate(TYPE C) {
  Polynomial<p_order + 1, TYPE, PRECISION> integ;
  for (int n = 1; n <= p_order + 1; n++) {
    integ[n] = coefficients_[n - 1] / (PRECISION) n;
  }
  integ[0] = C;
  return integ;
}

// Define [] to retrieve the coefficients
template<int p_order, typename TYPE, typename PRECISION>
PRECISION &Polynomial<p_order, TYPE, PRECISION>::operator[](int a) {
  return coefficients_.at(a);
}

// Define the iterators for easy loop
template<int p_order, typename TYPE, typename PRECISION>
typename std::array<PRECISION, p_order>::iterator Polynomial<p_order, TYPE, PRECISION>::begin() {
  return coefficients_.begin();
}

template<int p_order, typename TYPE, typename PRECISION>
typename std::array<PRECISION, p_order>::iterator Polynomial<p_order, TYPE, PRECISION>::end() {
  return coefficients_.end();
}

template<int p_order, typename TYPE, typename PRECISION>
int Polynomial<p_order, TYPE, PRECISION>::order() const {
  return p_order;
}

template<int p_order, typename TYPE, typename PRECISION>
int Polynomial<p_order, TYPE, PRECISION>::data_size() const {
  return data_size_;
}

template<int p_order, typename TYPE, typename PRECISION>
PRECISION Polynomial<p_order, TYPE, PRECISION>::residual() const {
  return residual_;
}

template<int p_order, typename TYPE, typename PRECISION>
void Polynomial<p_order, TYPE, PRECISION>::residual(PRECISION residual) {
  residual_ = residual;
}

template<int p_order, typename TYPE, typename PRECISION>
std::string Polynomial<p_order, TYPE, PRECISION>::DebugString() const {
  std::string expression;
  for (int n = p_order; n >= 0; n--) {
    TYPE k = coefficients_.at(n);
    switch (n) {
      case 0:
        expression = expression + std::to_string((float) k);
        break;
      case 1:
        expression = expression + std::to_string((float) k) + " * x + ";
        break;
      default:
        expression = expression + std::to_string((float) k) + " * x^" + std::to_string(n) + " + ";
        break;
    }
  }
  return expression;
}

template<int p_order, typename TYPE, typename PRECISION>
PRECISION Polynomial<p_order, TYPE, PRECISION>::avg_sqdif() const {
  return data_size_ > 0 ? residual_ / data_size_ : NAN;
}
