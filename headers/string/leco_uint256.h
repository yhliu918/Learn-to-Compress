#ifndef LECO_UINT256_H
#define LECO_UINT256_H

#include <type_traits>
const __uint128_t leco_uint128_0(0);
const __uint128_t leco_uint128_1(1);
const __uint128_t leco_uint128_64(64);
const __uint128_t leco_uint128_128(128);
const __uint128_t leco_uint128_256(256);

namespace std
{
  template <>
  struct is_arithmetic<__uint128_t> : std::true_type
  {
  };
  template <>
  struct is_integral<__uint128_t> : std::true_type
  {
  };
  template <>
  struct is_unsigned<__uint128_t> : std::true_type
  {
  };
}



inline uint32_t bits(__uint128_t v)
{
  uint32_t r(0);
  int length = sizeof(__uint128_t) * 8;
  if (length > 127 && v >= ((__uint128_t)1 << (uint8_t)127))
  {
    v >>= 128;
    r += 128;
  }
  if (length > 63 && v >= ((__uint128_t)1 << (uint8_t)63))
  {
    v >>= 64;
    r += 64;
  }
  if (length > 31 && v >= ((__uint128_t)1 << (uint8_t)31))
  {
    v >>= 32;
    r += 32;
  }
  if (length > 15 && v >= ((__uint128_t)1 << (uint8_t)15))
  {
    v >>= 16;
    r += 16;
  }
  if (length > 7 && v >= ((__uint128_t)1 << (uint8_t)7))
  {
    v >>= 8;
    r += 8;
  }
  if (length > 3 && v >= ((__uint128_t)1 << (uint8_t)3))
  {
    v >>= 4;
    r += 4;
  }
  if (length > 1 && v >= ((__uint128_t)1 << (uint8_t)1))
  {
    v >>= 2;
    r += 2;
  }
  if (v >= (__uint128_t)1)
  {
    r += 1;
  }

  return r;
}


class leco_uint256
{
public:
  // Constructor
  // leco_uint256() = default;
  leco_uint256(const leco_uint256& rhs) = default;
  leco_uint256(leco_uint256&& rhs) = default;
  leco_uint256()
  {
    UPPER = 0;
    LOWER = 0;
  }
  leco_uint256(__uint128_t upper, __uint128_t lower)
  {
    UPPER = upper;
    LOWER = lower;
  }
  uint32_t bit() const
  {
    return bits(UPPER) + bits(LOWER);
  }
  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256(const T& rhs) : LOWER(rhs), UPPER(0)
  {
    if (std::is_signed<T>::value)
    {
      if (rhs < 0)
      {
        UPPER = -1;
      }
    }
  }

  template <
    typename S, typename T,
    typename = typename std::enable_if<
    std::is_integral<S>::value&& std::is_integral<T>::value, void>::type>
  leco_uint256(const S& upper_rhs, const T& lower_rhs)
    : LOWER(lower_rhs), UPPER(upper_rhs) {}

  // Comparison Operators
  bool operator==(const __uint128_t& rhs) const;
  bool operator==(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator==(const T& rhs) const
  {
    return (!UPPER && (LOWER == __uint128_t(rhs)));
  }

  bool operator!=(const __uint128_t& rhs) const;
  bool operator!=(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator!=(const T& rhs) const
  {
    return ((bool)UPPER | (LOWER != __uint128_t(rhs)));
  }

  bool operator>(const __uint128_t& rhs) const;
  bool operator>(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator>(const T& rhs) const
  {
    return ((bool)UPPER | (LOWER > __uint128_t(rhs)));
  }

  bool operator<(const __uint128_t& rhs) const;
  bool operator<(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator<(const T& rhs) const
  {
    return (!UPPER) ? (LOWER < __uint128_t(rhs)) : false;
  }

  bool operator>=(const __uint128_t& rhs) const;
  bool operator>=(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator>=(const T& rhs) const
  {
    return ((*this > rhs) | (*this == rhs));
  }

  bool operator<=(const __uint128_t& rhs) const;
  bool operator<=(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  bool operator<=(const T& rhs) const
  {
    return ((*this < rhs) | (*this == rhs));
  }

  // Assignment Operator
  leco_uint256& operator=(const leco_uint256& rhs) = default;
  leco_uint256& operator=(leco_uint256&& rhs) = default;
  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator=(const T& rhs)
  {
    UPPER = leco_uint128_0;

    if (std::is_signed<T>::value)
    {
      if (rhs < 0)
      {
        UPPER = -1;
      }
    }

    LOWER = rhs;
    return *this;
  }
  leco_uint256& operator=(const bool& rhs)
  {
    UPPER = 0;
    LOWER = rhs;
    return *this;
  }

  // Arithmetic Operators
  leco_uint256 operator+(const leco_uint256& rhs) const
  {
    return leco_uint256(UPPER + rhs.UPPER + ((LOWER + rhs.LOWER) < LOWER),
      LOWER + rhs.LOWER);
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator+(const T& rhs) const
  {
    return leco_uint256(UPPER + ((LOWER + (uint64_t)rhs) < LOWER),
      LOWER + (uint64_t)rhs);
  }

  leco_uint256& operator+=(const leco_uint256& rhs)
  {
    UPPER += rhs.UPPER + ((LOWER + rhs.LOWER) < LOWER);
    LOWER += rhs.LOWER;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator+=(const T& rhs)
  {
    return *this += leco_uint256(rhs);
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator-(const T& rhs) const
  {
    return leco_uint256(UPPER - ((LOWER - rhs) > LOWER), LOWER - rhs);
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator-=(const T& rhs)
  {
    return *this = *this - leco_uint256(rhs);
  }

  leco_uint256 operator-(const __uint128_t& rhs) const
  {
    return *this - leco_uint256(rhs);
  }

  leco_uint256 operator-(const leco_uint256& rhs) const
  {
    return leco_uint256(UPPER - rhs.UPPER - ((LOWER - rhs.LOWER) > LOWER),
      LOWER - rhs.LOWER);
  }

  leco_uint256& operator-=(const __uint128_t& rhs)
  {
    return *this -= leco_uint256(rhs);
  }

  leco_uint256& operator-=(const leco_uint256& rhs)
  {
    *this = *this - rhs;
    return *this;
  }

  std::pair <leco_uint256, leco_uint256> divmod(const leco_uint256& lhs, const leco_uint256& rhs) const {
    // Save some calculations /////////////////////
    if (rhs == 0) {
      throw std::domain_error("Error: division or modulus by 0");
    }
    else if (rhs == 1) {
      return std::pair <leco_uint256, leco_uint256>(lhs, 0);
    }
    else if (lhs == rhs) {
      return std::pair <leco_uint256, leco_uint256>(1, 0);
    }
    else if ((lhs == 0) || (lhs < rhs)) {
      return std::pair <leco_uint256, leco_uint256>(0, lhs);
    }

    std::pair <leco_uint256, leco_uint256> qr(0, lhs);
    leco_uint256 copyd = rhs << (lhs.bit() - rhs.bit());
    leco_uint256 adder = 1 << (lhs.bit() - rhs.bit());
    if (copyd > qr.second) {
      copyd >>= 1;
      adder >>= 1;
    }
    while (qr.second >= rhs) {
      if (qr.second >= copyd) {
        qr.second -= copyd;
        qr.first |= adder;
      }
      copyd >>= 1;
      adder >>= 1;
    }
    return qr;
  }

  leco_uint256 operator/(const uint128_t& rhs) const {
    return *this / leco_uint256(rhs);
  }

  leco_uint256 operator/(const leco_uint256& rhs) const {
    return divmod(*this, rhs).first;
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator/(const T& rhs) const
  {
    return *this / leco_uint256(rhs);
  }

  leco_uint256 operator*(const leco_uint256& rhs) const
  {
    // split values into 4 64-bit parts
    __uint128_t top[4] = { UPPER >> 64, UPPER & 0xffffffffffffffffUL,
                          LOWER >> 64, LOWER & 0xffffffffffffffffUL };
    __uint128_t bottom[4] = { rhs.UPPER >> 64, rhs.UPPER & 0xffffffffffffffffUL,
                             rhs.LOWER >> 64, rhs.LOWER & 0xffffffffffffffffUL };
    __uint128_t products[4][4];

    // multiply each component of the values
    for (int y = 3; y > -1; y--)
    {
      for (int x = 3; x > -1; x--)
      {
        products[3 - y][x] = top[x] * bottom[y];
      }
    }

    // first row
    __uint128_t fourth64 = __uint128_t(products[0][3] & 0xffffffffffffffffUL);
    __uint128_t third64 = __uint128_t(products[0][2] & 0xffffffffffffffffUL) +
      __uint128_t(products[0][3] >> 64);
    __uint128_t second64 = __uint128_t(products[0][1] & 0xffffffffffffffffUL) +
      __uint128_t(products[0][2] >> 64);
    __uint128_t first64 = __uint128_t(products[0][0] & 0xffffffffffffffffUL) +
      __uint128_t(products[0][1] >> 64);

    // second row
    third64 += __uint128_t(products[1][3] & 0xffffffffffffffffUL);
    second64 += __uint128_t(products[1][2] & 0xffffffffffffffffUL) +
      __uint128_t(products[1][3] >> 64);
    first64 += __uint128_t(products[1][1] & 0xffffffffffffffffUL) +
      __uint128_t(products[1][2] >> 64);

    // third row
    second64 += __uint128_t(products[2][3] & 0xffffffffffffffffUL);
    first64 += __uint128_t(products[2][2] & 0xffffffffffffffffUL) +
      __uint128_t(products[2][3] >> 64);

    // fourth row
    first64 += __uint128_t(products[3][3] & 0xffffffffffffffffUL);

    // combines the values, taking care of carry over
    return leco_uint256(first64 << leco_uint128_64, leco_uint128_0) +
      leco_uint256(third64 >> 64, third64 << leco_uint128_64) +
      leco_uint256(second64, leco_uint128_0) + leco_uint256(fourth64);
  }

  leco_uint256& operator*=(const leco_uint256& rhs)
  {
    *this = *this * rhs;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator*(const T& rhs) const
  {
    return *this * leco_uint256(rhs);
  }

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator*=(const T& rhs)
  {
    return *this = *this * leco_uint256(rhs);
  }

  // Bitwise Operators
  leco_uint256 operator&(const __uint128_t& rhs) const;
  leco_uint256 operator&(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator&(const T& rhs) const
  {
    return leco_uint256(leco_uint128_0, LOWER & (__uint128_t)rhs);
  }

  leco_uint256& operator&=(const __uint128_t& rhs);
  leco_uint256& operator&=(const leco_uint256& rhs);

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator&=(const T& rhs)
  {
    UPPER = leco_uint128_0;
    LOWER &= rhs;
    return *this;
  }

  leco_uint256 operator|(const __uint128_t& rhs) const;
  leco_uint256 operator|(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator|(const T& rhs) const
  {
    return leco_uint256(UPPER, LOWER | __uint128_t(rhs));
  }

  leco_uint256& operator|=(const __uint128_t& rhs);
  leco_uint256& operator|=(const leco_uint256& rhs);

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator|=(const T& rhs)
  {
    LOWER |= (__uint128_t)rhs;
    return *this;
  }

  leco_uint256 operator^(const __uint128_t& rhs) const;
  leco_uint256 operator^(const leco_uint256& rhs) const;

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator^(const T& rhs) const
  {
    return leco_uint256(UPPER, LOWER ^ (__uint128_t)rhs);
  }

  leco_uint256& operator^=(const __uint128_t& rhs);
  leco_uint256& operator^=(const leco_uint256& rhs);

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator^=(const T& rhs)
  {
    LOWER ^= (__uint128_t)rhs;
    return *this;
  }

  leco_uint256 operator~() const;

  // Bit Shift Operators
  leco_uint256 operator<<(const __uint128_t& shift) const;
  leco_uint256 operator<<(const leco_uint256& shift) const;

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator<<(const T& rhs) const
  {
    return *this << leco_uint256(rhs);
  }

  leco_uint256& operator<<=(const __uint128_t& shift);
  leco_uint256& operator<<=(const leco_uint256& shift);

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator<<=(const T& rhs)
  {
    *this = *this << leco_uint256(rhs);
    return *this;
  }

  leco_uint256 operator>>(const __uint128_t& shift) const;
  leco_uint256 operator>>(const leco_uint256& shift) const;

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256 operator>>(const T& rhs) const
  {
    return *this >> leco_uint256(rhs);
  }

  leco_uint256& operator>>=(const __uint128_t& shift);
  leco_uint256& operator>>=(const leco_uint256& shift);

  template <typename T, typename = typename std::enable_if<
    std::is_integral<T>::value, T>::type>
  leco_uint256& operator>>=(const T& rhs)
  {
    *this = *this >> leco_uint256(rhs);
    return *this;
  }



  // Typecast operator
  operator bool() const { return (bool)(UPPER | LOWER); }

  operator uint8_t() const { return (uint8_t)LOWER; }

  operator uint16_t() const { return (uint16_t)LOWER; }

  operator uint32_t() const { return (uint32_t)LOWER; }

  operator uint64_t() const { return (uint64_t)LOWER; }

  // private:
  __uint128_t LOWER, UPPER;
};

// .cpp
const leco_uint256 leco_uint256_0 = 0;
const leco_uint256 leco_uint256_1 = 1;
// namespace std
// { // This is probably not a good idea
//   template <>
//   struct is_arithmetic<leco_uint256> : std::true_type
//   {
//   };
//   template <>
//   struct is_integral<leco_uint256> : std::true_type
//   {
//   };
//   template <>
//   struct is_unsigned<leco_uint256> : std::true_type
//   {
//   };
// } // namespace std

// comparison operators
bool leco_uint256::operator==(const __uint128_t& rhs) const
{
  return (*this == leco_uint256(rhs));
}

bool leco_uint256::operator==(const leco_uint256& rhs) const
{
  return ((UPPER == rhs.UPPER) && (LOWER == rhs.LOWER));
}

bool leco_uint256::operator!=(const __uint128_t& rhs) const
{
  return (*this != leco_uint256(rhs));
}

bool leco_uint256::operator!=(const leco_uint256& rhs) const
{
  return ((UPPER != rhs.UPPER) | (LOWER != rhs.LOWER));
}

bool leco_uint256::operator>(const __uint128_t& rhs) const
{
  return (*this > leco_uint256(rhs));
}

bool leco_uint256::operator>(const leco_uint256& rhs) const
{
  if (UPPER == rhs.UPPER)
  {
    return (LOWER > rhs.LOWER);
  }
  if (UPPER > rhs.UPPER)
  {
    return true;
  }
  return false;
}

bool leco_uint256::operator<(const __uint128_t& rhs) const
{
  return (*this < leco_uint256(rhs));
}

bool leco_uint256::operator<(const leco_uint256& rhs) const
{
  if (UPPER == rhs.UPPER)
  {
    return (LOWER < rhs.LOWER);
  }
  if (UPPER < rhs.UPPER)
  {
    return true;
  }
  return false;
}

bool leco_uint256::operator>=(const __uint128_t& rhs) const
{
  return (*this >= leco_uint256(rhs));
}

bool leco_uint256::operator>=(const leco_uint256& rhs) const
{
  return ((*this > rhs) | (*this == rhs));
}

bool leco_uint256::operator<=(const __uint128_t& rhs) const
{
  return (*this <= leco_uint256(rhs));
}

bool leco_uint256::operator<=(const leco_uint256& rhs) const
{
  return ((*this < rhs) | (*this == rhs));
}

leco_uint256 leco_uint256::operator&(const __uint128_t& rhs) const
{
  return leco_uint256(leco_uint128_0, LOWER & rhs);
}

leco_uint256 leco_uint256::operator&(const leco_uint256& rhs) const
{
  return leco_uint256(UPPER & rhs.UPPER, LOWER & rhs.LOWER);
}

leco_uint256& leco_uint256::operator&=(const __uint128_t& rhs)
{
  UPPER = leco_uint128_0;
  LOWER &= rhs;
  return *this;
}

leco_uint256& leco_uint256::operator&=(const leco_uint256& rhs)
{
  UPPER &= rhs.UPPER;
  LOWER &= rhs.LOWER;
  return *this;
}

leco_uint256 leco_uint256::operator|(const __uint128_t& rhs) const
{
  return leco_uint256(UPPER, LOWER | rhs);
}

leco_uint256 leco_uint256::operator|(const leco_uint256& rhs) const
{
  return leco_uint256(UPPER | rhs.UPPER, LOWER | rhs.LOWER);
}

leco_uint256& leco_uint256::operator|=(const __uint128_t& rhs)
{
  LOWER |= rhs;
  return *this;
}

leco_uint256& leco_uint256::operator|=(const leco_uint256& rhs)
{
  UPPER |= rhs.UPPER;
  LOWER |= rhs.LOWER;
  return *this;
}

leco_uint256 leco_uint256::operator^(const __uint128_t& rhs) const
{
  return leco_uint256(UPPER, LOWER ^ rhs);
}

leco_uint256 leco_uint256::operator^(const leco_uint256& rhs) const
{
  return leco_uint256(UPPER ^ rhs.UPPER, LOWER ^ rhs.LOWER);
}

leco_uint256& leco_uint256::operator^=(const __uint128_t& rhs)
{
  LOWER ^= rhs;
  return *this;
}

leco_uint256& leco_uint256::operator^=(const leco_uint256& rhs)
{
  UPPER ^= rhs.UPPER;
  LOWER ^= rhs.LOWER;
  return *this;
}

leco_uint256 leco_uint256::operator~() const
{
  return leco_uint256(~UPPER, ~LOWER);
}

leco_uint256 leco_uint256::operator<<(const __uint128_t& rhs) const
{
  return *this << leco_uint256(rhs);
}

leco_uint256 leco_uint256::operator<<(const leco_uint256& rhs) const
{
  const __uint128_t shift = rhs.LOWER;
  if (((bool)rhs.UPPER) || (shift >= leco_uint128_256))
  {
    return leco_uint256_0;
  }
  else if (shift == leco_uint128_128)
  {
    return leco_uint256(LOWER, leco_uint128_0);
  }
  else if (shift == leco_uint128_0)
  {
    return *this;
  }
  else if (shift < leco_uint128_128)
  {
    return leco_uint256(
      (UPPER << shift) + (LOWER >> (leco_uint128_128 - shift)),
      LOWER << shift);
  }
  else if ((leco_uint128_256 > shift) && (shift > leco_uint128_128))
  {
    return leco_uint256(LOWER << (shift - leco_uint128_128), leco_uint128_0);
  }
  else
  {
    return leco_uint256_0;
  }
}

leco_uint256& leco_uint256::operator<<=(const __uint128_t& shift)
{
  return *this <<= leco_uint256(shift);
}

leco_uint256& leco_uint256::operator<<=(const leco_uint256& shift)
{
  *this = *this << shift;
  return *this;
}

leco_uint256 leco_uint256::operator>>(const __uint128_t& rhs) const
{
  return *this >> leco_uint256(rhs);
}

leco_uint256 leco_uint256::operator>>(const leco_uint256& rhs) const
{
  const __uint128_t shift = rhs.LOWER;
  if (((bool)rhs.UPPER) | (shift >= leco_uint128_256))
  {
    return leco_uint256_0;
  }
  else if (shift == leco_uint128_128)
  {
    return leco_uint256(UPPER);
  }
  else if (shift == leco_uint128_0)
  {
    return *this;
  }
  else if (shift < leco_uint128_128)
  {
    return leco_uint256(UPPER >> shift, (UPPER << (leco_uint128_128 - shift)) +
      (LOWER >> shift));
  }
  else if ((leco_uint128_256 > shift) && (shift > leco_uint128_128))
  {
    return leco_uint256(UPPER >> (shift - leco_uint128_128));
  }
  else
  {
    return leco_uint256_0;
  }
}

leco_uint256& leco_uint256::operator>>=(const __uint128_t& shift)
{
  return *this >>= leco_uint256(shift);
}

leco_uint256& leco_uint256::operator>>=(const leco_uint256& shift)
{
  *this = *this >> shift;
  return *this;
}


#endif // LECO_UINT256_H