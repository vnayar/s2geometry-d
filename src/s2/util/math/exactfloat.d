// Copyright 2009 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.util.math.exactfloat;

import algorithm = std.algorithm;
import format = std.format;
import math = std.math;
import s2.util.hash.mix;
import std.bigint;
import std.range;
import traits = std.traits;

/**
 * ExactFloat is a multiple-precision floating point type based on the OpenSSL
 * Bignum library.  It has the same interface as the built-in "float" and
 * "double" types, but only supports the subset of operators and intrinsics
 * where it is possible to compute the result exactly.  So for example,
 * ExactFloat supports addition and multiplication but not division (since in
 * general, the quotient of two floating-point numbers cannot be represented
 * exactly).  Exact arithmetic is useful for geometric algorithms, especially
 * for disambiguating cases where ordinary double-precision arithmetic yields
 * an uncertain result.
 *
 * ExactFloat is a subset of the faster and more capable MPFloat class (which
 * is based on the GNU MPFR library).  The main reason to use this class
 * rather than MPFloat is that it is subject to a BSD-style license rather
 * than the much more restrictive LGPL license.
 *
 * It has the following features:
 *
 * - ExactFloat uses the same syntax as the built-in "float" and "double"
 * types, for example: x += 4 + fabs(2*y*y - z*z).  There are a few
 * differences (see below), but the syntax is compatible enough so that
 * ExactFloat can be used as a template argument to templatized classes
 * such as Vector2, VectorN, Matrix3x3, etc.
 *
 * - Results are not rounded; instead, precision is increased so that the
 * result can be represented exactly.  An inexact result is returned only
 * in the case of underflow or overflow (yielding signed zero or infinity
 * respectively), or if the maximum allowed precision is exceeded (yielding
 * NaN).  ExactFloat uses IEEE 754-2008 rules for handling infinities, NaN,
 * rounding to integers, etc.
 *
 * - ExactFloat only supports calculations where the result can be
 * represented exactly.  Therefore it supports intrinsics such as fabs()
 * but not transcendentals such as sin(), sqrt(), etc.
 *
 * Syntax Compatibility with "float" and "double"
 * ----------------------------------------------
 *
 * ExactFloat supports a subset of the operators and intrinsics for the
 * built-in "double" type.  (Thus it supports fabs() but not fabsf(), for
 * example.)  The syntax is different only in the following cases:
 *
 * - Casts and implicit conversions to built-in types (including "bool") are
 * not supported.  So for example, the following will not compile:
 *
 * ExactFloat x = 7.5;
 * double y = x;            // ERROR: use x.ToDouble() instead
 * long z = x;              // ERROR: use x.ToDouble() or lround(trunc(x))
 * q = static_cast<int>(x); // ERROR: use x.ToDouble() or lround(trunc(x))
 * if (x) { ... }           // ERROR: use (x != 0) instead
 *
 * - The glibc floating-point classification macros (fpclassify, isfinite,
 * isnormal, isnan, isinf) are not supported.  Instead there are
 * zero-argument methods:
 *
 * ExactFloat x;
 * if (isnan(x)) { ... }  // ERROR: use (x.is_nan()) instead
 * if (isinf(x)) { ... }  // ERROR: use (x.is_inf()) instead
 *
 * Using ExactFloat with Vector3, etc.
 * -----------------------------------
 *
 * ExactFloat can be used with templatized classes such as Vector2 and Vector3
 * (see "util/math/vector.h"), with the following limitations:
 *
 * - Cast() can be used to convert other vector types to an ExactFloat vector
 * type, but not the other way around.  This is because there are no
 * implicit conversions from ExactFloat to built-in types.  You can work
 * around this by calling an explicit conversion method such as
 * ToDouble().  For example:
 *
 * typedef Vector3<ExactFloat> Vector3_xf;
 * Vector3_xf x;
 * Vector3_d y;
 * x = Vector3_xf::Cast(y);   // This works.
 * y = Vector3_d::Cast(x);    // This doesn't.
 * y = Vector3_d(x[0].ToDouble(), x[1].ToDouble(), x[2].ToDouble()); // OK
 *
 * - IsNaN() is not supported because it calls isnan(), which is defined as a
 * macro in <math.h> and therefore can't easily be overrided.
 *
 * Precision Semantics
 * -------------------
 *
 * Unlike MPFloat, ExactFloat does not allow a maximum precision to be
 * specified (it is always unbounded).  Therefore it does not have any of the
 * corresponding constructors.
 *
 * The current precision of an ExactFloat (i.e., the number of bits in its
 * mantissa) is returned by prec().  The precision is increased as necessary
 * so that the result of every operation can be represented exactly.
 */
struct ExactFloat {
public:
  // The following limits are imposed by OpenSSL.

  /// The maximum exponent supported.  If a value has an exponent larger than
  /// this, it is replaced by infinity (with the appropriate sign).
  static immutable int MAX_EXP = 200*1000*1000;  // About 10**(60 million)

  /// The minimum exponent supported.  If a value has an exponent less than
  /// this, it is replaced by zero (with the appropriate sign).
  static immutable int MIN_EXP = -MAX_EXP;   // About 10**(-60 million)

  /**
   * The maximum number of mantissa bits supported.  If a value has more
   * mantissa bits than this, it is replaced with NaN.  (It is expected that
   * users of this class will never want this much precision.)
   */
  static immutable int MAX_PREC = 64 << 20;  // About 20 million digits

  /**
   * Rounding modes.  kRoundTiesToEven and kRoundTiesAwayFromZero both round
   * to the nearest representable value unless two values are equally close.
   * In that case kRoundTiesToEven rounds to the nearest even value, while
   * kRoundTiesAwayFromZero always rounds away from zero.
   */
  enum RoundingMode {
    ROUND_TIES_TO_EVEN,
    ROUND_TIES_AWAY_FROM_ZERO,
    ROUND_TOWARD_ZERO,
    ROUND_AWAY_FROM_ZERO,
    ROUND_TOWARD_POSITIVE,
    ROUND_TOWARD_NEGATIVE
  }

  /////////////////////////////////////////////////////////////////////////////
  // Constructors

  /// Copy constructor.
  this(in ExactFloat b) nothrow @nogc pure {
    _sign = b._sign;
    _bnExp = b._bnExp;
    _bn = b._bn;
  }

  this(T)(T v) nothrow pure {
    this = v;
  }

  /**
   * Construct an ExactFloat from a "double".  The constructor is implicit so
   * that this class can be used as a replacement for "float" or "double" in
   * templatized libraries.  (With an explicit constructor, code such as
   * "ExactFloat f = 2.5;" would not compile.)  All double-precision values are
   * supported, including denormalized numbers, infinities, and NaNs.
   */
  void opAssign(T)(T v) nothrow pure
  if (traits.isFloatingPoint!T) {
    _bn = BigInt();
    _sign = math.signbit(v) ? -1 : 1;
    if (math.isNaN(v)) {
      setNan();
    } else if (math.isInfinity(v)) {
      setInf(_sign);
    } else {
      // The following code is much simpler than messing about with bit masks,
      // has the advantage of handling denormalized numbers and zero correctly,
      // and is actually quite efficient (at least compared to the rest of this
      // code).  "f" is a fraction in the range [0.5, 1), so if we shift it left
      // by the number of mantissa bits in a double (53, including the leading
      // "1") then the result is always an integer.
      int exp;
      T f = math.frexp(math.fabs(v), exp);
      _bn = cast(ulong) math.ldexp(f, DOUBLE_MANTISSA_BITS);
      _bnExp = exp - DOUBLE_MANTISSA_BITS;
      canonicalize();
    }
  }

  /**
   * Construct an ExactFloat from an "int".  Note that in general, ints are
   * automatically converted to doubles and so would be handled by the
   * constructor above.  However, the particular argument (0) would be
   * ambiguous; the compiler wouldn't know whether to treat it as a "double" or
   * "const char*" (since 0 is a valid null pointer constant).  Adding an "int"
   * constructor solves this problem.
   *
   * We do not provide constructors for "unsigned", "long", "unsigned long",
   * "long long", or "unsigned long long", since these types are not typically
   * used in floating-point calculations and it is safer to require them to be
   * explicitly cast.
   */
  void opAssign(T)(T v) nothrow pure
  if (traits.isIntegral!T) {
    _sign = (v >= 0) ? 1 : -1;
    _bn = math.abs(v);
    _bnExp = 0;
    canonicalize();
  }

  /**
   * Construct an ExactFloat from a string (such as "1.2e50").  Requires that
   * the value is exactly representable as a floating-point number (so for
   * example, "0.125" is allowed but "0.1" is not).
   */
  void opAssign(T)(T s) nothrow pure
  if (traits.isSomeString!T) {
    ExactFloat.unimplemented();
  }

  /////////////////////////////////////////////////////////////////////
  // Constants
  //
  // As an alternative to the constants below, you can also just use the
  // constants defined in <math.h>, for example:
  //
  //   ExactFloat x = NAN, y = -INFINITY;

  /// Return an ExactFloat equal to positive zero (if sign >= 0) or
  /// negative zero (if sign < 0).
  static ExactFloat signedZero(int sign) nothrow pure {
    ExactFloat r;
    r.setZero(sign);
    return r;
  }

  /// Return an ExactFloat equal to positive infinity (if sign >= 0) or
  /// negative infinity (if sign < 0).
  static ExactFloat infinity(int sign) nothrow pure {
    ExactFloat r;
    r.setInf(sign);
    return r;
  }

  /// Return an ExactFloat that is NaN (Not-a-Number).
  static ExactFloat nan() nothrow pure {
    ExactFloat r;
    r.setNan();
    return r;
  }


  /////////////////////////////////////////////////////////////////////////////
  // Accessor Methods

  /// Return the maximum precision of the ExactFloat.  This method exists only
  /// for compatibility with MPFloat.
  @property
  int maxPrec() const nothrow @nogc pure {
    return MAX_PREC;
  }

  /// Return the actual precision of this ExactFloat (the current number of
  /// bits in its mantissa).  Returns 0 for non-normal numbers such as NaN.
  @property
  int prec() const nothrow @nogc pure {
    // TODO: BUG Fix this to get the exact number of used bits.
    int totalBits = cast(int) (_bn.uintLength() - 1) * 32;
    uint lastDigit = _bn.getDigit!uint(_bn.uintLength() - 1);
    while (lastDigit != 0) {
      lastDigit >>= 1;
      totalBits++;
    }
    return totalBits;
  }

  /**
   * Return the exponent of this ExactFloat given that the mantissa is in the
   * range [0.5, 1).  It is an error to call this method if the value is zero,
   * infinity, or NaN.
   */
  int exp() const nothrow @nogc pure
  in {
    assert(isNormal());
  } do {
    return _bnExp + prec();
  }

  /// Set the value of the ExactFloat to +0 (if sign >= 0) or -0 (if sign < 0).
  void setZero(int sign) nothrow pure {
    _sign = sign;
    _bnExp = EXP_ZERO;
    if (_bn != 0) {
      _bn = 0;
    }
  }

  /// Set the value of the ExactFloat to positive infinity (if sign >= 0) or
  /// negative infinity (if sign < 0).
  void setInf(int sign) nothrow pure {
    _sign = sign;
    _bnExp = EXP_INFINITY;
    if (_bn != 0) {
      _bn = 0;
    }
  }

  /// Set the value of the ExactFloat to NaN (Not-a-Number).
  void setNan() nothrow pure {
    _sign = 1;
    _bnExp = EXP_NAN;
    if (_bn != 0) {
      _bn = 0;
    }
  }

  // Unfortunately, isinf(x), isnan(x), isnormal(x), and isfinite(x) are
  // defined as macros in <math.h>.  Therefore we can't easily extend them
  // here.  Instead we provide methods with underscores in their names that do
  // the same thing: x.is_inf(), etc.
  //
  // These macros are not implemented: signbit(x), fpclassify(x).

  /// Return true if this value is zero (including negative zero).
  bool isZero() const nothrow @nogc pure {
    return _bnExp == EXP_ZERO;
  }

  /// Return true if this value is infinity (positive or negative).
  bool isInf() const nothrow @nogc pure {
    return _bnExp == EXP_INFINITY;
  }

  /// Return true if this value is NaN (Not-a-Number).
  bool isNan() const nothrow @nogc pure {
    return _bnExp == EXP_NAN;
  }

  /**
   * Return true if this value is a normal floating-point number.  Non-normal
   * values (zero, infinity, and NaN) often need to be handled separately
   * because they are represented using special exponent values and their
   * mantissa is not defined.
   */
  bool isNormal() const nothrow @nogc pure {
    return _bnExp < EXP_ZERO;
  }


  /// Return true if this value is a normal floating-point number or zero,
  /// i.e. it is not infinity or NaN.
  bool isFinite() const nothrow @nogc pure {
    return _bnExp <= EXP_ZERO;
  }


  /// Return true if the sign bit is set (this includes negative zero).
  @property
  bool signBit() const nothrow @nogc pure {
    return _sign < 0;
  }

  /**
   * Return +1 if this ExactFloat is positive, -1 if it is negative, and 0
   * if it is zero or NaN.  Note that unlike sign_bit(), sgn() returns 0 for
   * both positive and negative zero.
   */
  @property
  int sign() const nothrow @nogc pure {
    return (isNan() || isZero()) ? 0 : _sign;
  }


  /////////////////////////////////////////////////////////////////////////////
  // Conversion Methods
  //
  // Note that some conversions are defined as functions further below,
  // e.g. to convert to an integer you can use lround(), llrint(), etc.

  /**
   * Round to double precision.  Note that since doubles have a much smaller
   * exponent range than ExactFloats, very small values may be rounded to
   * (positive or negative) zero, and very large values may be rounded to
   * infinity.
   *
   * It is very important to make this a named method rather than an implicit
   * conversion, because otherwise there would be a silent loss of precision
   * whenever some desired operator or function happens not to be implemented.
   * For example, if fabs() were not implemented and "x" and "y" were
   * ExactFloats, then x = fabs(y) would silently convert "y" to a "double",
   * take its absolute value, and convert it back to an ExactFloat.
   */
  double toDouble() const nothrow pure {
    // If the mantissa has too many bits, we need to round it.
    if (prec() <= DOUBLE_MANTISSA_BITS) {
      return toDoubleHelper();
    } else {
      ExactFloat r = roundToMaxPrec(DOUBLE_MANTISSA_BITS, RoundingMode.ROUND_TIES_TO_EVEN);
      return r.toDoubleHelper();
    }
  }

  /**
   * Return a human-readable string such that if two values with the same
   * precision are different, then their string representations are different.
   * The format is similar to printf("%g"), except that the number of
   * significant digits depends on the precision (with a minimum of 10).
   * Trailing zeros are stripped (just like "%g").
   *
   * Note that if two values have different precisions, they may have the same
   * ToString() value even though their values are slightly different.  If you
   * need to distinguish such values, use ToUniqueString() intead.
   */
  string toString() const {
    int max_digits = algorithm.max(MIN_SIGNIFICANT_DIGITS,
        numSignificantDigitsForPrec(prec()));
    return toStringWithMaxDigits(max_digits);
  }

  /// Return a string formatted according to printf("%Ng") where N is the given
  /// maximum number of significant digits.
  string toStringWithMaxDigits(int max_digits) const
  in {
    assert(max_digits >= 0);
  } do {
    if (!isNormal()) {
      if (isNan()) return "nan";
      if (isZero()) return (_sign < 0) ? "-0" : "0";
      return (_sign < 0) ? "-inf" : "inf";
    }
    char[] digits;
    int exp10 = getDecimalDigits(max_digits, digits);
    string str;
    if (_sign < 0) str ~= "-";

    // We use the standard '%g' formatting rules.  If the exponent is less than
    // -4 or greater than or equal to the requested precision (i.e., max_digits)
    // then we use exponential notation.
    //
    // But since "exp10" is the base-10 exponent corresponding to a mantissa in
    // the range [0.1, 1), whereas the '%g' rules assume a mantissa in the range
    // [1.0, 10), we need to adjust these parameters by 1.
    if (exp10 <= -4 || exp10 > max_digits) {
      // Use exponential format.
      str ~= digits[0];
      if (digits.length > 1) {
        str ~= ".";
        str ~= digits[1 .. $];
      }
      char[20] exp_buf;
      format.sformat(exp_buf, "e%+02d", exp10 - 1);
      str ~= exp_buf;
    } else {
      // Use fixed format.  We split this into two cases depending on whether
      // the integer portion is non-zero or not.
      if (exp10 > 0) {
        if (exp10 >= digits.length) {
          str ~= digits;
          for (ulong i = exp10 - digits.length; i > 0; --i) {
            str ~= "0";
          }
        } else {
          str ~= digits[0 .. exp10];
          str ~= ".";
          str ~= digits[exp10 .. $];
        }
      } else {
        str ~= "0.";
        for (int i = exp10; i < 0; ++i) {
          str ~= "0";
        }
        str ~= digits;
      }
    }
    return str;
  }

  /**
   * Return a human-readable string such that if two ExactFloats have different
   * values, then their string representations are always different.  This
   * method is useful for debugging.  The string has the form "value<prec>",
   * where "prec" is the actual precision of the ExactFloat (e.g., "0.215<50>").
   */
  string toUniqueString() const {
    string precStr = format.format("<%d>", prec());
    return toString() ~ precStr;
  }

  size_t toHash() const nothrow @safe {
    return HashMix(_bn.toHash())
        .mix(cast(size_t) _bnExp)
        .mix(cast(size_t) _sign)
        .get();
  }

  /**
   * Return an upper bound on the number of significant digits required to
   * distinguish any two floating-point numbers with the given precision when
   * they are formatted as decimal strings in exponential format.
   */
  static int numSignificantDigitsForPrec(int prec) nothrow {
    // The simplest bound is
    //
    //    d <= 1 + ceil(prec * log10(2))
    //
    // The following bound is tighter by 0.5 digits on average, but requires
    // the exponent to be known as well:
    //
    //    d <= ceil(exp * log10(2)) - floor((exp - prec) * log10(2))
    //
    // Since either of these bounds can be too large by 0, 1, or 2 digits, we
    // stick with the simpler first bound.
    return cast(int) (1 + math.ceil(prec * (math.LN2 / math.LN10)));
  }

  /////////////////////////////////////////////////////////////////////////////
  // Other Methods

  /**
   * Round the ExactFloat so that its mantissa has at most "max_prec" bits
   * using the given rounding mode.  Requires "max_prec" to be at least 2
   * (since kRoundTiesToEven doesn't make sense with fewer bits than this).
   */
  ExactFloat roundToMaxPrec(int max_prec, RoundingMode mode) const nothrow pure
  in {
    // The "kRoundTiesToEven" mode requires at least 2 bits of precision
    // (otherwise both adjacent representable values may be odd).
    assert(max_prec >= 2);
    assert(max_prec <= MAX_PREC);
  } do {
    // The following test also catches zero, infinity, and NaN.
    int shift = prec() - max_prec;
    if (shift <= 0) {
      return this;
    }

    // Round by removing the appropriate number of bits from the mantissa.  Note
    // that if the value is rounded up to a power of 2, the high-order bit
    // position may increase, but in that case Canonicalize() will remove at
    // least one zero bit and so the output will still have prec() <= max_prec.
   return roundToPowerOf2(_bnExp + shift, mode);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Operators

  /// Unary plus.
  ExactFloat opUnary(string op)() const nothrow @nogc pure
  if (op == "+") {
    return *this;
  }

  /// Unary minus.
  ExactFloat opUnary(string op)() const nothrow @nogc pure
  if (op == "-") {
    return copyWithSign(-_sign);
  }

  /// Addition.
  ExactFloat opBinary(string op)(in ExactFloat b) const nothrow
  if (op == "+") {
    return signedSum(_sign, this, b._sign, b);
  }

  /// Subtraction.
  ExactFloat opBinary(string op)(in ExactFloat b) const nothrow
  if (op == "-") {
    return signedSum(_sign, this, -b._sign, b);
  }

  /// Multiplication.
  ExactFloat opBinary(string op)(in ExactFloat b) const nothrow
  if (op == "*") {
    int result_sign = _sign * b._sign;
    if (!isNormal() || !b.isNormal()) {
      // Handle zero, infinity, and NaN according to IEEE 754-2008.
      if (isNan()) return this;
      if (b.isNan()) return b;
      if (isInf()) {
        // Infinity times zero yields NaN.
        if (b.isZero()) return ExactFloat.nan();
        return ExactFloat.infinity(result_sign);
      }
      if (b.isInf()) {
        if (isZero()) return ExactFloat.nan();
        return ExactFloat.infinity(result_sign);
      }
      assert(isZero() || b.isZero());
      return ExactFloat.signedZero(result_sign);
    }
    ExactFloat r;
    r._sign = result_sign;
    r._bnExp = _bnExp + b._bnExp;
    r._bn = _bn * b._bn;
    r.canonicalize();
    return r;
  }

  /// Support operations with any convertable types.
  ExactFloat opBinary(string op, T)(in T b) const nothrow {
    return opBinary!op(ExactFloat(b));
  }

  /// Support operations with any convertable types.
  ExactFloat opBinaryRight(string op, T)(in T a) const nothrow {
    ExactFloat axf = ExactFloat(a);
    return axf.opBinary!op(this);
  }

  // Division is not implemented because the result cannot be represented
  // exactly in general.  Doing this properly would require extending all the
  // operations to support rounding to a specified precision.

  // Arithmetic assignment operators (+=, -=, *=).
  ref ExactFloat opOpAssign(string op)(in ExactFloat b) nothrow @nogc pure
  if (op == "+" || op == "-" || op == "*") {
    return mixin("this = this " ~ op ~ " b");
  }

  /// Comparison operators (==, !=).
  bool opEquals(in ExactFloat b) const nothrow @nogc pure {
    // NaN is not equal to anything, not even itself.
    if (isNan() || b.isNan()) return false;

    // Since Canonicalize() strips low-order zero bits, all other cases
    // (including non-normal values) require bn_exp_ to be equal.
    if (_bnExp != b._bnExp) return false;

    // Positive and negative zero are equal.
    if (isZero() && b.isZero()) return true;

    // Otherwise, the signs and mantissas must match.  Note that non-normal
    // values such as infinity have a mantissa of zero.
    return _sign == b._sign && _bn == b._bn;
  }

  /// Support operations with any convertable types.
  bool opEquals(T)(in T b) const nothrow @nogc pure {
    return opEquals(ExactFloat(b));
  }

  int scaleAndCompare(in ExactFloat b) const nothrow @nogc pure
  in {
    assert(isNormal() && b.isNormal() && _bnExp >= b._bnExp);
  } do {
    ExactFloat tmp = this;
    tmp._bnExp <<= _bnExp - b._bnExp;
    if (tmp._bn > b._bn) return 1;
    else if (tmp._bn < b._bn) return -1;
    else return 0;
  }

  bool unsignedLess(in ExactFloat b) const nothrow @nogc pure {
    // Handle the zero/infinity cases (NaN has already been done).
    if (isInf() || b.isZero()) return false;
    if (isZero() || b.isInf()) return true;
    // If the high-order bit positions differ, we are done.
    int cmp = exp() - b.exp();
    if (cmp != 0) return cmp < 0;
    // Otherwise shift one of the two values so that they both have the same
    // bn_exp_ and then compare the mantissas.
    return (_bnExp >= b._bnExp ?
        scaleAndCompare(b) < 0 : b.scaleAndCompare(this) > 0);
  }


  /// Comparison operators (<, <=, >, >=).
  int opCmp(in ExactFloat b) const nothrow @nogc pure {
    // NaN is unordered compared to everything, including itself.
    if (isNan() || b.isNan()) return -1;
    // Positive and negative zero are equal.
    if (isZero() && b.isZero()) return 0;
    // Otherwise, anything negative is less than anything positive.
    if (_sign != b._sign) return _sign > b._sign ? 1 : -1;
    // Now we just compare absolute values.
    bool isLess = (_sign > 0) ? unsignedLess(b) : b.unsignedLess(this);
    if (isLess) return -1;
    else if (this == b) return 0;
    return 1;
  }

  /// Support operations with any convertable types.
  int opCmp(T)(in T b) const nothrow pure {
    return opCmp(ExactFloat(b));
  }

  //////// Miscellaneous simple arithmetic functions.

  /// Absolute value.
  ExactFloat fabs() const nothrow @nogc pure {
    return abs();
  }

  ExactFloat abs() const nothrow @nogc pure {
    return copyWithSign(+1);
  }

private:
  // Numbers are always formatted with at least this many significant digits.
  // This prevents small integers from being formatted in exponential notation
  // (e.g. 1024 formatted as 1e+03), and also avoids the confusion of having
  // supposedly "high precision" numbers formatted with just 1 or 2 digits
  // (e.g. 1/512 == 0.001953125 formatted as 0.002).
  static const int MIN_SIGNIFICANT_DIGITS = 10;

  // Non-normal numbers are represented using special exponent values and a
  // mantissa of zero.  Do not change these values; methods such as
  // is_normal() make assumptions about their ordering.  Non-normal numbers
  // can have either a positive or negative sign (including zero and NaN).
  static immutable int EXP_NAN = int.max;
  static immutable int EXP_INFINITY = int.max - 1;
  static immutable int EXP_ZERO = int.max - 2;

  // Normal numbers are represented as (sign_ * bn_ * (2 ** bn_exp_)), where:
  //  - _sign is either +1 or -1
  //  - _bn is a BIGNUM with a positive value
  //  - _bnExp is the base-2 exponent applied to _bn.
  int _sign = 1;
  int _bnExp = EXP_ZERO;
  BigInt _bn = BigInt();

  /// A standard IEEE "double" has a 53-bit mantissa consisting of a 52-bit
  /// fraction plus an implicit leading "1" bit.
  static immutable int DOUBLE_MANTISSA_BITS = 53;

  /**
   * Convert an ExactFloat with no more than 53 bits in its mantissa to a
   * "double".  This method handles non-normal values (NaN, etc).
   */
  double toDoubleHelper() const nothrow @nogc pure
  in {
    assert(prec() <= DOUBLE_MANTISSA_BITS);
  } do {
    if (!isNormal()) {
      if (isZero()) {
        return math.copysign(0, cast(double) _sign);
      }
      if (isInf()) {
        return math.copysign(double.infinity, cast(double) _sign);
      }
      return math.copysign(double.nan, cast(double) _sign);
    }
    long d_mantissa = _bn.toLong();
    // We rely on ldexp() to handle overflow and underflow.  (It will return a
    // signed zero or infinity if the result is too small or too large.)
    return _sign * math.ldexp(cast(double) d_mantissa, _bnExp);
  }

  /**
   * Round an ExactFloat so that it is a multiple of (2 ** bit_exp), using the
   * given rounding mode.
   */
  ExactFloat roundToPowerOf2(int bit_exp, RoundingMode mode) const nothrow pure
  in {
    assert(bit_exp >= MIN_EXP - MAX_PREC);
    assert(bit_exp <= MAX_EXP);
  } do {
    // If the exponent is already large enough, or the value is zero, infinity,
    // or NaN, then there is nothing to do.
    int shift = bit_exp - _bnExp;
    if (shift <= 0) return this;
    assert(isNormal());

    // Convert rounding up/down to toward/away from zero, so that we don't need
    // to consider the sign of the number from this point onward.
    if (mode == RoundingMode.ROUND_TOWARD_POSITIVE) {
      mode = (_sign > 0) ? RoundingMode.ROUND_AWAY_FROM_ZERO : RoundingMode.ROUND_TOWARD_ZERO;
    } else if (mode == RoundingMode.ROUND_TOWARD_NEGATIVE) {
      mode = (_sign > 0) ? RoundingMode.ROUND_TOWARD_ZERO : RoundingMode.ROUND_AWAY_FROM_ZERO;
    }

    // Rounding consists of right-shifting the mantissa by "shift", and then
    // possibly incrementing the result (depending on the rounding mode, the
    // bits that were discarded, and sometimes the lowest kept bit).  The
    // following code figures out whether we need to increment.
    ExactFloat r;
    bool increment = false;
    if (mode == RoundingMode.ROUND_TOWARD_ZERO) {
      // Never increment.
    } else if (mode == RoundingMode.ROUND_TIES_AWAY_FROM_ZERO) {
      // Increment if the highest discarded bit is 1.
      if (isBitSet(_bn, shift - 1))
        increment = true;
    } else if (mode == RoundingMode.ROUND_AWAY_FROM_ZERO) {
      // Increment unless all discarded bits are zero.
      if (extCountLowZeroBits(_bn) < shift)
        increment = true;
    } else {
      assert(mode == RoundingMode.ROUND_TIES_TO_EVEN);
      // Let "w/xyz" denote a mantissa where "w" is the lowest kept bit and
      // "xyz" are the discarded bits.  Then using regexp notation:
      //    ./0.*       ->    Don't increment (fraction < 1/2)
      //    0/10*       ->    Don't increment (fraction = 1/2, kept part even)
      //    1/10*       ->    Increment (fraction = 1/2, kept part odd)
      //    ./1.*1.*    ->    Increment (fraction > 1/2)
      if (isBitSet(_bn, shift - 1)
          && ((isBitSet(_bn, shift)
              || extCountLowZeroBits(_bn) < shift - 1))) {
        increment = true;
      }
    }
    r._bnExp = _bnExp + shift;
    r._bn = _bn >> shift;
    if (increment) {
      r._bn += 1;
    }
    r._sign = _sign;
    r.canonicalize();
    return r;
  }

  private static bool isBitSet(in BigInt bn, int bitNum) nothrow @nogc pure {
    size_t digitNum = bitNum / (8 * ulong.sizeof);
    size_t shift = bitNum % (8 * ulong.sizeof);
    ulong digit = bn.getDigit!ulong(digitNum);
    return (digit & (1uL << shift)) != 0;
  }

  /**
   * Count the number of low-order zero bits in the given BIGNUM (ignoring its
   * sign).  Returns 0 if the argument is zero.
   */
  private static int extCountLowZeroBits(in BigInt bn) nothrow @nogc pure {
    int count = 0;
    for (int i = 0; i < bn.ulongLength(); ++i) {
      ulong w = bn.getDigit!ulong(i);
      if (w == 0) {
        count += 8 * ulong.sizeof;
      } else {
        for (; (w & 1) == 0; w >>= 1) {
          ++count;
        }
        break;
      }
    }
    return count;
  }

  /**
   * Convert the ExactFloat to a decimal value of the form 0.ddd * (10 ** x),
   * with at most "max_digits" significant digits (trailing zeros are removed).
   * Set (*digits) to the ASCII digits and return the decimal exponent "x".
   */
  int getDecimalDigits(int max_digits, out char[] digits) const
  in {
    assert(isNormal());
  } do {
    // Convert the value to the form (bn * (10 ** bn_exp10)) where "bn" is a
    // positive integer (BIGNUM).
    BigInt bn = BigInt();
    int bn_exp10;
    if (_bnExp >= 0) {
      // The easy case: bn = bn_ * (2 ** bn_exp_)), bn_exp10 = 0.
      bn = _bn << _bnExp;
      bn_exp10 = 0;
    } else {
      // Set bn = bn_ * (5 ** -bn_exp_) and bn_exp10 = bn_exp_.  This is
      // equivalent to the original value of (bn_ * (2 ** bn_exp_)).
      bn = 5;
      bn ^^= -_bnExp;
      bn *= _bn;
      bn_exp10 = _bnExp;
    }
    // Now convert "bn" to a decimal string.
    string all_digits = format.format("%d", bn);
    assert(all_digits != null);
    // Check whether we have too many digits and round if necessary.
    size_t num_digits = all_digits.length;
    if (num_digits <= max_digits) {
      digits = all_digits.dup;
    } else {
      digits = all_digits[0 .. max_digits].dup;
      // Standard "printf" formatting rounds ties to an even number.  This means
      // that we round up (away from zero) if highest discarded digit is '5' or
      // more, unless all other discarded digits are zero in which case we round
      // up only if the lowest kept digit is odd.
      if (all_digits[max_digits] >= '5' &&
          ((all_digits[max_digits-1] & 1) == 1 ||
              !algorithm.findAmong(all_digits[max_digits + 1 .. $], "123456789").empty)) {
        // This can increase the number of digits by 1, but in that case at
        // least one trailing zero will be stripped off below.
        incrementDecimalDigits(digits);
      }
      // Adjust the base-10 exponent to reflect the digits we have removed.
      bn_exp10 += num_digits - max_digits;
    }

    // Now strip any trailing zeros.
    assert(digits[0] != '0');
    size_t pos = digits.length - 1;
    while (digits[pos] == '0') {
      --pos;
    }
    if (pos < digits.length - 1) {
      bn_exp10 += digits.length - pos;
      digits.length = pos + 1;
    }
    assert(digits.length <= max_digits);

    // Finally, we adjust the base-10 exponent so that the mantissa is a
    // fraction in the range [0.1, 1) rather than an integer.
    return bn_exp10 + cast(int) digits.length;
  }

  /// Increment an unsigned integer represented as a string of ASCII digits.
  private static void incrementDecimalDigits(char[] digits) nothrow {
    size_t pos = digits.length;
    while (--pos >= 0) {
      if (digits[pos] < '9') {
        ++digits[pos];
        return;
      }
      digits[pos] = '0';
    }
    digits = "1" ~ digits;
  }

  /// Return a_sign * fabs(a) + b_sign * fabs(b).  Used to implement addition
  /// and subtraction.
  static ExactFloat signedSum(int a_sign, ExactFloat a, int b_sign, ExactFloat b) nothrow {
    if (!a.isNormal() || !b.isNormal()) {
      // Handle zero, infinity, and NaN according to IEEE 754-2008.
      if (a.isNan()) return a;
      if (b.isNan()) return b;
      if (a.isInf()) {
        // Adding two infinities with opposite sign yields NaN.
        if (b.isInf() && a_sign != b_sign) return ExactFloat.nan();
        return infinity(a_sign);
      }
      if (b.isInf()) return infinity(b_sign);
      if (a.isZero()) {
        if (!b.isZero()) return b.copyWithSign(b_sign);
        // Adding two zeros with the same sign preserves the sign.
        if (a_sign == b_sign) return signedZero(a_sign);
        // Adding two zeros of opposite sign produces +0.
        return signedZero(+1);
      }
      assert(b.isZero());
      return a.copyWithSign(a_sign);
    }

    // Swap the numbers if necessary so that "a" has the larger bn_exp_.
    if (a._bnExp < b._bnExp) {
      algorithm.swap(a_sign, b_sign);
      algorithm.swap(a, b);
    }
    // Shift "a" if necessary so that both values have the same bn_exp_.
    ExactFloat r;
    if (a._bnExp > b._bnExp) {
      r._bn = a._bn << (a._bnExp - b._bnExp);
      a = r;  // The only field of "a" used below is bn_.
    }
    r._bnExp = b._bnExp;
    if (a_sign == b_sign) {
      r._bn = a._bn + b._bn;
      r._sign = a_sign;
    } else {
      // Note that the BIGNUM documentation is out of date -- all methods now
      // allow the result to be the same as any input argument, so it is okay if
      // (a == &r) due to the shift above.
      r._bn = a._bn - b._bn;
      if (r._bn == 0) {
        r._sign = +1;
      } else if (r._bn < 0) {
        // The magnitude of "b" was larger.
        r._sign = b_sign;
        r._bn = -r._bn;
      } else {
        // They were equal, or the magnitude of "a" was larger.
        r._sign = a_sign;
      }
    }
    r.canonicalize();
    return r;
  }

  /**
   * Convert an ExactFloat to its canonical form.  Underflow results in signed
   * zero, overflow results in signed infinity, and precision overflow results
   * in NaN.  A zero mantissa is converted to the canonical zero value with
   * the given sign; otherwise the mantissa is normalized so that its low bit
   * is 1.  Non-normal numbers are left unchanged.
   */
  void canonicalize() nothrow pure {
    if (!isNormal()) return;

    // Underflow/overflow occurs if exp() is not in [kMinExp, kMaxExp].
    // We also convert a zero mantissa to signed zero.
    int my_exp = exp();
    if (my_exp < MIN_EXP || _bn == 0) {
      setZero(_sign);
    } else if (my_exp > MAX_EXP) {
      setInf(_sign);
    } else if (_bn % 2 == 1) {
      // Remove any low-order zero bits from the mantissa.
      assert(_bn != 0);
      int shift = countLowZeroBits(_bn);
      if (shift > 0) {
        _bn = _bn >> shift;
        _bnExp += shift;
      }
    }
    // If the mantissa has too many bits, we replace it by NaN to indicate
    // that an inexact calculation has occurred.
    if (prec() > MAX_PREC) {
      setNan();
    }
  }

  /**
   * Count the number of low-order zero bits in the given BIGNUM (ignoring its
   * sign).  Returns 0 if the argument is zero.
   */
  private static int countLowZeroBits(in BigInt bn) nothrow pure {
    int count = 0;
    for (int i = 0; i < bn.ulongLength(); ++i) {
      ulong w = bn.getDigit!ulong(i);
      if (w == 0) {
        count += 8 * ulong.sizeof;
      } else {
        for (; (w & 1) == 0; w >>= 1) {
          ++count;
        }
        break;
      }
    }
    return count;
  }

  /**
   * Return an ExactFloat with the magnitude of this ExactFloat and the given
   * sign.  (Similar to copysign, except that the sign is given explicitly
   * rather than being copied from another ExactFloat.)
   */
  ExactFloat copyWithSign(int sign) const nothrow @nogc pure {
    ExactFloat r = this;
    r._sign = sign;
    return r;
  }

  /**
   * Convert an ExactFloat to an integer of type "T" using the given rounding
   * mode.  The type "T" must be signed.  Returns the largest possible integer
   * for NaN, and clamps out of range values to the largest or smallest
   * possible values.
   */
  T toInteger(T)(RoundingMode mode) const nothrow
  if (traits.isIntegral!T && traits.isSigned!T) {
    ExactFloat r = roundToPowerOf2(0, mode);
    if (r.isNan()) return T.min;
    if (r.isZero()) return 0;
    if (!r.isInf()) {
      // If the unsigned value has more than 63 bits it is always clamped.
      if (r.exp() < 64) {
        long value = r._bn.toLong() << r._bnExp;
        if (r._sign < 0) value = -value;
        return algorithm.max(T.min, algorithm.min(T.max, value));
      }
    }
    return (r._sign < 0) ? T.min : T.max;
  }

  // Log a fatal error message (used for unimplemented methods).
  static ExactFloat unimplemented() nothrow {
    assert(false, "Unimplemented ExactFloat method called.");
  }
}


//////////////////////////////////////////////////////////////////////
// Math Intrinsics
//
// The math intrinsics currently supported by ExactFloat are listed below.
// Except as noted, they behave identically to the usual glibc intrinsics
// except that they have greater precision.  See the man pages for more
// information.
//
// They can be used in a consistent manner with functions used for double
// and real with Uniform Function Call Syntax.
//
// For example:
//   double dval = 5.4;
//   ExactFloat xfval = 5.4;  // Consistent, no problem.
//
//   import math = std.math;
//   math.fabs(dval);
//   math.fabs(xfval);  // Error: Does not exist!
//
//   import std.math;
//   import s2.util.math.exactfloat;
//   fabs(dval);   // OK, but symbol fabs is in the main namespace.
//   dval.fabs();  // OK, uses UFCS.
//   fabs(xfval);  // OK as well, but function comes from exactfloat.
//   xfval.fabs(); // OK, uses UFCS.


//////// Miscellaneous simple arithmetic functions.

/// Absolute value.
ExactFloat fabs(in ExactFloat a) nothrow @nogc pure {
  return abs(a);
}

ExactFloat abs(in ExactFloat a) nothrow @nogc pure {
  return a.copyWithSign(+1);
}

/// Maximum of two values.
ExactFloat fmax(in ExactFloat a, in ExactFloat b) nothrow {
  // If one argument is NaN, return the other argument.
  if (a.isNan()) return b;
  if (b.isNan()) return a;
  // Not required by IEEE 754, but we prefer +0 over -0.
  if (a.sign != b.sign) {
    return (a.sign() < b.sign()) ? b : a;
  }
  return (a < b) ? b : a;
}

/// Minimum of two values.
ExactFloat fmin(in ExactFloat a, in ExactFloat b) nothrow {
  // If one argument is NaN, return the other argument.
  if (a.isNan()) return b;
  if (b.isNan()) return a;
  // Not required by IEEE 754, but we prefer -0 over +0.
  if (a.sign != b.sign) {
    return (a.sign < b.sign) ? a : b;
  }
  return (a < b) ? a : b;
}

/// Positive difference: max(a - b, 0).
ExactFloat fdim(in ExactFloat a, in ExactFloat b) nothrow {
  // This formulation has the correct behavior for NaNs.
  return (a <= b) ? ExactFloat(0) : (a - b);
}

//////// Integer rounding functions that return ExactFloat values.

/// Round up to the nearest integer.
ExactFloat ceil(in ExactFloat a) nothrow {
  return a.roundToPowerOf2(0, ExactFloat.RoundingMode.ROUND_TOWARD_POSITIVE);
}

/// Round down to the nearest integer.
ExactFloat floor(in ExactFloat a) nothrow {
  return a.roundToPowerOf2(0, ExactFloat.RoundingMode.ROUND_TOWARD_NEGATIVE);
}

/// Round to the nearest integer not larger in absolute value.
/// For example: f(-1.9) = -1, f(2.9) = 2.
ExactFloat trunc(in ExactFloat a) nothrow {
  return a.roundToPowerOf2(0, ExactFloat.RoundingMode.ROUND_TOWARD_ZERO);
}

/// Round to the nearest integer, rounding halfway cases away from zero.
/// For example: f(-0.5) = -1, f(0.5) = 1, f(1.5) = 2, f(2.5) = 3.
ExactFloat round(in ExactFloat a) nothrow {
  return a.roundToPowerOf2(0, ExactFloat.RoundingMode.ROUND_TIES_AWAY_FROM_ZERO);
}

/// Round to the nearest integer, rounding halfway cases to an even integer.
/// For example: f(-0.5) = 0, f(0.5) = 0, f(1.5) = 2, f(2.5) = 2.
ExactFloat rint(in ExactFloat a) nothrow {
  return a.roundToPowerOf2(0, ExactFloat.RoundingMode.ROUND_TIES_TO_EVEN);
}

/// A synonym for rint().
ExactFloat nearbyint(in ExactFloat a) nothrow { return rint(a); }

//////// Integer rounding functions that return C++ integer types.

/// Like rint(), but rounds to the nearest "long" value.  Returns the
/// minimum/maximum possible integer if the value is out of range.
long lrint(in ExactFloat a) nothrow {
  return a.toInteger!long(ExactFloat.RoundingMode.ROUND_TIES_TO_EVEN);
}

/// Like round(), but rounds to the nearest "long" value.  Returns the
/// minimum/maximum possible integer if the value is out of range.
long lround(in ExactFloat a) nothrow {
  return a.toInteger!long(ExactFloat.RoundingMode.ROUND_TIES_AWAY_FROM_ZERO);
}

//////// Remainder functions.

/// The remainder of dividing "a" by "b", where the quotient is rounded toward
/// zero to the nearest integer.  Similar to (a - trunc(a / b) * b).
ExactFloat fmod(in ExactFloat a, in ExactFloat b) nothrow {
  // Note that it is possible to implement this operation exactly, it just
  // hasn't been done.
  return ExactFloat.unimplemented();
}

/// The remainder of dividing "a" by "b", where the quotient is rounded to the
/// nearest integer, rounding halfway cases to an even integer.  Similar to
/// (a - rint(a / b) * b).
ExactFloat remainder(in ExactFloat a, in ExactFloat b) nothrow {
  // Note that it is possible to implement this operation exactly, it just
  // hasn't been done.
  return ExactFloat.unimplemented();
}

/// A synonym for remainder().
ExactFloat drem(in ExactFloat a, in ExactFloat b) nothrow {
  return remainder(a, b);
}

/**
 * Break the argument "a" into integer and fractional parts, each of which
 * has the same sign as "a".  The fractional part is returned, and the
 * integer part is stored in the output parameter "i_ptr".  Both output
 * values are set to have the same maximum precision as "a".
 */
ExactFloat modf(in ExactFloat a, out ExactFloat i_ptr) nothrow {
  // Note that it is possible to implement this operation exactly, it just
  // hasn't been done.
  return ExactFloat.unimplemented();
}

//////// Floating-point manipulation functions.

/**
 * Return an ExactFloat with the magnitude of "a" and the sign bit of "b".
 * (Note that an IEEE zero can be either positive or negative.)
 */
ExactFloat copysign(in ExactFloat a, in ExactFloat b) nothrow @nogc pure {
  return a.copyWithSign(b.sign);
}

/**
 * Convert "a" to a normalized fraction in the range [0.5, 1) times a power
 * of two.  Return the fraction and set "exp" to the exponent.  If "a" is
 * zero, infinity, or NaN then return "a" and set "exp" to zero.
 */
ExactFloat frexp(in ExactFloat a, out int exp) nothrow {
  if (!a.isNormal()) {
    // If a == 0, exp should be zero.  If a.is_inf() or a.is_nan(), exp is not
    // defined but the glibc implementation returns zero.
    exp = 0;
    return a;
  }
  exp = a.exp();
  return ldexp(a, -a.exp());
}

/// Return "a" multiplied by 2 raised to the power "exp".
ExactFloat ldexp(in ExactFloat a, int exp) nothrow {
  if (!a.isNormal()) return a;

  // To prevent integer overflow, we first clamp "exp" so that
  // (kMinExp - 1) <= (a_exp + exp) <= (kMaxExp + 1).
  int a_exp = a.exp();
  exp = algorithm.min(ExactFloat.MAX_EXP + 1 - a_exp,
      algorithm.max(ExactFloat.MIN_EXP - 1 + a_exp, exp));

  // Now modify the exponent and check for overflow/underflow.
  ExactFloat r = a;
  r._bnExp += exp;
  r.canonicalize();
  return r;
}

/// A synonym for ldexp().
ExactFloat scalbn(in ExactFloat a, int exp) nothrow {
  return ldexp(a, exp);
}

/// A version of ldexp() where "exp" is a long integer.
ExactFloat scalbln(in ExactFloat a, long exp) nothrow {
  // Clamp the exponent to the range of "int" in order to avoid truncation.
  exp = algorithm.max(cast(long) int.min, algorithm.min(cast(long) int.max, exp));
  return ldexp(a, cast(int) exp);
}

/**
 * Convert "a" to a normalized fraction in the range [1,2) times a power of
 * two, and return the exponent value as an integer.  This is equivalent to
 * lrint(floor(log2(fabs(a)))) but it is computed more efficiently.  Returns
 * the constants documented in the man page for zero, infinity, or NaN.
 */
int ilogb(in ExactFloat a) nothrow {
  if (a.isZero()) return math.FP_ILOGB0;
  if (a.isInf()) return int.max;
  if (a.isNan()) return math.FP_ILOGBNAN;
  // a.exp() assumes the significand is in the range [0.5, 1).
  return a.exp() - 1;
}

/**
 * Convert "a" to a normalized fraction in the range [1,2) times a power of
 * two, and return the exponent value as an ExactFloat.  This is equivalent to
 * floor(log2(fabs(a))) but it is computed more efficiently.
 */
ExactFloat logb(in ExactFloat a) nothrow {
  if (a.isZero()) return ExactFloat.infinity(-1);
  if (a.isInf()) return ExactFloat.infinity(+1);  // Even if a < 0.
  if (a.isNan()) return a;
  // exp() assumes the significand is in the range [0.5,1).
  return ExactFloat(a.exp() - 1);
}


