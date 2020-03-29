#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_

#include <cmath>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;
#endif

namespace amrex { namespace Math {

#ifdef AMREX_USE_DPCPP

//using sycl::abs;
// We have to do this because somehow sycl::abs(int) return unsigned int.
template <typename T> T abs (T a) { return sycl::abs(a); }

using sycl::ceil;
using sycl::copysign;
using sycl::floor;
using sycl::round;

using sycl::exp;
using sycl::expm1;
using sycl::log;
using sycl::log10;
using sycl::log2;
using sycl::log1p;
using sycl::pow;
using sycl::sqrt;
using sycl::cbrt;
using sycl::hypot;
using sycl::sin;
using sycl::cos;
using sycl::tan;
using sycl::asin;
using sycl::acos;
using sycl::atan;
using sycl::atan2;
using sycl::sinh;
using sycl::cosh;
using sycl::tanh;
using sycl::asinh;
using sycl::acosh;
using sycl::atanh;

#elif defined (AMREX_USE_HIP)

using abs;
using ceil;
using copysign;
using floor;
using round;

using std::exp;
using std::expm1;
using std::log;
using std::log10;
using std::log2;
using std::log1p;
using std::pow;
using std::sqrt;
using std::cbrt;
using std::hypot;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;
using std::sinh;
using std::cosh;
using std::tanh;
using std::asinh;
using std::acosh;
using std::atanh;

#else

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;

using std::exp;
using std::expm1;
using std::log;
using std::log10;
using std::log2;
using std::log1p;
using std::pow;
using std::sqrt;
using std::cbrt;
using std::hypot;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;
using std::sinh;
using std::cosh;
using std::tanh;
using std::asinh;
using std::acosh;
using std::atanh;

#endif

}}

#endif
