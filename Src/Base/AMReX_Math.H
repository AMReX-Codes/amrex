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

#elif defined (AMREX_USE_HIP)

using abs;
using ceil;
using copysign;
using floor;
using round;

#else

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;

#endif

}}

#endif
