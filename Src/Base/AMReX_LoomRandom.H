#ifndef AMREX_LOOM_RANDOM_H_
#define AMREX_LOOM_RANDOM_H_
#include <AMReX_Config.H>

// Random numbers for AMReX_GPU_BACKEND=LOOM (AMReX_Random under LOOM).
//
// Device state: xorshift128+ seeded by splitmix64 from (seed, subsequence).
// Host generator: values drawn on the host with std::mt19937_64 and copied to
// device memory. Statistical quality suits AMReX's tests and particle
// initialization; the sequences are not rocRAND's or cuRAND's.

#ifdef AMREX_USE_LOOM

#include <AMReX_GpuQualifiers.H>  // AMREX_LOOM_DEVICE_PASS, AMREX_LOOM_DEVICE_INLINE

#include <cstddef>
#if !defined(AMREX_LOOM_DEVICE_PASS)
#include <cmath>
#endif

namespace amrex_loom {

struct RandState {
  unsigned long long s0;
  unsigned long long s1;
};

namespace rand_detail {
AMREX_LOOM_DEVICE_INLINE unsigned long long splitmix64(unsigned long long& x) {
  x += 0x9E3779B97F4A7C15ull;
  unsigned long long z = x;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
AMREX_LOOM_DEVICE_INLINE unsigned long long next(RandState* state) {
  unsigned long long s1 = state->s0;
  const unsigned long long s0 = state->s1;
  const unsigned long long result = s0 + s1;
  state->s0 = s0;
  s1 ^= s1 << 23;
  state->s1 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
  return result;
}
#if defined(AMREX_LOOM_DEVICE_PASS)
AMREX_LOOM_DEVICE_INLINE double log(double x) { return __builtin_log(x); }
AMREX_LOOM_DEVICE_INLINE double sqrt(double x) { return __builtin_sqrt(x); }
AMREX_LOOM_DEVICE_INLINE double cos(double x) { return __builtin_cos(x); }
AMREX_LOOM_DEVICE_INLINE double exp(double x) { return __builtin_exp(x); }
#else
inline double log(double x) { return std::log(x); }
inline double sqrt(double x) { return std::sqrt(x); }
inline double cos(double x) { return std::cos(x); }
inline double exp(double x) { return std::exp(x); }
#endif
}  // namespace rand_detail

AMREX_LOOM_DEVICE_INLINE void rand_init(unsigned long long seed, unsigned long long subsequence,
                                        unsigned long long offset, RandState* state) {
  unsigned long long x = seed ^ (subsequence * 0xD1B54A32D192ED03ull);
  RandState s;
  s.s0 = rand_detail::splitmix64(x);
  s.s1 = rand_detail::splitmix64(x);
  if (s.s0 == 0 && s.s1 == 0) s.s1 = 1;
  for (unsigned long long i = 0; i < offset; ++i) (void)rand_detail::next(&s);
  *state = s;
}

AMREX_LOOM_DEVICE_INLINE unsigned int rand_uint(RandState* state) {
  return static_cast<unsigned int>(rand_detail::next(state) >> 32);
}

// Uniform on (0, 1].
AMREX_LOOM_DEVICE_INLINE double rand_uniform_double(RandState* state) {
  unsigned long long bits = rand_detail::next(state) >> 11;  // 53 bits
  return (static_cast<double>(bits) + 1.0) * (1.0 / 9007199254740992.0);
}
AMREX_LOOM_DEVICE_INLINE float rand_uniform(RandState* state) {
  unsigned int bits = static_cast<unsigned int>(rand_detail::next(state) >> 40);  // 24 bits
  return (static_cast<float>(bits) + 1.0f) * (1.0f / 16777216.0f);
}

AMREX_LOOM_DEVICE_INLINE double rand_normal_double(RandState* state) {
  double u1 = rand_uniform_double(state);
  double u2 = rand_uniform_double(state);
  return rand_detail::sqrt(-2.0 * rand_detail::log(u1)) *
         rand_detail::cos(6.283185307179586476925 * u2);
}
AMREX_LOOM_DEVICE_INLINE float rand_normal(RandState* state) {
  return static_cast<float>(rand_normal_double(state));
}
AMREX_LOOM_DEVICE_INLINE unsigned int rand_poisson(RandState* state, double lambda) {
  double limit = rand_detail::exp(-lambda);
  double product = rand_uniform_double(state);
  unsigned int k = 0;
  while (product > limit) {
    ++k;
    product *= rand_uniform_double(state);
  }
  return k;
}

// --- host generator (AMReX_Loom.cpp) ------------------------------------------

struct RandGenerator;
RandGenerator* rand_generator_create (unsigned long long seed);
void rand_generator_destroy (RandGenerator* generator);
//! Fill n device values; these return once the values are in place.
void rand_generate_uniform (RandGenerator* generator, float* p, std::size_t n);
void rand_generate_uniform_double (RandGenerator* generator, double* p, std::size_t n);
void rand_generate_normal (RandGenerator* generator, float* p, std::size_t n, float mean,
                           float stddev);
void rand_generate_normal_double (RandGenerator* generator, double* p, std::size_t n,
                                  double mean, double stddev);

}  // namespace amrex_loom

#endif  // AMREX_USE_LOOM
#endif  // AMREX_LOOM_RANDOM_H_
