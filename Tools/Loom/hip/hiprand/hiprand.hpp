// Minimal hiprand replacement for AMReX_GPU_BACKEND=LOOM.
//
// Device state: xorshift128+ seeded by splitmix64 from (seed, subsequence).
// Host generator API: values drawn on the host and copied to device memory.
// Statistical quality suits AMReX's tests and particle initialization; it is
// not a drop-in for rocRAND's specific sequences.
#ifndef AMREX_LOOM_HIPRAND_HPP_
#define AMREX_LOOM_HIPRAND_HPP_

#include <hip/hip_runtime.h>

#include <cstddef>
#include <cstdint>
#if !defined(__loom__)
#include <cmath>
#include <random>
#include <vector>
#endif

struct hiprandState_t {
  unsigned long long s0;
  unsigned long long s1;
};

namespace amrex_loom_rand {
AMREX_LOOM_DEVICE_INLINE unsigned long long splitmix64(unsigned long long& x) {
  x += 0x9E3779B97F4A7C15ull;
  unsigned long long z = x;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
AMREX_LOOM_DEVICE_INLINE unsigned long long next(hiprandState_t* state) {
  unsigned long long s1 = state->s0;
  const unsigned long long s0 = state->s1;
  const unsigned long long result = s0 + s1;
  state->s0 = s0;
  s1 ^= s1 << 23;
  state->s1 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
  return result;
}
}  // namespace amrex_loom_rand

AMREX_LOOM_DEVICE_INLINE void hiprand_init(unsigned long long seed,
                                           unsigned long long subsequence,
                                           unsigned long long offset,
                                           hiprandState_t* state) {
  unsigned long long x = seed ^ (subsequence * 0xD1B54A32D192ED03ull);
  hiprandState_t s;
  s.s0 = amrex_loom_rand::splitmix64(x);
  s.s1 = amrex_loom_rand::splitmix64(x);
  if (s.s0 == 0 && s.s1 == 0) s.s1 = 1;
  for (unsigned long long i = 0; i < offset; ++i) (void)amrex_loom_rand::next(&s);
  *state = s;
}

AMREX_LOOM_DEVICE_INLINE unsigned int hiprand(hiprandState_t* state) {
  return static_cast<unsigned int>(amrex_loom_rand::next(state) >> 32);
}

// Uniform on (0, 1].
AMREX_LOOM_DEVICE_INLINE double hiprand_uniform_double(hiprandState_t* state) {
  unsigned long long bits = amrex_loom_rand::next(state) >> 11;  // 53 bits
  return (static_cast<double>(bits) + 1.0) * (1.0 / 9007199254740992.0);
}
AMREX_LOOM_DEVICE_INLINE float hiprand_uniform(hiprandState_t* state) {
  unsigned int bits = static_cast<unsigned int>(amrex_loom_rand::next(state) >> 40);  // 24 bits
  return (static_cast<float>(bits) + 1.0f) * (1.0f / 16777216.0f);
}

#if defined(__loom__)
#define AMREX_LOOM_LOG(x) __builtin_log(x)
#define AMREX_LOOM_SQRT(x) __builtin_sqrt(x)
#define AMREX_LOOM_COS(x) __builtin_cos(x)
#define AMREX_LOOM_EXP(x) __builtin_exp(x)
#else
#define AMREX_LOOM_LOG(x) std::log(x)
#define AMREX_LOOM_SQRT(x) std::sqrt(x)
#define AMREX_LOOM_COS(x) std::cos(x)
#define AMREX_LOOM_EXP(x) std::exp(x)
#endif

AMREX_LOOM_DEVICE_INLINE double hiprand_normal_double(hiprandState_t* state) {
  double u1 = hiprand_uniform_double(state);
  double u2 = hiprand_uniform_double(state);
  return AMREX_LOOM_SQRT(-2.0 * AMREX_LOOM_LOG(u1)) *
         AMREX_LOOM_COS(6.283185307179586476925 * u2);
}
AMREX_LOOM_DEVICE_INLINE float hiprand_normal(hiprandState_t* state) {
  return static_cast<float>(hiprand_normal_double(state));
}
AMREX_LOOM_DEVICE_INLINE unsigned int hiprand_poisson(hiprandState_t* state,
                                                      double lambda) {
  double limit = AMREX_LOOM_EXP(-lambda);
  double product = hiprand_uniform_double(state);
  unsigned int k = 0;
  while (product > limit) {
    ++k;
    product *= hiprand_uniform_double(state);
  }
  return k;
}

// ---- Host generator API ----------------------------------------------------
enum hiprandStatus_t { HIPRAND_STATUS_SUCCESS = 0, HIPRAND_STATUS_FAILURE = 1 };
enum hiprandRngType_t { HIPRAND_RNG_PSEUDO_DEFAULT = 400 };
struct amrex_loom_rand_generator;
typedef amrex_loom_rand_generator* hiprandGenerator_t;

#if !defined(__loom__)
struct amrex_loom_rand_generator {
  std::mt19937_64 engine{20260108ull};
};
inline hiprandStatus_t hiprandCreateGenerator(hiprandGenerator_t* g, hiprandRngType_t) {
  *g = new amrex_loom_rand_generator;
  return HIPRAND_STATUS_SUCCESS;
}
inline hiprandStatus_t hiprandDestroyGenerator(hiprandGenerator_t g) {
  delete g;
  return HIPRAND_STATUS_SUCCESS;
}
inline hiprandStatus_t hiprandSetPseudoRandomGeneratorSeed(hiprandGenerator_t g,
                                                           unsigned long long seed) {
  g->engine.seed(seed);
  return HIPRAND_STATUS_SUCCESS;
}
inline hiprandStatus_t hiprandSetStream(hiprandGenerator_t, hipStream_t) {
  return HIPRAND_STATUS_SUCCESS;
}
template <class T, class D>
inline hiprandStatus_t amrex_loom_generate(hiprandGenerator_t g, T* p, std::size_t n, D dist) {
  std::vector<T> host(n);
  for (auto& v : host) v = static_cast<T>(dist(g->engine));
  hipMemcpy(p, host.data(), n * sizeof(T), hipMemcpyHostToDevice);
  return HIPRAND_STATUS_SUCCESS;
}
inline hiprandStatus_t hiprandGenerateUniform(hiprandGenerator_t g, float* p, std::size_t n) {
  return amrex_loom_generate(g, p, n, std::uniform_real_distribution<float>(0.0f, 1.0f));
}
inline hiprandStatus_t hiprandGenerateUniformDouble(hiprandGenerator_t g, double* p,
                                                    std::size_t n) {
  return amrex_loom_generate(g, p, n, std::uniform_real_distribution<double>(0.0, 1.0));
}
inline hiprandStatus_t hiprandGenerateNormal(hiprandGenerator_t g, float* p, std::size_t n,
                                             float mean, float stddev) {
  return amrex_loom_generate(g, p, n, std::normal_distribution<float>(mean, stddev));
}
inline hiprandStatus_t hiprandGenerateNormalDouble(hiprandGenerator_t g, double* p,
                                                   std::size_t n, double mean,
                                                   double stddev) {
  return amrex_loom_generate(g, p, n, std::normal_distribution<double>(mean, stddev));
}
#else
hiprandStatus_t hiprandCreateGenerator(hiprandGenerator_t*, hiprandRngType_t);
hiprandStatus_t hiprandDestroyGenerator(hiprandGenerator_t);
hiprandStatus_t hiprandSetPseudoRandomGeneratorSeed(hiprandGenerator_t, unsigned long long);
hiprandStatus_t hiprandSetStream(hiprandGenerator_t, hipStream_t);
hiprandStatus_t hiprandGenerateUniform(hiprandGenerator_t, float*, std::size_t);
hiprandStatus_t hiprandGenerateUniformDouble(hiprandGenerator_t, double*, std::size_t);
hiprandStatus_t hiprandGenerateNormal(hiprandGenerator_t, float*, std::size_t, float, float);
hiprandStatus_t hiprandGenerateNormalDouble(hiprandGenerator_t, double*, std::size_t, double,
                                            double);
#endif

#endif  // AMREX_LOOM_HIPRAND_HPP_
