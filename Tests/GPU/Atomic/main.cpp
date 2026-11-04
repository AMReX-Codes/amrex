#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

// Tests for Gpu::Atomic and HostDevice::Atomic. Each kernel has many
// threads update a few addresses, and the result is checked against the
// same arithmetic done serially on the host. amrex::For is used instead
// of ParallelFor because the loop bodies are not SIMD safe on the host.

using namespace amrex;

namespace {

constexpr int N = 1 << 16;
constexpr int nbins = 7;

// Mixed-sign values with magnitudes beyond 32 bits for the 64-bit types.
// Floating-point values are multiples of 0.25 so that sums are exact.
template <typename T>
AMREX_GPU_HOST_DEVICE T value (int i)
{
    if constexpr (std::is_floating_point_v<T>) {
        return T((i % 2) ? -0.5 : 0.25) * T(i % 7 + 1);
    } else {
        constexpr int shift = (sizeof(T) > 4) ? 40 : 20;
        constexpr T big = T(1) << shift;
        if constexpr (std::is_signed_v<T>) {
            return (i % 2 == 0) ? (big + T(i)) : -(big + T(i));
        } else {
            return big + T(i);
        }
    }
}

template <typename T>
std::vector<T> to_host (Gpu::DeviceVector<T> const& d)
{
    std::vector<T> h(d.size());
    Gpu::copyAsync(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    Gpu::streamSynchronize();
    return h;
}

template <typename T>
void test_add ()
{
    Gpu::DeviceVector<T> d(nbins, T(0));
    T* p = d.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        T const v = value<T>(i);
        Gpu::Atomic::Add(p + i%nbins, v);
        Gpu::Atomic::AddNoRet(p + i%nbins, v);
        HostDevice::Atomic::Add(p + i%nbins, v);
    });
    auto h = to_host(d);
    std::vector<T> expected(nbins, T(0));
    for (int i = 0; i < N; ++i) { expected[i%nbins] += T(3) * value<T>(i); }
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(h[b] == expected[b]);
    }
}

// Add returns the old value: incrementing one counter from N threads must
// hand out every value in [0, N) exactly once.
template <typename T>
void test_fetch_add ()
{
    Gpu::DeviceVector<T> d(2, T(0));
    Gpu::DeviceVector<T> olds(2*N);
    T* p = d.data();
    T* po = olds.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        po[i]   = Gpu::Atomic::Add(p, T(1));
        po[N+i] = HostDevice::Atomic::FetchAdd(p+1, T(1));
    });
    auto h = to_host(d);
    auto ho = to_host(olds);
    AMREX_ALWAYS_ASSERT(h[0] == T(N) && h[1] == T(N));
    for (int k = 0; k < 2; ++k) {
        std::sort(ho.begin() + k*N, ho.begin() + (k+1)*N);
        for (int i = 0; i < N; ++i) {
            AMREX_ALWAYS_ASSERT(ho[k*N+i] == T(i));
        }
    }
}

template <typename T>
void test_min_max ()
{
    std::vector<T> hv(N);
    for (int i = 0; i < N; ++i) { hv[i] = value<T>(i); }
    if constexpr (std::is_floating_point_v<T>) {
        hv[3] = T(-0.0);
        hv[5] = T(0.0);
        hv[9] = std::numeric_limits<T>::denorm_min();
        hv[11] = -std::numeric_limits<T>::denorm_min();
    }
    Gpu::DeviceVector<T> dv(N);
    Gpu::copyAsync(Gpu::hostToDevice, hv.begin(), hv.end(), dv.begin());
    T const* pv = dv.data();

    // Start from the extreme values, and from zero so that memory often
    // holds the opposite sign of the incoming value.
    for (int pass = 0; pass < 2; ++pass) {
        T const init_min = (pass == 0) ? std::numeric_limits<T>::max() : T(0);
        T const init_max = (pass == 0) ? std::numeric_limits<T>::lowest() : T(0);
        Gpu::DeviceVector<T> dmin(nbins, init_min);
        Gpu::DeviceVector<T> dmax(nbins, init_max);
        T* pmin = dmin.data();
        T* pmax = dmax.data();
        amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            Gpu::Atomic::Min(pmin + i%nbins, pv[i]);
            Gpu::Atomic::Max(pmax + i%nbins, pv[i]);
        });
        auto hmin = to_host(dmin);
        auto hmax = to_host(dmax);
        std::vector<T> emin(nbins, init_min), emax(nbins, init_max);
        for (int i = 0; i < N; ++i) {
            emin[i%nbins] = std::min(emin[i%nbins], hv[i]);
            emax[i%nbins] = std::max(emax[i%nbins], hv[i]);
        }
        for (int b = 0; b < nbins; ++b) {
            AMREX_ALWAYS_ASSERT(hmin[b] == emin[b]);   // == treats -0.0 and 0.0 as equal
            AMREX_ALWAYS_ASSERT(hmax[b] == emax[b]);
        }
    }

    // Signed zeros must not be mistaken for the most negative or positive
    // value. Only equality to zero is checked, so either sign may come back.
    if constexpr (std::is_floating_point_v<T>) {
        std::vector<T> hz = {T(-0.0), T(-5.0), T(-0.0), T(-3.0), T(0.0), T(-1.0)};
        Gpu::DeviceVector<T> dz(hz.size());
        Gpu::copyAsync(Gpu::hostToDevice, hz.begin(), hz.end(), dz.begin());
        Gpu::DeviceVector<T> dr(2, T(-0.0));
        T const* pz = dz.data();
        T* pr = dr.data();
        amrex::For(int(hz.size()), [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            Gpu::Atomic::Min(pr, pz[i]);
            Gpu::Atomic::Max(pr+1, pz[i]);
        });
        auto hr = to_host(dr);
        AMREX_ALWAYS_ASSERT(hr[0] == T(-5.0));
        AMREX_ALWAYS_ASSERT(hr[1] == T(0.0));
    }
}

void test_logical ()
{
    Gpu::DeviceVector<int> dor(nbins, 0);
    Gpu::DeviceVector<int> dand(nbins, 1);
    int* por = dor.data();
    int* pand = dand.data();
    constexpr int special = 12345;
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        Gpu::Atomic::LogicalOr(por + i%nbins, (i == special) ? 7 : 0);  // any nonzero is true
        Gpu::Atomic::LogicalAnd(pand + i%nbins, (i == special) ? 0 : 3);
    });
    auto hor = to_host(dor);
    auto hand = to_host(dand);
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT((hor[b] != 0) == (b == special % nbins));
        AMREX_ALWAYS_ASSERT((hand[b] != 0) == (b != special % nbins));
    }
}

// The old values handed out by Exch plus the final contents must be exactly
// the initial contents plus everything that was written.
template <typename T>
void test_exch ()
{
    Gpu::DeviceVector<T> d(nbins, T(0));
    Gpu::DeviceVector<T> olds(N);
    T* p = d.data();
    T* po = olds.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        po[i] = Gpu::Atomic::Exch(p + i%nbins, T(i+1));
    });
    auto h = to_host(d);
    auto ho = to_host(olds);
    ho.insert(ho.end(), h.begin(), h.end());
    std::sort(ho.begin(), ho.end());
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(ho[b] == T(0));
    }
    for (int i = 0; i < N; ++i) {
        AMREX_ALWAYS_ASSERT(ho[nbins+i] == T(i+1));
    }
}

template <typename T>
void test_cas ()
{
    Gpu::DeviceVector<T> d(nbins, T(0));
    T* p = d.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        T* a = p + i%nbins;
        T old = *a, assumed;
        do {
            assumed = old;
            old = Gpu::Atomic::CAS(a, assumed, assumed + T(1));
        } while (assumed != old);
    });
    auto h = to_host(d);
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(h[b] == T(N/nbins + (b < N%nbins ? 1 : 0)));
    }
}

// Every thread multiplies (divides) by 1, except for a few that use 2, so
// that the results stay exact for all types.
template <typename T>
void test_mul_div ()
{
    constexpr int stride = 1000;
    Gpu::DeviceVector<T> dmul(nbins, T(1));
    Gpu::DeviceVector<T> ddiv(nbins, T(1 << 20));
    T* pmul = dmul.data();
    T* pdiv = ddiv.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        T const f = (i % stride == 0) ? T(2) : T(1);
        Gpu::Atomic::Multiply(pmul + i%nbins, f);
        Gpu::Atomic::Divide(pdiv + i%nbins, f);
    });
    auto hmul = to_host(dmul);
    auto hdiv = to_host(ddiv);
    std::vector<T> emul(nbins, T(1)), ediv(nbins, T(1 << 20));
    for (int i = 0; i < N; i += stride) {
        emul[i%nbins] *= T(2);
        ediv[i%nbins] /= T(2);
    }
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(hmul[b] == emul[b]);
        AMREX_ALWAYS_ASSERT(hdiv[b] == ediv[b]);
    }
}

// A value with random bits. For floating-point types NaNs are skipped, since
// Min/Max are unspecified for them, and so are denormals, so that the host
// comparison does not depend on flush-to-zero settings. Zeros of both signs
// and infinities are included.
template <typename T>
T random_value (std::mt19937_64& rng)
{
    using U = std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>;
    constexpr int ebits = (sizeof(T) == 4) ? 8 : 11;
    constexpr int mbits = (sizeof(T) == 4) ? 23 : 52;
    for (;;) {
        U const u = U(rng());
        if constexpr (std::is_floating_point_v<T>) {
            U const e = (u >> mbits) & ((U(1) << ebits) - 1);
            U const m = u & ((U(1) << mbits) - 1);
            if (m != 0 && (e == 0 || e == ((U(1) << ebits) - 1))) { continue; }
            T x;
            std::memcpy(&x, &u, sizeof(T));
            return x;
        } else {
            return T(u);
        }
    }
}

template <typename T>
void test_min_max_random (std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::vector<T> hv(N), emin(nbins), emax(nbins);
    for (auto& x : hv) { x = random_value<T>(rng); }
    for (int b = 0; b < nbins; ++b) {
        emin[b] = random_value<T>(rng);   // random starting contents too
        emax[b] = random_value<T>(rng);
    }
    Gpu::DeviceVector<T> dv(N), dmin(nbins), dmax(nbins);
    Gpu::copyAsync(Gpu::hostToDevice, hv.begin(), hv.end(), dv.begin());
    Gpu::copyAsync(Gpu::hostToDevice, emin.begin(), emin.end(), dmin.begin());
    Gpu::copyAsync(Gpu::hostToDevice, emax.begin(), emax.end(), dmax.begin());
    T const* pv = dv.data();
    T* pmin = dmin.data();
    T* pmax = dmax.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        Gpu::Atomic::Min(pmin + i%nbins, pv[i]);
        Gpu::Atomic::Max(pmax + i%nbins, pv[i]);
    });
    auto hmin = to_host(dmin);
    auto hmax = to_host(dmax);
    for (int i = 0; i < N; ++i) {
        emin[i%nbins] = std::min(emin[i%nbins], hv[i]);
        emax[i%nbins] = std::max(emax[i%nbins], hv[i]);
    }
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(hmin[b] == emin[b]);
        AMREX_ALWAYS_ASSERT(hmax[b] == emax[b]);
    }
}

// Random integers, mixed sign, small enough that the sums cannot overflow.
template <typename T>
void test_add_random (std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::vector<T> hv(N), expected(nbins, T(0));
    for (int i = 0; i < N; ++i) {
        std::uint64_t const u = rng();
        auto v = T(u >> ((sizeof(T) == 4) ? 50 : 24));
        if constexpr (std::is_signed_v<T>) { if (u & 1) { v = -v; } }
        hv[i] = v;
        expected[i%nbins] += v;
    }
    Gpu::DeviceVector<T> dv(N), d(nbins, T(0));
    Gpu::copyAsync(Gpu::hostToDevice, hv.begin(), hv.end(), dv.begin());
    T const* pv = dv.data();
    T* p = d.data();
    amrex::For(N, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        Gpu::Atomic::AddNoRet(p + i%nbins, pv[i]);
    });
    auto h = to_host(d);
    for (int b = 0; b < nbins; ++b) {
        AMREX_ALWAYS_ASSERT(h[b] == expected[b]);
    }
}

template <typename T>
struct IsNonnegative {
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator () (T const val) const noexcept {
        return val >= T(0.);
    }
};

template <typename T>
void test_if ()
{
    constexpr int M = 8000;
    Gpu::DeviceVector<T> v_d(M);
    auto v_d_ptr = v_d.data();
    amrex::For(M, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        v_d_ptr[i] = T(i);
    });

    for (int j = 0; j < 1000; ++j) {
        amrex::For(M, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            // add -1.0 to v_d_ptr[i] if the result would be non-negative
            Gpu::Atomic::If(&v_d_ptr[i], T(-1.0), amrex::Plus<T>(), IsNonnegative<T>());
            // can also use +1.0 amrex::Minus<T>
            Gpu::Atomic::If(&v_d_ptr[M-i-1], T(1.0), amrex::Minus<T>(), IsNonnegative<T>());
        });

        amrex::For(M, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            // lambdas also work
            Gpu::Atomic::If(&v_d_ptr[M-i-1], T(1.0), amrex::Minus<T>(),
                            [=] AMREX_GPU_DEVICE (T tmp) noexcept
                            {
                                return (tmp >= T(0.0));
                            });
            Gpu::Atomic::If(&v_d_ptr[i], T(1.0), amrex::Minus<T>(),
                            [=] AMREX_GPU_DEVICE (T tmp) noexcept
                            {
                                return (tmp >= T(0.0));
                            });
        });
    }
    auto v_h = to_host(v_d);

    // The first half should all be 0, the second half should go from 0 to M/2-1
    for (int i = 0; i < M/2; ++i) {
        AMREX_ALWAYS_ASSERT(v_h[i] == 0);
    }
    for (int i = M/2, j = 0; i < M; ++i, ++j) {
        AMREX_ALWAYS_ASSERT(v_h[i] == j);
    }
}

void test_host ()
{
    int sum = 0;
    Long lsum = 0;
    double dsum = 0.;
    std::vector<int> olds(N);
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; ++i) {
        HostDevice::Atomic::Add(&sum, 1);
        HostDevice::Atomic::Add(&lsum, Long(i));
        HostDevice::Atomic::Add(&dsum, 0.25);
        olds[i] = HostDevice::Atomic::FetchAdd(&sum, 0);  // snapshot
    }
    AMREX_ALWAYS_ASSERT(sum == N);
    AMREX_ALWAYS_ASSERT(lsum == Long(N) * Long(N-1) / 2);
    AMREX_ALWAYS_ASSERT(dsum == 0.25 * N);
    for (int i = 0; i < N; ++i) {
        AMREX_ALWAYS_ASSERT(olds[i] >= 1 && olds[i] <= N);
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        test_add<int>();
        test_add<unsigned int>();
        test_add<Long>();
        test_add<unsigned long long>();
        test_add<float>();
        test_add<double>();

        test_fetch_add<int>();
        test_fetch_add<unsigned int>();
        test_fetch_add<Long>();
        test_fetch_add<unsigned long long>();
        test_fetch_add<float>();
        test_fetch_add<double>();

        test_min_max<int>();
        test_min_max<unsigned int>();
        test_min_max<Long>();
        test_min_max<unsigned long long>();
        test_min_max<float>();
        test_min_max<double>();

        // Random seeds, printed so that a failure can be reproduced with
        // the runtime parameter "seed", which then runs a single round.
        std::random_device rd;
        long long seed_in = -1;
        ParmParse pp;
        pp.query("seed", seed_in);
        int const nrounds = (seed_in >= 0) ? 1 : 8;
        for (int round = 0; round < nrounds; ++round) {
            auto const seed = (seed_in >= 0) ? static_cast<std::uint64_t>(seed_in)
                                             : (((std::uint64_t(rd()) << 32) | rd()) >> 1);
            amrex::AllPrint() << "Proc " << ParallelDescriptor::MyProc()
                              << ": seed = " << seed << "\n";
            test_min_max_random<int>(seed);
            test_min_max_random<unsigned int>(seed);
            test_min_max_random<Long>(seed);
            test_min_max_random<unsigned long long>(seed);
            test_min_max_random<float>(seed);
            test_min_max_random<double>(seed);

            test_add_random<int>(seed);
            test_add_random<unsigned int>(seed);
            test_add_random<Long>(seed);
            test_add_random<unsigned long long>(seed);
        }

        test_logical();

        test_exch<int>();
        test_exch<unsigned int>();
        test_exch<unsigned long long>();
        test_exch<float>();

        test_cas<int>();
        test_cas<unsigned int>();
        test_cas<unsigned long long>();

        test_mul_div<int>();
        test_mul_div<Long>();
        test_mul_div<float>();
        test_mul_div<double>();

        test_if<int>();
        test_if<Long>();
        test_if<float>();
        test_if<double>();

        test_host();

        amrex::Print() << "All atomic tests passed.\n";
    }
    amrex::Finalize();
}
