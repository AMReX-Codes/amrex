#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>

#include <cmath>

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void halve (double& x)
{
    x *= 0.5;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double f (double x, double y)
{
    return std::cos(x*y);
}

struct S {
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double operator() (double a, double b, double c) const {
        return std::cos(a*b*c);
    }
};

struct P {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (double a, double b, double c) const {
        return std::cos(a*b*c);
    }
};

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        double pi = amrex::Math::pi<double>();
        double one = 1.0;

        auto g = [] AMREX_GPU_DEVICE (double a, double b, double c)
        {
            return std::sin(a*b*c);
        };

        S s{};
        P p{};

        amrex::Gpu::HostVector<double> ones(2);
        amrex::Gpu::HostVector<double> zeroes(3);
        auto* p1 = ones.data();
        auto* p0 = zeroes.data();

        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
        {
            p1[0] = callNoinline([] (double a) { return std::sin(a); }, pi*one*0.5);
            p1[1] = callNoinline(g, pi, one, 0.5);

            auto half = one;
#ifdef AMREX_USE_SYCL
            callNoinline([] (double& x) { halve(x); }, half);
            p0[0] = callNoinline([] (double a, double b) { return f(a,b); }, pi, half);
#else
            callNoinline(halve, half);
            p0[0] = callNoinline(f, pi, half);
#endif
            auto half2 = one;
            callNoinline([] (double& a) { a *= 0.5; }, half2);
            p0[1] = callNoinline(s, pi, one, half2);
            p0[2] = callNoinline(p, pi, one, half2);
        });
        Gpu::streamSynchronize();

        zeroes.push_back(callNoinline(f, pi, 0.5));
        zeroes.push_back(callNoinline(P{}, pi, one, 0.5));

        amrex::Print() << "ones: " << amrex::ToString(ones) << "\n"
                       << "zeroes: " << amrex::ToString(zeroes) << "\n";

        for (auto x : ones) {
            AMREX_ALWAYS_ASSERT(almostEqual(x, 1.0, 10));
        }

        for (auto x : zeroes) {
            AMREX_ALWAYS_ASSERT(std::abs(x) < 1.e-15);
        }
    }
    amrex::Finalize();
}
