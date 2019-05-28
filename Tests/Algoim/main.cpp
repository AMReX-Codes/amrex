#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_algoim.H>
#include <AMReX_algoim_K.H>

#include <algoim_quad.hpp>

using namespace amrex;
using namespace amrex::algoim;

namespace orig_algoim = Algoim;

static_assert(AMREX_SPACEDIM == 3, "3d only");

namespace {
    GpuArray<Real,3> normal (Real x, Real y, Real z) {
        Real norminv = 1./std::sqrt(x*x+y*y+z*z);
        return {x*norminv, y*norminv, z*norminv};
    }
}

namespace {
struct EBshape
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,3>& x) const
    {
        return (x(0) - cent[0])*(norm[0]) + (x(1) - cent[1])*(norm[1]) + (x(2) - cent[2])*(norm[2]);
    }

    template<typename T>
    blitz::TinyVector<T,3> grad(const blitz::TinyVector<T,3>& x) const
    {
        return blitz::TinyVector<double,3>(norm[0],norm[1],norm[2]);
    }

    EBshape (GpuArray<Real,3> const& c, GpuArray<Real,3> const& n)
        : cent(c), norm(c)
        {}

    EBshape (algoim::EBPlane const& rhs)
        : cent(rhs.cent), norm(rhs.norm)
        {}

    GpuArray<Real,3> cent{};
    GpuArray<Real,3> norm{};
};
}

void test_algoim (algoim::EBPlane& ebmax, Real& smax, algoim::EBPlane const& p);
Real test_algoim_perf (Vector<algoim::EBPlane> const& planes, Real& tnew, Real& told);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        long ntry = 1000;
        long nperf = 1000000;
        {
            ParmParse pp;
            pp.query("ntry", ntry);
            pp.query("nperf", nperf);
        }

        algoim::EBPlane ebmax;
        Real smax = -1.0;

        Vector<algoim::EBPlane> special_cases;
        for (int ibc = -1; ibc <= 1; ++ibc) {
        for (int jbc = -1; jbc <= 1; ++jbc) {
        for (int kbc = -1; kbc <= 1; ++kbc) {
        for (int ibn = -1; ibn <= 1; ++ibn) {
        for (int jbn = -1; jbn <= 1; ++jbn) {
        for (int kbn = -1; kbn <= 1; ++kbn) {
            if (ibn == 0 and jbn == 0 and kbn == 0) continue;
            special_cases.emplace_back(GpuArray<Real,3>{0.5*ibc,0.5*jbc,0.5*kbc},
                                       normal(ibn, jbn, kbn));
        }}}}}}

        for (auto const& p : special_cases)
        {
            test_algoim(ebmax, smax, p);
        }

        ResetRandomSeed(time(0));
        for (long itry = 0; itry < ntry; ++itry) {
            Real d1 = amrex::Random()-0.5;
            Real d2 = amrex::Random()-0.5;
            Real d3 = amrex::Random()-0.5;
            Real d4 = amrex::Random()-0.5;
            Real d5 = amrex::Random()-0.5;
            Real d6 = amrex::Random()-0.5;
            if (d4 != 0.0 or d5 != 0.0 or d6 != 0.0) {
                algoim::EBPlane p(GpuArray<Real,3>{d1,d2,d3},
                                  normal(d4,d5,d6));
                test_algoim(ebmax, smax, p);
            }
        }

        Vector<algoim::EBPlane> planes;
        for (long iperf = 0; iperf < nperf; ++iperf) {
            Real d1 = amrex::Random()-0.5;
            Real d2 = amrex::Random()-0.5;
            Real d3 = amrex::Random()-0.5;
            Real d4 = amrex::Random()-0.5;
            Real d5 = amrex::Random()-0.5;
            Real d6 = amrex::Random()-0.5;
            if (d4 != 0.0 or d5 != 0.0 or d6 != 0.0) {
                planes.emplace_back(GpuArray<Real,3>{d1,d2,d3},
                                    normal(d4,d5,d6));
            }
        }

        Real tnew, told;
        Real total = test_algoim_perf(planes, tnew, told);
        if (total != 0) amrex::Print() << "\n";
        amrex::Print() << "New and old performance times: " << tnew << " "
                       << told << "\n";

        amrex::Print().SetPrecision(17)
            << "\nMax diff. " << smax << " with\n    centroid("
            << ebmax.cent[0] << "," << ebmax.cent[1] << ","
            << ebmax.cent[2] << ")\n    normal  (" << ebmax.norm[0]
            << "," << ebmax.norm[1] << ","
            << ebmax.norm[2] << ")\n\n";
    }

    amrex::Finalize();
}

AMREX_FORCE_INLINE void
test_algoim_new (QuadratureRule const& q, Real& vol,
                 GpuArray<Real,algoim::numIntgs>& intg)
{
    vol             = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return 1.0;});
    intg[i_S_x    ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x; });
    intg[i_S_y    ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y; });
    intg[i_S_z    ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return z; });
    intg[i_S_x2   ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*x; });
    intg[i_S_y2   ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y*y; });
    intg[i_S_z2   ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return z*z; });
    intg[i_S_x_y  ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*y; });
    intg[i_S_x_z  ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*z; });
    intg[i_S_y_z  ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y*z; });
    intg[i_S_x2_y ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*x*y; });
    intg[i_S_x2_z ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*x*z; });
    intg[i_S_x_y2 ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*y*y; });
    intg[i_S_y2_z ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y*y*z; });
    intg[i_S_x_z2 ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*z*z; });
    intg[i_S_y_z2 ] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y*z*z; });
    intg[i_S_x2_y2] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*x*y*y;});
    intg[i_S_x2_z2] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return x*x*z*z;});
    intg[i_S_y2_z2] = q([] AMREX_GPU_HOST_DEVICE (Real x, Real y, Real z) {return y*y*z*z;});
}

AMREX_FORCE_INLINE void
test_algoim_old (orig_algoim::QuadratureRule<3> const& q, Real& vol,
                 GpuArray<Real,algoim::numIntgs>& intg)
{
    vol             = q([](const auto& w) {return 1.0;});
    intg[i_S_x    ] = q([](const auto& w) {return w[0]; });
    intg[i_S_y    ] = q([](const auto& w) {return w[1]; });
    intg[i_S_z    ] = q([](const auto& w) {return w[2]; });
    intg[i_S_x2   ] = q([](const auto& w) {return w[0]*w[0]; });
    intg[i_S_y2   ] = q([](const auto& w) {return w[1]*w[1]; });
    intg[i_S_z2   ] = q([](const auto& w) {return w[2]*w[2]; });
    intg[i_S_x_y  ] = q([](const auto& w) {return w[0]*w[1]; });
    intg[i_S_x_z  ] = q([](const auto& w) {return w[0]*w[2]; });
    intg[i_S_y_z  ] = q([](const auto& w) {return w[1]*w[2]; });
    intg[i_S_x2_y ] = q([](const auto& w) {return w[0]*w[0]*w[1]; });
    intg[i_S_x2_z ] = q([](const auto& w) {return w[0]*w[0]*w[2]; });
    intg[i_S_x_y2 ] = q([](const auto& w) {return w[0]*w[1]*w[1]; });
    intg[i_S_y2_z ] = q([](const auto& w) {return w[1]*w[1]*w[2]; });
    intg[i_S_x_z2 ] = q([](const auto& w) {return w[0]*w[2]*w[2]; });
    intg[i_S_y_z2 ] = q([](const auto& w) {return w[1]*w[2]*w[2]; });
    intg[i_S_x2_y2] = q([](const auto& w) {return w[0]*w[0]*w[1]*w[1];});
    intg[i_S_x2_z2] = q([](const auto& w) {return w[0]*w[0]*w[2]*w[2];});
    intg[i_S_y2_z2] = q([](const auto& w) {return w[1]*w[1]*w[2]*w[2];});
}

void
test_algoim (algoim::EBPlane& ebmax, Real& smax, algoim::EBPlane const& p)
{
    const QuadratureRule q = quadGen(p);
    Real vol;
    GpuArray<Real,algoim::numIntgs> intg;
    test_algoim_new(q, vol, intg);

    EBshape phi = p;
    const auto q2 = orig_algoim::quadGen<3>(phi, orig_algoim::BoundingBox<Real,3>({-0.5,-0.5,-0.5},{0.5,0.5,0.5}), -1, -1, 4);
    Real vol2;
    GpuArray<Real,algoim::numIntgs> intg2;
    test_algoim_old(q2, vol2, intg2);

    Real lsmax = std::abs(vol-vol2);
    for (int i = 0; i < numIntgs; ++i) {
        lsmax = std::max(lsmax, std::abs(intg[i]-intg2[i]));
    }

    if (lsmax > smax) {
        smax = lsmax;
        ebmax = p;
    }
}

Real test_algoim_perf (Vector<algoim::EBPlane> const& planes, Real& tnew, Real& told)
{
    Real total = 0.0;

    Real t0 = amrex::second();

    for (auto const& p : planes)
    {
        const QuadratureRule q = quadGen(p);
        Real vol;
        GpuArray<Real,algoim::numIntgs> intg;
        test_algoim_new(q, vol, intg);
        total += vol;
        for (int i = 0; i < intg.size(); ++i) {
            total += intg[i];
        }
    }

    Real t1 = amrex::second();

    for (auto const& p : planes)
    {
        EBshape phi = p;
        const auto q = orig_algoim::quadGen<3>(phi, orig_algoim::BoundingBox<Real,3>({-0.5,-0.5,-0.5},{0.5,0.5,0.5}), -1, -1, 4);
        Real vol;
        GpuArray<Real,algoim::numIntgs> intg;
        test_algoim_old(q, vol, intg);
        total += vol;
        for (int i = 0; i < intg.size(); ++i) {
            total += intg[i];
        }
    }

    Real t2 = amrex::second();

    tnew = t1-t0;
    told = t2-t1;
    return total;
}
