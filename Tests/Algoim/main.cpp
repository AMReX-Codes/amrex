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
    Array<Real,3> normal (Real x, Real y, Real z) {
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

    EBshape (Array<Real,3> const& c, Array<Real,3> const& n)
        : cent(c), norm(c)
        {}

    EBshape (algoim::EBPlane const& rhs)
        : cent(rhs.cent), norm(rhs.norm)
        {}

    Array<Real,3> cent{};
    Array<Real,3> norm{};
};
}

void test_algoim (algoim::EBPlane& ebmax, Real& smax, algoim::EBPlane const& p);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        long ntry = 1000;
        {
            ParmParse pp;
            pp.query("ntry", ntry);
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
            special_cases.emplace_back(Array<Real,3>{0.5*ibc,0.5*jbc,0.5*kbc},
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
                algoim::EBPlane p(Array<Real,3>{d1,d2,d3},
                                  normal(d4,d5,d6));
                test_algoim(ebmax, smax, p);
            }
        }

        amrex::Print() << "\nMax diff. " << smax << " with\n    centroid("
                       << ebmax.cent[0] << "," << ebmax.cent[1] << ","
                       << ebmax.cent[2] << ")\n    normal  (" << ebmax.norm[0]
                       << "," << ebmax.norm[1] << ","
                       << ebmax.norm[2] << ")\n\n";
    }

    amrex::Finalize();
}

void
test_algoim (algoim::EBPlane& ebmax, Real& smax, algoim::EBPlane const& p)
{
    const QuadratureRule q = quadGen(p);
    Array<Real,algoim::numIntgs> intg;
    Real vol        = q([](Real x, Real y, Real z) {return 1.0;});
    intg[i_S_x    ] = q([](Real x, Real y, Real z) {return x; });
    intg[i_S_y    ] = q([](Real x, Real y, Real z) {return y; });
    intg[i_S_z    ] = q([](Real x, Real y, Real z) {return z; });
    intg[i_S_x2   ] = q([](Real x, Real y, Real z) {return x*x; });
    intg[i_S_y2   ] = q([](Real x, Real y, Real z) {return y*y; });
    intg[i_S_z2   ] = q([](Real x, Real y, Real z) {return z*z; });
    intg[i_S_x_y  ] = q([](Real x, Real y, Real z) {return x*y; });
    intg[i_S_x_z  ] = q([](Real x, Real y, Real z) {return x*z; });
    intg[i_S_y_z  ] = q([](Real x, Real y, Real z) {return y*z; });
    intg[i_S_x2_y ] = q([](Real x, Real y, Real z) {return x*x*y; });
    intg[i_S_x2_z ] = q([](Real x, Real y, Real z) {return x*x*z; });
    intg[i_S_x_y2 ] = q([](Real x, Real y, Real z) {return x*y*y; });
    intg[i_S_y2_z ] = q([](Real x, Real y, Real z) {return y*y*z; });
    intg[i_S_x_z2 ] = q([](Real x, Real y, Real z) {return x*z*z; });
    intg[i_S_y_z2 ] = q([](Real x, Real y, Real z) {return y*z*z; });
    intg[i_S_x2_y2] = q([](Real x, Real y, Real z) {return x*x*y*y;});
    intg[i_S_x2_z2] = q([](Real x, Real y, Real z) {return x*x*z*z;});
    intg[i_S_y2_z2] = q([](Real x, Real y, Real z) {return y*y*z*z;});

    EBshape phi = p;
    const auto q2 = orig_algoim::quadGen<3>(phi, orig_algoim::BoundingBox<Real,3>({-0.5,-0.5,-0.5},{0.5,0.5,0.5}), -1, -1, 4);
    Array<Real,algoim::numIntgs> intg2;
    Real vol2        = q2([](const auto& w) {return 1.0;});
    intg2[i_S_x    ] = q2([](const auto& w) {return w[0]; });
    intg2[i_S_y    ] = q2([](const auto& w) {return w[1]; });
    intg2[i_S_z    ] = q2([](const auto& w) {return w[2]; });
    intg2[i_S_x2   ] = q2([](const auto& w) {return w[0]*w[0]; });
    intg2[i_S_y2   ] = q2([](const auto& w) {return w[1]*w[1]; });
    intg2[i_S_z2   ] = q2([](const auto& w) {return w[2]*w[2]; });
    intg2[i_S_x_y  ] = q2([](const auto& w) {return w[0]*w[1]; });
    intg2[i_S_x_z  ] = q2([](const auto& w) {return w[0]*w[2]; });
    intg2[i_S_y_z  ] = q2([](const auto& w) {return w[1]*w[2]; });
    intg2[i_S_x2_y ] = q2([](const auto& w) {return w[0]*w[0]*w[1]; });
    intg2[i_S_x2_z ] = q2([](const auto& w) {return w[0]*w[0]*w[2]; });
    intg2[i_S_x_y2 ] = q2([](const auto& w) {return w[0]*w[1]*w[1]; });
    intg2[i_S_y2_z ] = q2([](const auto& w) {return w[1]*w[1]*w[2]; });
    intg2[i_S_x_z2 ] = q2([](const auto& w) {return w[0]*w[2]*w[2]; });
    intg2[i_S_y_z2 ] = q2([](const auto& w) {return w[1]*w[2]*w[2]; });
    intg2[i_S_x2_y2] = q2([](const auto& w) {return w[0]*w[0]*w[1]*w[1];});
    intg2[i_S_x2_z2] = q2([](const auto& w) {return w[0]*w[0]*w[2]*w[2];});
    intg2[i_S_y2_z2] = q2([](const auto& w) {return w[1]*w[1]*w[2]*w[2];});

    Real lsmax = std::abs(vol-vol2);
    for (int i = 0; i < numIntgs; ++i) {
        lsmax = std::max(lsmax, std::abs(intg[i]-intg2[i]));
    }

    if (lsmax > smax) {
        smax = lsmax;
        ebmax = p;
    }
}
