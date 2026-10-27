/**
 * \file main.cpp
 * \brief Sign test for amrex::FillSignedDistance.
 *
 * The sign of every node of the signed distance is compared with the EB2
 * level set (negative = fluid, zero or positive = body) and with the cell
 * flags around it (all covered = body, all regular = fluid).  Cases: walls a roundoff
 * distance inside the domain faces, as MFIX-Exa builds its box, with and
 * without a refined level-set grid; a sphere with fluid inside; a sphere with
 * fluid outside in a periodic domain; a sphere centered on a node with a
 * radius through nodes, where the sign is also checked against the exact
 * distance.
 */

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_utils.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <functional>
#include <limits>

using namespace amrex;

namespace {

// exact, if given, is the exact signed distance (positive in the fluid).
int check (std::string const& name, Geometry const& geom, Geometry const& ls_geom,
           EB2::Level const& eb_level, EB2::Level const& ls_level, int refratio,
           std::function<Real(RealArray const&)> const& exact)
{
    const int ng_eb = 4;
    const int ng_ls = 2*refratio;

    BoxArray ba(geom.Domain());
    ba.maxSize(16);
    DistributionMapping dm(ba);
    EBFArrayBoxFactory fac(eb_level, geom, ba, dm, {ng_eb, ng_eb, ng_eb}, EBSupport::full);

    BoxArray lsba = amrex::convert(ba, IntVect(1));
    lsba.refine(refratio);
    MultiFab sd(lsba, dm, 1, ng_ls);
    FillSignedDistance(sd, ls_level, fac, refratio, true);

    MultiFab lsn(lsba, dm, 1, ng_ls);
    ls_level.fillLevelSet(lsn, ls_geom);

    auto const& flags = fac.getMultiEBCellFlagFab();

    int nerr = 0;
    for (MFIter mfi(sd); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.growntilebox();
        FArrayBox hs(sd[mfi].box(), 1, The_Pinned_Arena());
        FArrayBox hl(lsn[mfi].box(), 1, The_Pinned_Arena());
        BaseFab<EBCellFlag> hf(flags[mfi].box(), 1, The_Pinned_Arena());
        Gpu::dtoh_memcpy_async(hs.dataPtr(), sd[mfi].dataPtr(), hs.nBytes());
        Gpu::dtoh_memcpy_async(hl.dataPtr(), lsn[mfi].dataPtr(), hl.nBytes());
        Gpu::dtoh_memcpy_async(hf.dataPtr(), flags[mfi].dataPtr(), hf.nBytes());
        Gpu::streamSynchronize();
        auto const& s = hs.const_array();
        auto const& l = hl.const_array();
        auto const& f = hf.const_array();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            // EB2 level set: negative is fluid, zero or positive is body.  A zero
            // distance is accepted either way: the EB of a cell with two walls on
            // its faces is one facet through the cell center.
            bool bad = (l(i,j,k) < 0.0_rt) ? (s(i,j,k) < 0.0_rt) : (s(i,j,k) > 0.0_rt);
            // cells around the node on the EB grid
            IntVect iv(AMREX_D_DECL(i,j,k));
            IntVect ic = amrex::coarsen(iv, refratio);
            bool all_covered = true, all_regular = true;
            Box const nb(ic - 1, ic);
            amrex::LoopOnCpu(nb, [&] (int ii, int jj, int kk)
            {
                all_covered = all_covered && f(ii,jj,kk).isCovered();
                all_regular = all_regular && f(ii,jj,kk).isRegular();
            });
            if (all_covered && s(i,j,k) > 0.0_rt) { bad = true; }
            if (all_regular && s(i,j,k) < 0.0_rt) { bad = true; }
            // Nodes clearly off the surface must have the exact sign, not zero.
            if (exact) {
                RealArray p;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    p[d] = ls_geom.ProbLo(d) + Real(iv[d])*ls_geom.CellSize(d);
                }
                Real const e = exact(p);
                Real const tol = 0.05_rt*ls_geom.CellSize(0);
                if ((e > tol && s(i,j,k) <= 0.0_rt) || (e < -tol && s(i,j,k) >= 0.0_rt)) {
                    bad = true;
                }
            }
            if (bad) {
                if (nerr < 10) {
                    amrex::AllPrint() << name << ": node " << iv << " level set " << l(i,j,k)
                                      << " signed distance " << s(i,j,k)
                                      << " all_covered " << all_covered
                                      << " all_regular " << all_regular << "\n";
                }
                ++nerr;
            }
        });
    }
    ParallelDescriptor::ReduceIntSum(nerr);
    amrex::Print() << name << ": " << nerr << " nodes with the wrong sign\n";
    return nerr;
}

template <class F>
int run_case (std::string const& name, F const& f, Geometry const& geom, int refratio,
              std::function<Real(RealArray const&)> const& exact = {})
{
    auto gshop = EB2::makeShop(f);
    EB2::Build(gshop, geom, 0, 100);
    EB2::Level const& eb_level = EB2::IndexSpace::top().getLevel(geom);
    if (refratio == 1) {
        return check(name, geom, geom, eb_level, eb_level, 1, exact);
    }
    Geometry ls_geom = amrex::refine(geom, refratio);
    EB2::Build(gshop, ls_geom, 0, 100);
    EB2::Level const& ls_level = EB2::IndexSpace::top().getLevel(ls_geom);
    return check(name, geom, ls_geom, eb_level, ls_level, refratio, exact);
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 8;
        Real offset = Real(4.0) * std::numeric_limits<Real>::epsilon();
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("offset", offset);
        }
        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox rb({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});

        int nerr = 0;

        // Walls a roundoff distance inside the domain faces, fluid inside.
        {
            Geometry geom(domain, rb, 0, {AMREX_D_DECL(0,0,0)});
            EB2::PlaneIF plo_x({AMREX_D_DECL(offset,0._rt,0._rt)}, {AMREX_D_DECL(-1._rt,0._rt,0._rt)});
            EB2::PlaneIF phi_x({AMREX_D_DECL(1._rt-offset,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,0._rt,0._rt)});
            EB2::PlaneIF plo_y({AMREX_D_DECL(0._rt,offset,0._rt)}, {AMREX_D_DECL(0._rt,-1._rt,0._rt)});
            EB2::PlaneIF phi_y({AMREX_D_DECL(0._rt,1._rt-offset,0._rt)}, {AMREX_D_DECL(0._rt,1._rt,0._rt)});
#if (AMREX_SPACEDIM == 3)
            EB2::PlaneIF plo_z({0._rt,0._rt,offset}, {0._rt,0._rt,-1._rt});
            EB2::PlaneIF phi_z({0._rt,0._rt,1._rt-offset}, {0._rt,0._rt,1._rt});
            auto walls = EB2::makeUnion(plo_x, phi_x, plo_y, phi_y, plo_z, phi_z);
#else
            auto walls = EB2::makeUnion(plo_x, phi_x, plo_y, phi_y);
#endif
            nerr += run_case("walls", walls, geom, 1);
            nerr += run_case("walls refratio 2", walls, geom, 2);
        }

        // Sphere with fluid inside.
        {
            Geometry geom(domain, rb, 0, {AMREX_D_DECL(0,0,0)});
            EB2::SphereIF sphere(0.3_rt, {AMREX_D_DECL(0.5_rt,0.5_rt,0.5_rt)}, true);
            nerr += run_case("sphere fluid inside", sphere, geom, 1);
        }

        // Sphere with fluid outside, periodic domain.
        {
            Geometry geom(domain, rb, 0, {AMREX_D_DECL(1,1,1)});
            EB2::SphereIF sphere(0.3_rt, {AMREX_D_DECL(0.5_rt,0.5_rt,0.5_rt)}, false);
            nerr += run_case("sphere fluid outside", sphere, geom, 1);
        }

        // Sphere centered on a node with a radius of 3 cells, so that it
        // touches nodes from inside.  Fluid outside.
        {
            Geometry geom(domain, rb, 0, {AMREX_D_DECL(0,0,0)});
            Real const r = Real(3.0)/Real(n_cell);
            EB2::SphereIF sphere(r, {AMREX_D_DECL(0.5_rt,0.5_rt,0.5_rt)}, false);
            auto exact = [=] (RealArray const& p) {
                Real d2 = 0.0_rt;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { d2 += (p[d]-0.5_rt)*(p[d]-0.5_rt); }
                return std::sqrt(d2) - r;
            };
            nerr += run_case("sphere through nodes", sphere, geom, 1, exact);
        }

        if (nerr == 0) {
            amrex::Print() << "SignedDistance test PASSED\n";
        } else {
            amrex::Abort("SignedDistance test FAILED");
        }
    }
    amrex::Finalize();
    return 0;
}
