#include "MyTest.H"

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_EB2.H>
#include <AMReX_HypreSolver.H>

#include <cmath>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    initData();
}

void
MyTest::solve ()
{
    // In this example, we assume the domain boundary is homegeneous Dirichlet.

    auto const& nddom = amrex::surroundingNodes(geom.Domain());
    const auto dxi = geom.InvCellSizeArray();
    GpuArray<HYPRE_Real,AMREX_SPACEDIM> fac
        {AMREX_D_DECL(static_cast<HYPRE_Real>(dxi[0]*dxi[0]),
                      static_cast<HYPRE_Real>(dxi[1]*dxi[1]),
                      static_cast<HYPRE_Real>(dxi[2]*dxi[2]))};

    if (factory->isAllRegular())
    {
        HYPRE_Real fac0 = HYPRE_Real(-2.)*(AMREX_D_TERM(fac[0],+fac[1],+fac[2]));

        // Is variable n at (i,j,k) in Box boxno (local index) valid?
        // (i.e., not exactly on Dirichlet boundary)
        auto marker = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int /*n*/)
            -> bool
        {
            return nddom.strictly_contains(i,j,k);
        };

        // For variable n at (i,j,k) in Box boxno (local index), fill its row in
        // the matrix.
        // [in ] gid : gid[n] is the id for variable n at (i,j,k)
        // [out] ncols: # of columns in this row.
        // [out] cols: column indices in this row.
        // [out] mat : matrix elemens in this row.
        auto filler = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int n,
                                            Array4<HYPRE_Int const> const* gid,
                                            HYPRE_Int& ncols, HYPRE_Int* cols,
                                            HYPRE_Real* mat)
        {
            ncols = 0;
            if (i > nddom.smallEnd(0)+1) {
                cols[ncols] = gid[n](i-1,j,k);
                mat [ncols] = fac[0];
                ++ncols;
            }
            if (i < nddom.bigEnd(0)-1) {
                cols[ncols] = gid[n](i+1,j,k);
                mat [ncols] = fac[0];
                ++ncols;
            }
            if (j > nddom.smallEnd(1)+1) {
                cols[ncols] = gid[n](i,j-1,k);
                mat [ncols] = fac[1];
                ++ncols;
            }
            if (j < nddom.bigEnd(1)-1) {
                cols[ncols] = gid[n](i,j+1,k);
                mat [ncols] = fac[1];
                ++ncols;
            }
#if (AMREX_SPACEDIM > 2)
            if (k > nddom.smallEnd(2)+1) {
                cols[ncols] = gid[n](i,j,k-1);
                mat [ncols] = fac[2];
                ++ncols;
            }
            if (k < nddom.bigEnd(2)-1) {
                cols[ncols] = gid[n](i,j,k+1);
                mat [ncols] = fac[2];
                ++ncols;
            }
#endif
            cols[ncols] = gid[n](i,j,k);
            mat [ncols] = fac0;
            ++ncols;
        };

        constexpr int max_stencil_size = 2*AMREX_SPACEDIM+1;
        HypreSolver<max_stencil_size> hypre_solver
            ({IndexType::TheNodeType()}, IntVect(1), geom, grids, dmap,
             marker, filler, verbose);

        hypre_solver.solve(Vector<MultiFab*>{&phi}, Vector<MultiFab const*>{&rhs},
                           reltol, 0.0, max_iter);
    }
    else
    {
        auto const& levset_mf = factory->getLevelSet();
        auto const& levset_a = levset_mf.const_arrays();
        auto const& edgecent = factory->getEdgeCent();
        AMREX_D_TERM(auto const& ecx_a = edgecent[0]->const_arrays();,
                     auto const& ecy_a = edgecent[1]->const_arrays();,
                     auto const& ecz_a = edgecent[2]->const_arrays());
        auto const& rhs_a = rhs.arrays();

        // Is variable n at (i,j,k) in Box boxno (local index) valid?
        // (i.e., not exactly on Dirichlet boundary, not covered)
        auto marker = [=] AMREX_GPU_DEVICE (int boxno, int i, int j, int k, int /*n*/)
            -> bool
        {
            // level set < 0 means the node is in the fluid and valid.
            return nddom.strictly_contains(i,j,k) && (levset_a[boxno](i,j,k) < Real(0.));
        };

        auto peb = phi_eb;
        auto pdom = phi_domain;
        auto dlo = amrex::lbound(nddom);
        auto dhi = amrex::ubound(nddom);

        // For variable n at (i,j,k) in Box boxno (local index), fill its row in
        // the matrix.
        // [in ] gid    gid[n] is the id for variable n at (i,j,k)
        // [out] ncols  # of non-zero columns in this row.
        // [out] cols   array of indices of columns with a non-zero matrix element in this row.
        // [out] mat    array of (non-zero) matrix elements in this row.
        auto filler = [=] AMREX_GPU_DEVICE (int boxno, int i, int j, int k, int n,
                                            Array4<HYPRE_Int const> const* gid,
                                            HYPRE_Int& ncols, HYPRE_Int* cols,
                                            HYPRE_Real* mat)
        {
            ncols = 0;

            AMREX_D_TERM(auto const& ecx = ecx_a[boxno];,
                         auto const& ecy = ecy_a[boxno];,
                         auto const& ecz = ecz_a[boxno]);
            auto const& levset = levset_a[boxno];

            Real hp, hm;
            Real scale = Real(1.0);
            Real s0 = Real(0.0);
            Real srhs = Real(0.0);

            hp = (!ecx || ecx(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecx(i,j,k));
            hm = (!ecx || ecx(i-1,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecx(i-1,j,k));
            scale = std::min({scale, hp, hm});
            Real s = fac[0]*Real(2.0)/(hp+hm);

            if (levset(i+1,j,k) < Real(0.0)) {
                s0 -= s;
                if (i+1 < dhi.x) {
                    cols[ncols] = gid[n](i+1,j,k);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hp;
                srhs -= s/hp * peb;
            }

            if (levset(i-1,j,k) < Real(0.0)) {
                s0 -= s;
                if (i-1 > dlo.x) {
                    cols[ncols] = gid[n](i-1,j,k);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hm;
                srhs -= s/hm * peb;
            }

            hp = (!ecy || ecy(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecy(i,j,k));
            hm = (!ecy || ecy(i,j-1,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecy(i,j-1,k));
            scale = std::min({scale, hp, hm});
            s = fac[1]*Real(2.0)/(hp+hm);

            if (levset(i,j+1,k) < Real(0.0)) {
                s0 -= s;
                if (j+1 < dhi.y) {
                    cols[ncols] = gid[n](i,j+1,k);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hp;
                srhs -= s/hp * peb;
            }

            if (levset(i,j-1,k) < Real(0.0)) {
                s0 -= s;
                if (j-1 > dlo.y) {
                    cols[ncols] = gid[n](i,j-1,k);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hm;
                srhs -= s/hm * peb;
            }

#if (AMREX_SPACEDIM > 2)
            hp = (!ecz || ecz(i,j,k) == Real(1.0)) ? Real(1.0) : (Real(1.0)+Real(2.)*ecz(i,j,k));
            hm = (!ecz || ecz(i,j,k-1) == Real(1.0)) ? Real(1.0) : (Real(1.0)-Real(2.)*ecz(i,j,k-1));
            scale = std::min({scale, hp, hm});
            s = fac[2]*Real(2.0)/(hp+hm);

            if (levset(i,j,k+1) < Real(0.0)) {
                s0 -= s;
                if (k+1 < dhi.z) {
                    cols[ncols] = gid[n](i,j,k+1);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hp;
                srhs -= s/hp * peb;
            }

            if (levset(i,j,k-1) < Real(0.0)) {
                s0 -= s;
                if (k-1 > dlo.z) {
                    cols[ncols] = gid[n](i,j,k-1);
                    mat[ncols] = s;
                    ++ncols;
                } else {
                    srhs -= s * pdom;
                }
            } else {
                s0   -= s/hm;
                srhs -= s/hm * peb;
            }
#endif

            for (int icol = 0; icol < ncols; ++icol) {
                mat[icol] *= scale;
            }

            cols[ncols] = gid[n](i,j,k);
            mat [ncols] = s0*scale;
            ++ncols;

            rhs_a[boxno](i,j,k,n) += srhs;
            rhs_a[boxno](i,j,k,n) *= scale;
        };

        constexpr int max_stencil_size = 2*AMREX_SPACEDIM+1;
        HypreSolver<max_stencil_size> hypre_solver
            ({IndexType::TheNodeType()}, IntVect(1), geom, grids, dmap,
             marker, filler, verbose);

        hypre_solver.solve(Vector<MultiFab*>{&phi}, Vector<MultiFab const*>{&rhs},
                           reltol, 0.0, max_iter);
    }

    amrex::VisMF::Write(phi, "phi");
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("phi_domain", phi_domain);
    pp.query("phi_eb", phi_eb);

    pp.query("verbose", verbose);
    pp.query("max_iter", max_iter);
    pp.query("reltol", reltol);
}

void
MyTest::initGrids ()
{
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(0,0,0)};
    geom.define(domain0, &rb, 0, isperiodic.data());

    grids.define(domain0);
    grids.maxSize(max_grid_size);
}

void
MyTest::initData ()
{
    dmap.define(grids);
    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    factory = std::make_unique<EBFArrayBoxFactory>
        (eb_level, geom, grids, dmap, Vector<int>{1,1,1}, EBSupport::full);

    BoxArray const& nba = amrex::convert(grids, IntVect(1));

    phi.define(nba, dmap, 1, 0, MFInfo(), *factory);
    rhs.define(nba, dmap, 1, 0, MFInfo(), *factory);

    phi.setVal(0.0);
    rhs.setVal(1.0);
}
