
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids);
void init_prob (const Vector<Geometry>& geom, Vector<MultiFab>& alpha, Vector<MultiFab>& beta,
                Vector<MultiFab>& rhs, Vector<MultiFab>& exact);
void solve_with_mlmg (const Vector<Geometry>& geom,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs);
void write_plotfile (const Vector<Geometry>& geom, int rr,
                     const Vector<MultiFab>& soln, const Vector<MultiFab>& exact,
                     const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                     const Vector<MultiFab>& rhs);

namespace {
    static int max_level     = 1;
    static int nlevels       = 2;
    static int n_cell        = 64;
    static int max_grid_size = 32;
    static int ref_ratio     = 4;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    {
        ParmParse pp;
        pp.query("max_level", max_level);
        nlevels = max_level + 1;

        Vector<Geometry> geom(nlevels);
        Vector<BoxArray> grids(nlevels);
        build_geometry_and_grids(geom, grids);
        
        Vector<MultiFab> soln(nlevels);
        Vector<MultiFab> exact(nlevels);
        Vector<MultiFab> alpha(nlevels);
        Vector<MultiFab> beta(nlevels);
        Vector<MultiFab> rhs(nlevels);
        
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            DistributionMapping dm{grids[ilev]};
            soln [ilev].define(grids[ilev], dm, 1, 1);  // 1 ghost cell to store boundary conditions
            exact[ilev].define(grids[ilev], dm, 1, 0);
            alpha[ilev].define(grids[ilev], dm, 1, 0);
            beta [ilev].define(grids[ilev], dm, 1, 1);  // 1 ghost cell for averaging to faces
            rhs  [ilev].define(grids[ilev], dm, 1, 0);
        }

        init_prob (geom, alpha, beta, rhs, exact);

        for (auto& mf : soln) {
            mf.setVal(0.0); // initial guess
        }
        
        solve_with_mlmg (geom, soln, alpha, beta, rhs);

        write_plotfile (geom, ref_ratio, soln, exact, alpha, beta, rhs);
    }

    amrex::Finalize();
}

void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids)
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("ref_ratio", ref_ratio);

    IntVect dom0_lo {IntVect::TheZeroVector()};
    IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};

    Box dom0 {dom0_lo, dom0_hi};
    BoxArray ba0{dom0};
    
    grids[0] = ba0;
    grids[0].maxSize(max_grid_size);
    
    for (int ilev=1, n=grids.size(); ilev < n; ++ilev)
    {
        ba0.grow(-n_cell/4);
        ba0.refine(ref_ratio);
        grids[ilev] = ba0;
        grids[ilev].maxSize(max_grid_size);
    }

    std::array<Real,AMREX_SPACEDIM> prob_lo{AMREX_D_DECL(0.,0.,0.)};
    std::array<Real,AMREX_SPACEDIM> prob_hi{AMREX_D_DECL(1.,1.,1.)};
    RealBox real_box{prob_lo, prob_hi};

    const int coord = 0;  // Cartesian coordinates
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};  // Non-periodic for now

    geom[0].define(dom0, &real_box, coord, is_periodic.data());
    for (int ilev=1, n=grids.size(); ilev < n; ++ilev)
    {
        dom0.refine(ref_ratio);
        geom[ilev].define(dom0, &real_box, coord, is_periodic.data());
    }
}
