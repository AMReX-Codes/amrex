
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Machine.H>

#ifdef AMREX_SOFT_PERF_COUNTERS
// need for perf counters
#include <AMReX_MLCellLinOp.H>
#endif

#include <prob_par.H>

using namespace amrex;

void init_prob_parms ();
void init_prob (const Vector<Geometry>& geom, Vector<MultiFab>& alpha, Vector<MultiFab>& beta,
                Vector<MultiFab>& rhs, Vector<MultiFab>& exact);
void solve_with_mlmg (const Vector<Geometry>& geom, int rr,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs, const Vector<MultiFab>& exact);
void write_plotfile (const Vector<Geometry>& geom, int rr,
                     const Vector<MultiFab>& soln, const Vector<MultiFab>& exact,
                     const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                     const Vector<MultiFab>& rhs);

namespace {
    void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids);
    BoxArray readBoxList (const std::string& file, Box& domain);
}

namespace {
    int max_level     = 1;
    int nlevels       = 2;
    int n_cell        = 64;
    int max_grid_size = 32;
    int ref_ratio     = 2;
    std::string boxes_file;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    {
        BL_PROFILE("main()");

        init_prob_parms();

        Vector<Geometry> geom;
        Vector<BoxArray> grids;
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
        
        solve_with_mlmg (geom, ref_ratio, soln, alpha, beta, rhs, exact);

#ifdef AMREX_SOFT_PERF_COUNTERS
        MLCellLinOp::perf_counters.print();
#endif

        write_plotfile (geom, ref_ratio, soln, exact, alpha, beta, rhs);
    }

    amrex::Finalize();
}

namespace {

void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids)
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_level", max_level);
    pp.query("max_grid_size", max_grid_size);
    pp.query("ref_ratio", ref_ratio);
    pp.query("boxes", boxes_file);

    if (!boxes_file.empty())
    {
        Box dmn;
        const BoxArray& ba = readBoxList(boxes_file, dmn);

        const BoxArray& uncovered = amrex::complementIn(dmn, ba);

        if (uncovered.empty() && dmn.isSquare() && dmn.smallEnd() == IntVect::TheZeroVector())
        {
            max_level = 0;
            nlevels = max_level + 1;
            n_cell = dmn.longside();

            geom.resize(nlevels);
            grids.resize(nlevels);

            grids[0] = ba;
        }
        else
        {
            amrex::Print() << "Add a coarse level\n";
            max_level = 1;
            nlevels = max_level + 1;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dmn.coarsenable(ref_ratio), "Domain must be coarsenable");
            dmn.coarsen(ref_ratio);
            dmn.setSmall(IntVect::TheZeroVector());
            n_cell = dmn.longside();

            geom.resize(nlevels);
            grids.resize(nlevels);

            IntVect dom0_lo {IntVect::TheZeroVector()};
            IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};
            Box dom0 {dom0_lo, dom0_hi};
            grids[0] = BoxArray{dom0};
            grids[0].maxSize(max_grid_size);

            grids[1] = ba;
        }
    }
    else
    {
        nlevels = max_level + 1;
        geom.resize(nlevels);
        grids.resize(nlevels);
        
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
    }
    
    std::array<Real,AMREX_SPACEDIM> prob_lo{AMREX_D_DECL(0.,0.,0.)};
    std::array<Real,AMREX_SPACEDIM> prob_hi{AMREX_D_DECL(1.,1.,1.)};
    RealBox real_box{prob_lo, prob_hi};
    
    const int coord = 0;  // Cartesian coordinates
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    if (prob::bc_type == MLLinOp::BCType::Periodic)
    {
        std::fill(is_periodic.begin(), is_periodic.end(), 1);
    }

    IntVect dom0_lo {IntVect::TheZeroVector()};
    IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};
    Box dom0 {dom0_lo, dom0_hi};
    
    geom[0].define(dom0, &real_box, coord, is_periodic.data());
    for (int ilev=1, n=grids.size(); ilev < n; ++ilev)
    {
        dom0.refine(ref_ratio);
        geom[ilev].define(dom0, &real_box, coord, is_periodic.data());
    }
}

BoxArray
readBoxList (const std::string& file, Box& domain)
{
    BoxList retval;

    std::ifstream boxspec;

    boxspec.open(file.c_str(), std::ios::in);

    if( !boxspec )
    {
        std::string msg = "readBoxList: unable to open ";
        msg += file;
        amrex::Error(msg.c_str());
    }
    boxspec >> domain;
    
    int numbox = 0;
    boxspec >> numbox;

    for ( int i=0; i<numbox; i++ )
    {
        Box tmpbox;
        boxspec >> tmpbox;
        if( !domain.contains(tmpbox) )
	{
            std::cerr << "readBoxList: bogus box " << tmpbox << '\n';
            exit(1);
        }
        retval.push_back(tmpbox);
    }

    return BoxArray(retval);
}

}

