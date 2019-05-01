
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

#include "Poisson.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int verbose = 1;
        int n_cell = 128;
        int max_grid_size = 32;
        int is_periodic = 1;

        // read parameters
        {
            ParmParse pp;
            pp.query("verbose", verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});

            // periodicity
            // note there is currently an AMReX issue with enclosed domains that do not touch the boundary
            // you must specifiy Dirichlet domain boundary conditions so AMReX doesn't think the
            // problem is singular
            Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

            Geometry::Setup(&rb, 0, is_periodic.data());

            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
            geom.define(domain);
            
            grids.define(domain); // define the BoxArray to be a single grid
            grids.maxSize(max_grid_size); // chop domain up into boxes with length max_Grid_size

            dmap.define(grids); // create a processor distribution mapping given the BoxARray
        }

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

        // charge density and electric potential
        MultiFab q  (grids, dmap, 1, 0, MFInfo(), factory);
        MultiFab phi(grids, dmap, 1, 0, MFInfo(), factory);

        q.setVal(0.0);
        InitData(q);

        LPInfo info;

        MLEBABecLap mlebabec({geom},{grids},{dmap},info,{&factory});

        // define array of LinOpBCType for domain boundary conditions
        std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
        std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;
	for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // see comment about boundary conditions above
            bc_lo[idim] = LinOpBCType::Dirichlet;
            bc_hi[idim] = LinOpBCType::Dirichlet;
        }

        // Boundary of the whole domain. This functions must be called,
        // and must be called before other bc functions.
        mlebabec.setDomainBC(bc_lo,bc_hi);

        // see AMReX_MLLinOp.H for an explanation
        mlebabec.setLevelBC(0, NULL);
    
        // operator looks like (ACoef - div BCoef grad) phi = rhs

        // set ACoef to zero
        MultiFab acoef(grids, dmap, 1, 0, MFInfo(), factory);
        acoef.setVal(0.);
        mlebabec.setACoeffs(0, acoef);
        
        // set BCoef to 1.0 (and array of face-centered coefficients)
        Array<MultiFab,AMREX_SPACEDIM> bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0, MFInfo(), factory);
            bcoef[idim].setVal(1.0);
        }
        mlebabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

        // scaling factors; these multiply ACoef and BCoef
        Real ascalar = 0.0;
        Real bscalar = 1.0;
        mlebabec.setScalars(ascalar, bscalar);

        // think of this beta as the "BCoef" associated with an EB face
        MultiFab beta(grids, dmap, 1, 0, MFInfo(), factory);
        beta.setVal(1.);

        // set homogeneous Dirichlet BC for EB
        mlebabec.setEBHomogDirichlet(0,beta);

        MLMG mlmg(mlebabec);

        // relative and absolute tolerances for linear solve
        const Real tol_rel = 1.e-10;
        const Real tol_abs = 0.0;

        mlmg.setVerbose(verbose);
        
        // Solve linear system
        phi.setVal(0.0); // initial guess for phi
        mlmg.solve({&phi}, {&q}, tol_rel, tol_abs);
        
        // store plotfile variables; q and phi
        MultiFab plotfile_mf(grids, dmap, 2, 0, MFInfo(), factory);
        MultiFab::Copy(plotfile_mf,  q,0,0,1,0);
        MultiFab::Copy(plotfile_mf,phi,0,1,1,0);
        
        EB_WriteSingleLevelPlotfile("plt", plotfile_mf, {"q", "phi"}, geom, 0.0, 0);
    }

    amrex::Finalize();
}
