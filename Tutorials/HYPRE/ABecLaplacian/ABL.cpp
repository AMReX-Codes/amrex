
#include <ABL.H>
#include <ABL_F.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

#include <string>

using namespace amrex;

ABL::ABL ()
{
    static_assert(AMREX_SPACEDIM == 3, "3D only");

    // runtime parameters
    {
        ParmParse pp;

        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);

        std::string bc_type_s{"Dirichlet"};
        pp.query("bc_type", bc_type_s);
        std::transform(bc_type_s.begin(), bc_type_s.end(), bc_type_s.begin(), ::tolower);
        if (bc_type_s == "dirichlet") {
            bc_type = amrex::LinOpBCType::Dirichlet;
        } else if (bc_type_s == "neumann") {
            bc_type = amrex::LinOpBCType::Neumann;            
        } else if (bc_type_s == "periodic") {
            bc_type = amrex::LinOpBCType::interior;
        } else {
            amrex::Abort("Unknown bc_type: "+bc_type_s);
        }
        
        pp.query("bc_value", bc_value);

        pp.query("tol_rel", tol_rel);
        pp.query("tol_abs", tol_abs);
        pp.query("maxiter", maxiter);
        
        pp.query("verbose", verbose);
    }
    
    BoxArray ba;
    {
        IntVect dom_lo(0,0,0);
        IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
        Box domain(dom_lo,dom_hi);
        ba.define(domain);
        ba.maxSize(max_grid_size);

        RealBox real_box({0.0,0.0,0.0}, {1.0,1.0,1.0});
        std::array<int,3> is_periodic {0,0,0};
        if (bc_type == amrex::LinOpBCType::interior) {
            is_periodic = {1,1,1};
        }
        geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
    }

    DistributionMapping dmap{ba};

    rhs.define(ba, dmap, 1, 0);
    init_rhs();

    alpha.define(ba, dmap, 1, 0);
    for (int idim = 0; idim < 3; ++idim) {
        BoxArray nba = ba;
        nba.surroundingNodes(idim);
        beta[idim].define(nba, dmap, 1, 0);
    }
    init_coeffs();

    soln.define(ba, dmap, 1, 0);
    the_soln.define(ba, dmap, 1, 0);

    comp_the_solution();
}

void
ABL::init_rhs ()
{
    const int ibnd = static_cast<int>(bc_type);
    const Real* dx = geom.CellSize();

    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_init_rhs(BL_TO_FORTRAN_BOX(tbx),
                      BL_TO_FORTRAN_ANYD(rhs[mfi]),
                      dx, &a, &b, &sigma, &w, &ibnd);
    }
}

void
ABL::init_coeffs ()
{
    alpha.setVal(1.0);

    const Real* dx = geom.CellSize();

    MultiFab betacc(alpha.boxArray(), alpha.DistributionMap(), 1, 1);

    for (MFIter mfi(betacc,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.growntilebox();
        fort_init_cc_coef(BL_TO_FORTRAN_BOX(tbx),
                          BL_TO_FORTRAN_ANYD(betacc[mfi]),
                          dx, &sigma, &w);
    }

    amrex::average_cellcenter_to_face({&beta[0], &beta[1], &beta[2]},
                                      betacc, geom);
}

void
ABL::comp_the_solution ()
{
    const int ibnd = static_cast<int>(bc_type);
    const Real* dx = geom.CellSize();

    for (MFIter mfi(the_soln); mfi.isValid(); ++mfi)
    {
        fort_comp_asol(BL_TO_FORTRAN_ANYD(the_soln[mfi]),
                       dx, &ibnd);
    }
}

void
ABL::solve ()
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();
        
    {
        BL_PROFILE("__solve__");    
        Hypre hypre_solver(ba, dm, geom, ParallelDescriptor::Communicator());
        hypre_solver.setScalars(a, b);
        hypre_solver.setACoeffs(alpha);
        hypre_solver.setBCoeffs({&beta[0],&beta[1],&beta[2]});
        hypre_solver.setVerbose(verbose);
        hypre_solver.solve(soln, rhs, tol_rel, tol_abs, maxiter, bc_type, bc_value);
    }        

    {
        MultiFab diff(ba, dm, 1, 0);
        MultiFab::Copy(diff, soln, 0, 0, 1, 0);
        MultiFab::Subtract(diff, the_soln, 0, 0, 1, 0);
        amrex::Print() << "\nMax-norm of the error is " << diff.norm0()
                       << "\nMaximum absolute value of the solution is " << the_soln.norm0()
                       << "\nMaximum absolute value of the rhs is " << rhs.norm0()
                       << "\n";
    }
}
