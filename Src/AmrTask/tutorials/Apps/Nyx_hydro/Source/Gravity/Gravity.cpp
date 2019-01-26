#include <cmath>

#include <AMReX_ParmParse.H>
#include "Gravity.H"
#include "Nyx.H"
#include <Gravity_F.H>
#include <Nyx_F.H>

#include <AMReX_MultiGrid.H>
#include <AMReX_Laplacian.H>
#include <AMReX_MacBndry.H>
#include <AMReX_LO_BCTYPES.H>

#ifdef USEHPGMG
#include <BL_HPGMG.H>
#endif

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

// MAX_LEV defines the maximum number of AMR levels allowed by the parent "Amr" object
#define MAX_LEV 15

// Give this a bogus default value to force user to define in inputs file
std::string Gravity::gravity_type = "fill_me_please";
int  Gravity::verbose       = 0;
int  Gravity::no_sync       = 0;
int  Gravity::no_composite  = 0;
int  Gravity::dirichlet_bcs = 0;
int  Gravity::solve_with_hpgmg = 0;
int  Gravity::solve_with_mlmg = 1;
int  Gravity::mlmg_max_fmg_iter = 0;
int  Gravity::mlmg_agglomeration = 0;
int  Gravity::mlmg_consolidation = 0;
Real Gravity::sl_tol        = 1.e-12;
Real Gravity::ml_tol        = 1.e-12;
Real Gravity::delta_tol     = 1.e-12;
Real Gravity::mass_offset   = 0;

extern "C"
{void fort_get_grav_const(Real* Gconst);}

// ************************************************************************** //

// Ggravity is defined as -4 * pi * G, where G is the gravitational constant.
// G is defined as Gconst in `fParallel/extern/constants/nyx_constants.f90` if
// NYX is defined in the GNUmakefile. G is defined as Gconst in
// `fParallel/extern/constants/constants.f90` if NYX is not defined in the
// GNUmakefile

// In CGS, this constant is currently
//      Gconst   =  6.67428e-8           cm^3/g/s^2 , which results in
//      Ggravity = -83.8503442814844e-8  cm^3/g/s^2

// In cosmological units, this constant is currently
//      Gconst   =   4.3   e-9  (km/s)^2 Mpc / Msun, which results in
//      Ggravity = -54.0354e-9  (km/s)^2 Mpc / Msun

// ************************************************************************** //

static Real Ggravity = 0;

Gravity::Gravity (Amr*   Parent,
                  int    _finest_level,
                  BCRec* _phys_bc,
                  int    _density)
  :
    parent(Parent),
    LevelData(MAX_LEV),
    grad_phi_curr(MAX_LEV),
    grad_phi_prev(MAX_LEV),
    phi_flux_reg(MAX_LEV),
    grids(Parent->boxArray()),
    dmap(Parent->DistributionMap()),
    level_solver_resnorm(MAX_LEV),
    phys_bc(_phys_bc)
{
     density = _density;
     read_params();
     finest_level_allocated = -1;
#ifdef CGRAV
     if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav" || gravity_type == "StaticGrav")
          make_mg_bc();
#else
     if(gravity_type == "PoissonGrav") make_mg_bc();
#endif
}


Gravity::~Gravity ()
{
    // nothing to see here.
}

void
Gravity::read_params ()
{
    static bool done = false;

    if (!done)
    {
        ParmParse pp("gravity");
        pp.get("gravity_type", gravity_type);

#ifdef CGRAV
        if (gravity_type != "PoissonGrav" && gravity_type != "CompositeGrav" && gravity_type != "StaticGrav")
        {
            std::cout << "Sorry -- dont know this gravity type" << std::endl;
            amrex::Abort("Options are PoissonGrav, CompositeGrav and StaticGrav");
        }
#else
        if (gravity_type != "PoissonGrav")
        {
            std::cout << "Sorry -- dont know this gravity type" << std::endl;
            amrex::Abort("Options are PoissonGrav");
        }
#endif

        pp.query("v", verbose);
        pp.query("no_sync", no_sync);
        pp.query("no_composite", no_composite);

        pp.query("dirichlet_bcs", dirichlet_bcs);

        pp.query("solve_with_hpgmg", solve_with_hpgmg);
        pp.query("solve_with_mlmg", solve_with_mlmg);
        pp.query("mlmg_max_fmg_iter", mlmg_max_fmg_iter);
        pp.query("mlmg_agglomeration", mlmg_agglomeration);
        pp.query("mlmg_consolidation", mlmg_consolidation);

        const int nflags = static_cast<int>(solve_with_hpgmg)
            +              static_cast<int>(solve_with_mlmg);
        if (nflags >= 2) {
          amrex::Error("Multiple gravity solvers selected.");
        }

#ifndef USEHPGMG
        if (solve_with_hpgmg)
          amrex::Error("To use the HPGMG solver you must compile with USE_HPGMG = TRUE");
#endif

        // Allow run-time input of solver tolerances
        pp.query("ml_tol", ml_tol);
        pp.query("sl_tol", sl_tol);
        pp.query("delta_tol", delta_tol);

        Real Gconst;
        fort_get_grav_const(&Gconst);
        Ggravity = -4.0 * M_PI * Gconst;
        if (verbose > 0)
        {
            amrex::Print() << "Getting Gconst from nyx_constants: " << Gconst
                      << '\n';
            amrex::Print() << "Using " << Ggravity << " for 4 pi G in Gravity.cpp "
                      << '\n';
        }
        done = true;
    }
}

void
Gravity::install_level (int       level,
                        AmrLevel* level_data_to_install)
{
    if (verbose > 1)
        amrex::Print() << "Installing Gravity level " << level << '\n';

    LevelData[level] = level_data_to_install;

    level_solver_resnorm[level] = 0;

#ifdef CGRAV
    if (gravity_type != "StaticGrav")
    {
#endif

    const auto& dm = level_data_to_install->DistributionMap();

    grad_phi_prev[level].resize(BL_SPACEDIM);
    for (int n=0; n<BL_SPACEDIM; ++n)
    {
        grad_phi_prev[level][n].reset(new MultiFab(level_data_to_install->getEdgeBoxArray(n),dm,1,1));
        grad_phi_prev[level][n]->setVal(0.);
    }

    grad_phi_curr[level].resize(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        grad_phi_curr[level][n].reset(new MultiFab(level_data_to_install->getEdgeBoxArray(n),dm,1,1));
        grad_phi_curr[level][n]->setVal(0.);
    }

    if (level > 0)
    {
        IntVect crse_ratio = parent->refRatio(level-1);
        phi_flux_reg[level].reset(new FluxRegister(level_data_to_install->boxArray(),
						   dm, crse_ratio, level, 1));
    }

#ifdef CGRAV
    }
#endif

    finest_level_allocated = level;
}

std::string
Gravity::get_gravity_type ()
{
    return gravity_type;
}

int
Gravity::get_no_sync ()
{
    return no_sync;
}

int
Gravity::get_no_composite ()
{
    return no_composite;
}

Vector<MultiFab*>
Gravity::get_grad_phi_prev (int level)
{
    return amrex::GetVecOfPtrs(grad_phi_prev[level]);
}

Vector<MultiFab*>
Gravity::get_grad_phi_curr (int level)
{
    return amrex::GetVecOfPtrs(grad_phi_curr[level]);
}

void
Gravity::plus_grad_phi_curr (int level, const Vector<MultiFab*>& addend)
{
    for (int n = 0; n < BL_SPACEDIM; n++)
        grad_phi_curr[level][n]->plus(*addend[n], 0, 1, 0);
}

void
Gravity::swap_time_levels (int level)
{

#ifdef CGRAV
    if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav")
#else
    if (gravity_type == "PoissonGrav")
#endif
    {
        for (int n=0; n < BL_SPACEDIM; n++)
        {
	    std::swap(grad_phi_prev[level][n], grad_phi_curr[level][n]);
            grad_phi_curr[level][n].reset(new MultiFab(BoxArray(grids[level]).surroundingNodes(n), 
						       dmap[level], 1, 1));
            grad_phi_curr[level][n]->setVal(1.e50);
        }
    }
}

void
Gravity::zero_phi_flux_reg (int level)
{
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    phi_flux_reg[level]->setVal(0);
}

void
Gravity::solve_for_old_phi (int               level,
                            MultiFab&         phi,
                            const Vector<MultiFab*>& grad_phi,
                            int               ngrow_for_solve,
                            int               fill_interior)
{
    BL_PROFILE("Gravity::solve_for_old_phi()");
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    if (verbose)
        amrex::Print() << "Gravity ... single level solve for old phi at level "
                  << level << std::endl;
    MultiFab Rhs(grids[level], dmap[level], 1, 0);
    Rhs.setVal(0.0);

#ifndef NO_HYDRO
    if (Nyx::Do_Hydro() == 1)
    {
       MultiFab&  S_old = LevelData[level]->get_old_data(State_Type);
       MultiFab::Copy(Rhs, S_old, density, 0, 1, 0);
    }
#endif

    AddParticlesToRhs(level,Rhs,ngrow_for_solve);

    // We shouldn't need to use virtual or ghost particles for old phi solves.

    const Real time  = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Gravity::solve_for_new_phi (int               level,
                            MultiFab&         phi,
                            const Vector<MultiFab*>& grad_phi,
                            int               fill_interior,
                            int               ngrow_for_solve)
{
    BL_PROFILE("Gravity::solve_for_new_phi()");
#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    if (verbose)
        amrex::Print() << "Gravity ... single level solve for new phi at level "
                  << level << std::endl;

    MultiFab Rhs(grids[level], dmap[level], 1, 0);
    Rhs.setVal(0.0);

#ifndef NO_HYDRO
    if (Nyx::Do_Hydro() == 1)
    {
       MultiFab& S_new = LevelData[level]->get_new_data(State_Type);
       MultiFab::Copy(Rhs, S_new, density, 0, 1, 0);
    }
#endif

    AddParticlesToRhs(level,Rhs,ngrow_for_solve);
    AddVirtualParticlesToRhs(level,Rhs,ngrow_for_solve);
    AddGhostParticlesToRhs(level,Rhs);

    const Real time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Gravity::solve_for_phi (int               level,
                        MultiFab&         Rhs,
                        MultiFab&         phi,
                        const Vector<MultiFab*>& grad_phi,
                        Real              time,
                        int               fill_interior)

{
    BL_PROFILE("Gravity::solve_for_phi()");
    if (verbose)
        amrex::Print() << " ... solve for phi at level " << level << '\n';

    // This is a correction for fully periodic domains only
    if (Geometry::isAllPeriodic())
        CorrectRhsUsingOffset(level,Rhs);

    Rhs.mult(Ggravity);

    Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(level));

    BL_ASSERT(cs != 0);

    // Here we divide by a for the Poisson solve.
    Rhs.mult(1 / cs->get_comoving_a(time));

#ifndef NDEBUG
    if (Rhs.contains_nan(0,1,0))
    {
        std::cout << "Rhs in solve_for_phi at level " << level << " has NaNs" << std::endl;
        amrex::Abort("");
    }
#endif

    // Need to set the boundary values here so they can get copied into "bndry"
    if (dirichlet_bcs) set_dirichlet_bcs(level,&phi);

    if (solve_with_mlmg)
    {
        solve_for_phi_with_mlmg(level, Rhs, phi, grad_phi, time);
        return;
    }

    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level], dmap[level], 1, geom);

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                     : IntVect::TheZeroVector();
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    if (level == 0)
    {
        bndry.setBndryValues(phi, src_comp, dest_comp, num_comp, *phys_bc);
    }
    else
    {
        MultiFab c_phi;
        get_crse_phi(level, c_phi, time);
        BoxArray crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes, dmap[level],
			      in_rad, out_rad, extent_rad, num_comp);
        crse_br.copyFrom(c_phi, c_phi.nGrow(), src_comp, dest_comp, num_comp);
        bndry.setBndryValues(crse_br, src_comp, phi, src_comp, dest_comp,
                             num_comp, crse_ratio, *phys_bc);
    }

    Vector<BoxArray> bav(1);
    bav[0] = phi.boxArray();
    Vector<DistributionMapping> dmv(1);
    dmv[0] = Rhs.DistributionMap();
    Vector<Geometry> fgeom(1);
    fgeom[0] = geom;

    Vector< Vector<Real> > xa(1);
    Vector< Vector<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0;
            xb[0][i] = 0;
        }
    }
    else
    {
        const Real* dx_crse = parent->Geom(level-1).CellSize();
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0.5 * dx_crse[i];
            xb[0][i] = 0.5 * dx_crse[i];
        }
    }

    if ( Geometry::isAllPeriodic() )
    {
        if ( parent->Geom(level).Domain().numPts() == grids[level].numPts() )
        {
            Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(level)));

            Real sum = nyx_level->vol_weight_sum(Rhs,false);

            Rhs.plus(-sum, 0, 1, 0);

            if (verbose)
                amrex::Print() << " ... subtracting " << sum
                               << " to ensure solvability " << '\n';
        }
    }

    Vector<MultiFab*> phi_p = { &phi };
    Vector<MultiFab*> Rhs_p = { &Rhs };

    const Real  tol     = sl_tol;
    const Real  abs_tol = 0.;

#ifdef USEHPGMG
    if (solve_with_hpgmg)
    {
        solve_with_HPGMG(level, phi, grad_phi, Rhs, tol, abs_tol);
    }
#endif
}

void
Gravity::solve_for_delta_phi (int                        crse_level,
                              int                        fine_level,
                              MultiFab&                  crse_rhs,
                              const Vector<MultiFab*>&         delta_phi,
                              const Vector<Vector<MultiFab*> >& grad_delta_phi)
{
    BL_PROFILE("Gravity::solve_for_delta_phi()");
    solve_for_delta_phi_with_mlmg(crse_level,fine_level,crse_rhs,delta_phi,grad_delta_phi);
}

void
Gravity::gravity_sync (int crse_level, int fine_level, int iteration, int ncycle,
                       const MultiFab& drho_and_drhoU, const MultiFab& dphi,
                       const Vector<MultiFab*>& grad_delta_phi_cc)
{
    BL_PROFILE("Gravity::gravity_sync()");
    BL_ASSERT(parent->finestLevel()>crse_level);

    if (verbose)
    {
        amrex::Print() << " ... gravity_sync at crse_level " << crse_level << '\n';
        amrex::Print() << " ...     up to finest_level     " << fine_level << '\n';
    }

    // Build Rhs for solve for delta_phi
    MultiFab crse_rhs(grids[crse_level], dmap[crse_level], 1, 0);
    MultiFab::Copy(crse_rhs, drho_and_drhoU, 0, 0, 1, 0);
    crse_rhs.mult(Ggravity);
    crse_rhs.plus(dphi, 0, 1, 0);

    const Geometry& crse_geom   = parent->Geom(crse_level);
    const Box&      crse_domain = crse_geom.Domain();

    // In the all-periodic case we enforce that CrseRhsSync sums to zero.
    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()))
    {
        Real local_correction = 0;
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:local_correction)
#endif
        for (MFIter mfi(crse_rhs,true); mfi.isValid(); ++mfi)
            local_correction += crse_rhs[mfi].sum(mfi.tilebox(), 0, 1);
        ParallelDescriptor::ReduceRealSum(local_correction);

        local_correction /= grids[crse_level].numPts();

        if (verbose)
            amrex::Print() << "WARNING: Adjusting RHS in gravity_sync solve by " << local_correction << '\n';

        crse_rhs.plus(-local_correction,0,1,0);
    }

    // delta_phi needs a ghost cell for the solve
    Vector<std::unique_ptr<MultiFab> >  delta_phi(fine_level - crse_level + 1);
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        delta_phi[lev-crse_level].reset(new MultiFab(grids[lev], dmap[lev], 1, 1));
        delta_phi[lev-crse_level]->setVal(0);
    }

    Vector <Vector<std::unique_ptr<MultiFab> > > ec_gdPhi(fine_level - crse_level + 1);
    for (int lev = crse_level; lev <= fine_level; lev++) {
        Nyx* Nyx_lev = dynamic_cast<Nyx*>(&parent->getLevel(lev));
        ec_gdPhi[lev-crse_level].resize(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
           ec_gdPhi[lev-crse_level][n].reset(new MultiFab(Nyx_lev->getEdgeBoxArray(n),
							  Nyx_lev->DistributionMap(),
							  1,0));
    }

    // Do multi-level solve for delta_phi
    solve_for_delta_phi(crse_level, fine_level, crse_rhs, 
			amrex::GetVecOfPtrs(delta_phi),
			amrex::GetVecOfVecOfPtrs(ec_gdPhi));

    crse_rhs.clear();

    // In the all-periodic case we enforce that delta_phi averages to zero.
    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()) ) {
       Real local_correction = 0.0;
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:local_correction)
#endif
       for (MFIter mfi(*delta_phi[0],true); mfi.isValid(); ++mfi) {
           local_correction += (*delta_phi[0])[mfi].sum(mfi.tilebox(),0,1);
       }
       ParallelDescriptor::ReduceRealSum(local_correction);

       local_correction = local_correction / grids[crse_level].numPts();

       for (int lev = crse_level; lev <= fine_level; lev++) {
           delta_phi[lev-crse_level]->plus(-local_correction,0,1,1);
       }
    }

    // Add delta_phi to phi_new, and grad(delta_phi) to grad(delta_phi_curr) on each level
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        LevelData[lev]->get_new_data(PhiGrav_Type).plus(*delta_phi[lev-crse_level], 0, 1, 0);
        for (int n = 0; n < BL_SPACEDIM; n++)
            grad_phi_curr[lev][n]->plus(*ec_gdPhi[lev-crse_level][n], 0, 1, 0);
    }

    int is_new = 1;

    // Average phi_new from fine to coarse level
    for (int lev = fine_level - 1; lev >= crse_level; lev--)
    {
        const IntVect& ratio = parent->refRatio(lev);
        amrex::average_down(LevelData[lev+1]->get_new_data(PhiGrav_Type),
                             LevelData[lev  ]->get_new_data(PhiGrav_Type),
                             0, 1, ratio);
    }

    // Average the edge-based grad_phi from finer to coarser level
    for (int lev = fine_level-1; lev >= crse_level; lev--)
        average_fine_ec_onto_crse_ec(lev, is_new);

    // Add the contribution of grad(delta_phi) to the flux register below if necessary.
    if (crse_level > 0 && iteration == ncycle)
    {
        for (MFIter mfi(*delta_phi[0]); mfi.isValid(); ++mfi)
            for (int n = 0; n < BL_SPACEDIM; ++n)
                phi_flux_reg[crse_level]->FineAdd((*ec_gdPhi[0][n])[mfi], n,
						  mfi.index(), 0, 0, 1, 1);
    }

    for (int lev = crse_level; lev <= fine_level; lev++) {
        grad_delta_phi_cc[lev-crse_level]->setVal(0.0);
        const Geometry& geom = parent->Geom(lev);
        amrex::average_face_to_cellcenter(*grad_delta_phi_cc[lev-crse_level],
                                           amrex::GetVecOfConstPtrs(ec_gdPhi[lev-crse_level]),
                                           geom);
    }
}

void
Gravity::get_crse_phi (int       level,
                       MultiFab& phi_crse,
                       Real      time)
{
    BL_PROFILE("Gravity::get_crse_phi()");
    BL_ASSERT(level != 0);

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);

    phi_crse.clear();
    phi_crse.define(grids[level-1], dmap[level-1], 1, 1);

    // BUT NOTE we don't trust phi's ghost cells.
    FArrayBox phi_crse_temp;

    // Note that we must do these cases separately because it's possible to do a
    //   new solve after a regrid when the old data on the coarse grid may not yet
    //   be defined.
    for (MFIter mfi(phi_crse,true); mfi.isValid(); ++mfi)
    {
        const Box& gtbx = mfi.growntilebox();

        phi_crse_temp.resize(gtbx,1);

        if (fabs(alpha-1.0) < 1.e-15)
        {
            phi_crse[mfi].copy(LevelData[level-1]->get_new_data(PhiGrav_Type)[mfi]);
        }
        else if (fabs(alpha) < 1.e-15)
        {
            phi_crse[mfi].copy(LevelData[level-1]->get_old_data(PhiGrav_Type)[mfi]);
        }
        else
        {
            phi_crse_temp.copy(LevelData[level-1]->get_old_data(PhiGrav_Type)[mfi]);
            Real omalpha = 1.0 - alpha;
            phi_crse_temp.mult(omalpha);

            phi_crse[mfi].copy(LevelData[level-1]->get_new_data(PhiGrav_Type)[mfi],gtbx);
            phi_crse[mfi].mult(alpha,gtbx);
            phi_crse[mfi].plus(phi_crse_temp);
        }
    }

    const Geometry& geom = parent->Geom(level-1);
    phi_crse.FillBoundary(geom.periodicity());
}

void
Gravity::get_crse_grad_phi (int               level,
                            Vector<std::unique_ptr<MultiFab> >& grad_phi_crse,
                            Real              time)
{
    BL_PROFILE("Gravity::get_crse_grad_phi()");
    BL_ASSERT(level!=0);

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);
    const Real omalpha = 1.0 - alpha;

    Nyx* Nyx_crse_lev = dynamic_cast<Nyx*>(&parent->getLevel(level-1));

    BL_ASSERT(grad_phi_crse.size() == BL_SPACEDIM);

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        BL_ASSERT(!grad_phi_crse[i]);
        grad_phi_crse[i].reset(new MultiFab(Nyx_crse_lev->getEdgeBoxArray(i), 
					    Nyx_crse_lev->DistributionMap(),
					    1, 0));

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox grad_phi_crse_temp;
            for (MFIter mfi(*grad_phi_crse[i],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                grad_phi_crse_temp.resize(tbx,1);

                grad_phi_crse_temp.copy((*grad_phi_prev[level-1][i])[mfi]);
                grad_phi_crse_temp.mult(omalpha);

                (*grad_phi_crse[i])[mfi].copy((*grad_phi_curr[level-1][i])[mfi],tbx);
                (*grad_phi_crse[i])[mfi].mult(alpha,tbx);
                (*grad_phi_crse[i])[mfi].plus(grad_phi_crse_temp);
            }
        }
    }
}

void
Gravity::multilevel_solve_for_new_phi (int level,
                                       int finest_level,
                                       int ngrow_for_solve,
                                       int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::multilevel_solve_for_new_phi()");
    if (verbose)
        amrex::Print() << "Gravity ... multilevel solve for new phi at base level " << level
                       << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_curr[lev].size()==BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_curr[lev][n].reset(new MultiFab(eba, dmap[lev], 1, 1));
        }
    }

    int is_new = 1;
    actual_multilevel_solve(level, finest_level, 
			    amrex::GetVecOfVecOfPtrs(grad_phi_curr),
                            is_new, ngrow_for_solve, use_previous_phi_as_guess);
}

void
Gravity::multilevel_solve_for_old_phi (int level,
                                       int finest_level,
                                       int ngrow,
                                       int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::multilevel_solve_for_old_phi()");
    if (verbose)
        amrex::Print() << "Gravity ... multilevel solve for old phi at base level " << level
                       << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_prev[lev].size() == BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_prev[lev][n].reset(new MultiFab(eba, dmap[lev], 1, 1));
        }
    }

    int is_new  = 0;
    actual_multilevel_solve(level, finest_level,
			    amrex::GetVecOfVecOfPtrs(grad_phi_prev),
                            is_new, ngrow, use_previous_phi_as_guess);
}

void
Gravity::actual_multilevel_solve (int                       level,
                                  int                       finest_level,
                                  const Vector<Vector<MultiFab*> >& grad_phi,
                                  int                       is_new,
                                  int                       ngrow_for_solve,
                                  int                       use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::actual_multilevel_solve()");

    const int num_levels = finest_level - level + 1;

    Vector<MultiFab*> phi_p(num_levels);
    Vector<std::unique_ptr<MultiFab> > Rhs_p(num_levels);

    Vector<std::unique_ptr<MultiFab> > Rhs_particles(num_levels);
    for (int lev = 0; lev < num_levels; lev++)
    {
	Rhs_particles[lev].reset(new MultiFab(grids[level+lev], dmap[level+lev], 1, 0));
        Rhs_particles[lev]->setVal(0.);
    }

    const auto& rpp = amrex::GetVecOfPtrs(Rhs_particles);
    AddParticlesToRhs(level,finest_level,ngrow_for_solve,rpp);
    AddGhostParticlesToRhs(level,rpp);
    AddVirtualParticlesToRhs(finest_level,rpp);

    Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(level));

    BL_ASSERT(cs != 0);

    Real time = 0;

    if (is_new == 1)
    {
        time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    }
    else if (is_new == 0)
    {
        time = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    }

    // Here we get comoving_a b/c the RHS should be 4 * pi * G * density / a
    const Real a_inverse = 1. / (cs->get_comoving_a(time));

// *****************************************************************************

    for (int lev = 0; lev < num_levels; lev++)
    {
        if (is_new == 0)
        {
           // Working in result data structure directly
           phi_p[lev] = &LevelData[level+lev]->get_old_data(PhiGrav_Type);
        }
        else
        {
           // Working in result data structure directly
           phi_p[lev] = &LevelData[level+lev]->get_new_data(PhiGrav_Type);
        }

        if (!use_previous_phi_as_guess)
            phi_p[lev]->setVal(0);

        // Need to set the boundary values before "bndry" is defined so they get copied in
        if (dirichlet_bcs) set_dirichlet_bcs(level+lev,phi_p[lev]);

        Rhs_p[lev].reset(new MultiFab(grids[level+lev], dmap[level+lev], 1, 0));
        Rhs_p[lev]->setVal(0.0);

#ifndef NO_HYDRO
        if (Nyx::Do_Hydro() == 1)
        {
            if (is_new == 1)
            {
                MultiFab::Copy(*Rhs_p[lev], LevelData[level+lev]->get_new_data(State_Type), density, 0, 1, 0);
            }
            else if (is_new == 0)
            {
                MultiFab::Copy(*Rhs_p[lev], LevelData[level+lev]->get_old_data(State_Type), density, 0, 1, 0);
            }
        }
#endif
        MultiFab::Add(*Rhs_p[lev], *Rhs_particles[lev], 0, 0, 1, 0);
    }

    // Average phi from fine to coarse level before the solve.
    for (int lev = num_levels-1; lev > 0; lev--)
    {
        amrex::average_down(*phi_p[lev], *phi_p[lev-1],
                             0, 1, parent->refRatio(level+lev-1));
    }

// *****************************************************************************

    // This correction is for fully periodic domains only.
    if (Geometry::isAllPeriodic())
    {
        if (verbose)
            amrex::Print() << " ... subtracting average density " << mass_offset
                           << " from RHS at each level " << '\n';

        for (int lev = 0; lev < num_levels; lev++)
            for (MFIter mfi(*Rhs_p[lev]); mfi.isValid(); ++mfi)
                (*Rhs_p[lev])[mfi].plus(-mass_offset);

       // This is used to enforce solvability if appropriate.
       if ( parent->Geom(level).Domain().numPts() == grids[level].numPts() )
       {
           Real sum = 0;
           for (int lev = 0; lev < num_levels; lev++)
           {
               Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(level+lev)));
               sum += nyx_level->vol_weight_sum(*Rhs_p[lev],true);
           }

            sum /= parent->Geom(0).ProbSize();

            const Real eps = 1.e-10 * std::abs(mass_offset);
            if (std::abs(sum) > eps)
            {
               amrex::Print() << " ... current avg differs from mass_offset by " << sum << " " << '\n';
               amrex::Print() << " ... Gravity::actual_multilevel_solve -- total mass has changed!" << '\n';;
            }

            if (verbose)
                amrex::Print() << " ... subtracting " << sum << " to ensure solvability " << '\n';

            for (int lev = 0; lev < num_levels; lev++)
                (*Rhs_p[lev]).plus(-sum, 0, 1, 0);
       }
    }

// *****************************************************************************

    for (int lev = 0; lev < num_levels; lev++)
    {
        Rhs_p[lev]->mult(Ggravity, 0, 1);
        Rhs_p[lev]->mult(a_inverse);
    }

// *****************************************************************************

    if (solve_with_mlmg)
    {
        const MultiFab* crse_bcdata = nullptr;
        MultiFab CPhi;
        if (level > 0) {
            get_crse_phi(level, CPhi, time);
            crse_bcdata = &CPhi;
        }
        Real rel_eps = ml_tol;
        Real abs_eps = 0.;
        Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phi_aa;
        for (int amrlev = level; amrlev <= finest_level; ++amrlev) {
            grad_phi_aa.push_back({AMREX_D_DECL(grad_phi[amrlev][0],
                                                grad_phi[amrlev][1],
                                                grad_phi[amrlev][2])});
        }
        solve_with_MLMG(level, finest_level, phi_p, amrex::GetVecOfConstPtrs(Rhs_p),
                        grad_phi_aa, crse_bcdata, rel_eps, abs_eps);
    }

    // Average grad_phi from fine to coarse level
    for (int lev = finest_level; lev > level; lev--)
        average_fine_ec_onto_crse_ec(lev-1,is_new);
}

void
Gravity::get_old_grav_vector (int       level,
                              MultiFab& grav_vector,
                              Real      time)
{
    BL_PROFILE("Gravity::get_old_grav_vector()");
    // Set to zero to fill ghost cells.
    grav_vector.setVal(0);

#ifdef CGRAV
    if (gravity_type == "StaticGrav")
    {
        make_prescribed_grav(level,time,grav_vector,0);
        grav_vector.FillBoundary();
    }
    else
    {
#endif

    // Fill grow cells in grad_phi, will need to compute grad_phi_cc in 1 grow cell
    const Geometry& geom = parent->Geom(level);
#if 0
    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM ; i++)
           grad_phi_prev[level][i]->setBndry(0.);
    }
    else
    {
        Vector<std::unique_ptr<MultiFab> > crse_grad_phi(BL_SPACEDIM);
        get_crse_grad_phi(level, crse_grad_phi, time);
        fill_ec_grow(level, amrex::GetVecOfPtrs(grad_phi_prev[level]),
	                    amrex::GetVecOfPtrs(crse_grad_phi));
    }
#endif

    // Fill boundary values at the current level
    for (int i = 0; i < BL_SPACEDIM ; i++)
       grad_phi_prev[level][i]->FillBoundary(geom.periodicity());

    // Average edge-centered gradients to cell centers.
    amrex::average_face_to_cellcenter(grav_vector, 
				      amrex::GetVecOfConstPtrs(grad_phi_prev[level]),
				      geom);

#ifdef CGRAV
    if (gravity_type == "CompositeGrav")
    {
        make_prescribed_grav(level,time,grav_vector,1);
    }
#endif

    grav_vector.FillBoundary(geom.periodicity());

#ifdef CGRAV
    }
#endif

    MultiFab& G_old = LevelData[level]->get_old_data(Gravity_Type);

    // Fill G_old from grav_vector
    MultiFab::Copy(G_old, grav_vector, 0, 0, BL_SPACEDIM, 0);

    // This is a hack-y way to fill the ghost cell values of grav_vector
    //   before returning it
    // Note that this fills ghost cells over the coarse grid from interpolation, 
    //  not from the ghost cell values previously filled after the fill_ec_grow stuff.
    AmrLevel* amrlev = &parent->getLevel(level);
    int ng = grav_vector.nGrow();
    AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
}

void
Gravity::get_new_grav_vector (int       level,
                              MultiFab& grav_vector,
                              Real      time)
{
    BL_PROFILE("Gravity::get_new_grav_vector()");
#ifdef CGRAV
    if (gravity_type == "PoissonGrav" || gravity_type == "CompositeGrav")
#else
    if (gravity_type == "PoissonGrav")
#endif
    {
        // Set to zero to fill ghost cells
        grav_vector.setVal(0);

        // Fill grow cells in `grad_phi`, will need to compute `grad_phi_cc` in
        // 1 grow cell
        const Geometry& geom = parent->Geom(level);
#if 0
        if (level == 0)
        {
            for (int i = 0; i < BL_SPACEDIM ; i++)
            {
                grad_phi_curr[level][i]->setBndry(0);
                grad_phi_curr[level][i]->FillBoundary(geom.periodicity());
            }
        }
        else
        {
            Vector<std::unique_ptr<MultiFab> > crse_grad_phi(BL_SPACEDIM);
            get_crse_grad_phi(level, crse_grad_phi, time);
            fill_ec_grow(level, amrex::GetVecOfPtrs(grad_phi_curr[level]),
                                amrex::GetVecOfPtrs(crse_grad_phi));
        }
#endif

        for (int i = 0; i < BL_SPACEDIM ; i++)
            grad_phi_curr[level][i]->FillBoundary(geom.periodicity());

        // Average edge-centered gradients to cell centers, excluding grow cells
        amrex::average_face_to_cellcenter(grav_vector,
					   amrex::GetVecOfConstPtrs(grad_phi_curr[level]),
					   geom);

#ifdef CGRAV
        if (gravity_type == "CompositeGrav")
        {
          make_prescribed_grav(level,time,grav_vector,1);
        }
#endif

        grav_vector.FillBoundary(geom.periodicity());
    }

#ifdef CGRAV
    else if ( gravity_type == "StaticGrav")
    {
      make_prescribed_grav(level,time,grav_vector,0);
      grav_vector.FillBoundary();
    }
#endif

    MultiFab& G_new = LevelData[level]->get_new_data(Gravity_Type);

    // Fill G_new from grav_vector
    MultiFab::Copy(G_new, grav_vector, 0, 0, BL_SPACEDIM, 0);

    // This is a hack-y way to fill the ghost cell values of grav_vector
    //   before returning it
    // Note that this fills ghost cells over the coarse grid from interpolation, 
    //  not from the ghost cell values previously filled after the fill_ec_grow stuff.
    AmrLevel* amrlev = &parent->getLevel(level) ;
    int ng = grav_vector.nGrow();
    AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
}

void
Gravity::add_to_fluxes(int level, int iteration, int ncycle)
{
    BL_PROFILE("Gravity::add_to_fluxes()");
    const int       finest_level      = parent->finestLevel();
    FluxRegister*   phi_fine          = (level<finest_level ? phi_flux_reg[level+1].get() : nullptr);
    FluxRegister*   phi_current       = (level>0 ? phi_flux_reg[level].get() : nullptr);
    const Geometry& geom              = parent->Geom(level);
    const Real*     dx                = geom.CellSize();
    const Real      area[BL_SPACEDIM] = { dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1] };

#ifdef CGRAV
    if (gravity_type == "StaticGrav")
        return;
#endif

    int ngrow;
    fort_get_method_params(&ngrow);

    if (phi_fine)
    {
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            BoxArray ba = grids[level];
            ba.surroundingNodes(n);
            MultiFab fluxes(ba, dmap[level], 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(fluxes,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                FArrayBox& gphi_flux = fluxes[mfi];
                gphi_flux.copy((*grad_phi_curr[level][n])[mfi],tbx);
                gphi_flux.mult(area[n],tbx,0,1);
            }
            phi_fine->CrseInit(fluxes, n, 0, 0, 1, -1);
        }
    }

    if (phi_current && (iteration == ncycle))
    {
      MultiFab& phi_curr = LevelData[level]->get_new_data(PhiGrav_Type);
      for (MFIter mfi(phi_curr); mfi.isValid(); ++mfi)
      {
         for (int n=0; n<BL_SPACEDIM; ++n)
	     phi_current->FineAdd((*grad_phi_curr[level][n])[mfi],n,mfi.index(),0,0,1,area[n]);
      }
    }
}

void
Gravity::average_fine_ec_onto_crse_ec(int level, int is_new)
{
    BL_PROFILE("Gravity::average_fine_ec_to_crse_ec()");
    //
    // NOTE: this is called with level == the coarser of the two levels involved.
    //
    if (level == parent->finestLevel()) return;
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_gphi_fine_BA(grids[level+1].size());

    IntVect fine_ratio = parent->refRatio(level);

    for (int i = 0; i < crse_gphi_fine_BA.size(); ++i)
        crse_gphi_fine_BA.set(i, amrex::coarsen(grids[level+1][i],
                                                 fine_ratio));

    Vector<std::unique_ptr<MultiFab> > crse_gphi_fine(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        const BoxArray eba = BoxArray(crse_gphi_fine_BA).surroundingNodes(n);
        crse_gphi_fine[n].reset(new MultiFab(eba, dmap[level+1], 1, 0));
    }

    auto& grad_phi = (is_new) ? grad_phi_curr : grad_phi_prev;

    amrex::average_down_faces(amrex::GetVecOfConstPtrs(grad_phi[level+1]),
			       amrex::GetVecOfPtrs(crse_gphi_fine), fine_ratio);

    const Geometry& cgeom = parent->Geom(level);

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
	grad_phi[level][n]->copy(*crse_gphi_fine[n], cgeom.periodicity());
    }
}

void
Gravity::reflux_phi (int       level,
                     MultiFab& dphi)
{
    BL_PROFILE("Gravity::reflux()");
    const Geometry& geom = parent->Geom(level);
    dphi.setVal(0);
    phi_flux_reg[level+1]->Reflux(dphi, 1.0, 0, 0, 1, geom);
}

void
Gravity::fill_ec_grow (int                     level,
                       const Vector<MultiFab*>& ecF,
                       const Vector<MultiFab*>& ecC) const
{
    BL_PROFILE("Gravity::fill_ec_grow()");
    //
    // Fill grow cells of the edge-centered mfs.  Assume
    // ecF built on edges of grids at this amr level, and ecC
    // is build on edges of the grids at amr level-1
    //
    BL_ASSERT(ecF.size() == BL_SPACEDIM);

    const int nGrow = ecF[0]->nGrow();

    if (nGrow == 0 || level == 0) return;

    BL_ASSERT(nGrow == ecF[1]->nGrow());
    BL_ASSERT(nGrow == ecF[2]->nGrow());

    const BoxArray& fgrids = grids[level];
    const Geometry& fgeom  = parent->Geom(level);

    BoxList bl = amrex::GetBndryCells(fgrids, 1);

    BoxArray f_bnd_ba(bl);

    bl.clear();

    BoxArray c_bnd_ba   = BoxArray(f_bnd_ba.size());
    IntVect  crse_ratio = parent->refRatio(level-1);

    for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
        c_bnd_ba.set(i, Box(f_bnd_ba[i]).coarsen(crse_ratio));
        f_bnd_ba.set(i, Box(c_bnd_ba[i]).refine(crse_ratio));
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        //
        // crse_src & fine_src must have same parallel distribution.
        // We'll use the KnapSack distribution for the fine_src_ba.
        // Since fine_src_ba should contain more points, this'll lead
        // to a better distribution.
        //
        BoxArray crse_src_ba(c_bnd_ba);
        BoxArray fine_src_ba(f_bnd_ba);

        crse_src_ba.surroundingNodes(n);
        fine_src_ba.surroundingNodes(n);

        Vector<long> wgts(fine_src_ba.size());

        for (unsigned int i = 0; i < wgts.size(); i++)
        {
            wgts[i] = fine_src_ba[i].numPts();
        }
        DistributionMapping dm;
        //
        // This call doesn't invoke the MinimizeCommCosts() stuff.
        // There's very little to gain with these types of coverings
        // of trying to use SFC or anything else.
        // This also guarantees that these DMs won't be put into the
        // cache, as it's not representative of that used for more
        // usual MultiFabs.
        //
        dm.KnapSackProcessorMap(wgts, ParallelDescriptor::NProcs());

        MultiFab crse_src; crse_src.define(crse_src_ba, dm, 1, 0);
        MultiFab fine_src; fine_src.define(fine_src_ba, dm, 1, 0);

        crse_src.setVal(1.e200);
        fine_src.setVal(1.e200);
        //
        // We want to fill crse_src from ecC[n].
        //
	crse_src.copy(*ecC[n]);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
        {
            const int nComp = 1;
            const Box box = crse_src[mfi].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_PC_EDGE_INTERP, fort_pc_edge_interp)
                (box.loVect(), box.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(crse_src[mfi]), BL_TO_FORTRAN(fine_src[mfi]));
        }

        crse_src.clear();
        //
        // Replace pc-interpd fine data with preferred u_mac data at
        // this level u_mac valid only on surrounding faces of valid
        // region - this op will not fill grow region.
        //
        fine_src.copy(*ecF[n]); // parallel copy

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
        {
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine edges overlaying coarse edges.
            //
            const int nComp = 1;
            const Box& fbox = fine_src[mfi.index()].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_EDGE_INTERP, fort_edge_interp)
                (fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(fine_src[mfi]));
        }
        //
        // Build a mf with no grow cells on ecF grown boxes, do parallel copy into.
        //
        BoxArray fgridsG = ecF[n]->boxArray();
        fgridsG.grow(ecF[n]->nGrow());

        MultiFab ecFG(fgridsG, ecF[n]->DistributionMap(), 1, 0);

        ecFG.copy(fine_src); // Parallel copy

        fine_src.clear();

        ecFG.copy(*ecF[n]);   // Parallel copy

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*ecF[n]); mfi.isValid(); ++mfi)
            (*ecF[n])[mfi].copy(ecFG[mfi]);
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        ecF[n]->FillBoundary();
	ecF[n]->EnforcePeriodicity(fgeom.periodicity());
    }
}

void
Gravity::make_mg_bc ()
{
    BL_PROFILE("Gravity::make_mg_bc()");
    const Geometry& geom = parent->Geom(0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (phys_bc->lo(idim) == Symmetry) {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            }
            if (phys_bc->hi(idim) == Symmetry) {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            }
        }
    }
}

void
Gravity::set_mass_offset (Real time)
{
    BL_PROFILE("Gravity::set_mass_offset()");
    Real old_mass_offset = 0;

    int flev = parent->finestLevel();

    while (parent->getAmrLevels()[flev] == nullptr)
        flev--;

    if (flev > 0) old_mass_offset = mass_offset;

    mass_offset = 0;

    const Geometry& geom = parent->Geom(0);

    if (geom.isAllPeriodic())
    {
        //
        // Note: we must loop over levels because the volWgtSum routine zeros out
        //       crse regions under fine regions.
        //
        for (int lev = 0; lev <= flev; lev++)
        {
            Nyx* cs = dynamic_cast<Nyx*>(&parent->getLevel(lev));

            BL_ASSERT(cs != 0);

#ifndef NO_HYDRO
            if (Nyx::Do_Hydro() == 1)
                mass_offset += cs->vol_weight_sum("density",time,true);
#endif
            for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                mass_offset += Nyx::theActiveParticles()[i]->sumParticleMass(lev);
        }

        mass_offset /= geom.ProbSize();

        if (verbose)
            amrex::Print() << "Gravity ... defining average density in Gravity::set_mass_offset to be " 
                           << mass_offset << '\n';
    }

    const Real diff = std::abs(mass_offset - old_mass_offset);
    const Real eps  = 1.e-10 * std::abs(old_mass_offset);
    if (diff > eps && old_mass_offset > 0)
    {
        amrex::Print() << " ... new vs old mass_offset " << mass_offset << " "
                       << old_mass_offset << " ... diff is " << diff <<  '\n';
        amrex::Print() << " ... Gravity::set_mass_offset -- total mass has changed!"
                       << '\n';;
    }
}

void
Gravity::set_dirichlet_bcs (int       level,
                            MultiFab* phi)
{
    BL_PROFILE("Gravity::set_dirichlet_bcs()");
    const Real* dx        = parent->Geom(level).CellSize();
    const int*  domain_lo = parent->Geom(level).Domain().loVect();
    const int*  domain_hi = parent->Geom(level).Domain().hiVect();

    // Set phi to zero on all the ghost cells outside the domain.
    // If homogeneous bc's then we stop here
    for (MFIter mfi(*phi); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        BL_FORT_PROC_CALL(FORT_SET_HOMOG_BCS, fort_set_homog_bcs)
            (lo, hi, domain_lo, domain_hi, BL_TO_FORTRAN((*phi)[mfi]), dx);
    }
}

#ifdef CGRAV
void
Gravity::make_prescribed_grav (int       level,
                               Real      time,
                               MultiFab& grav_vector,
                               int       addToExisting)
{
    BL_PROFILE("Gravity::make_prescribed_grav()");
    const Geometry& geom = parent->Geom(level);
    const Real*     dx   = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();

       BL_FORT_PROC_CALL(FORT_PRESCRIBE_GRAV,fort_prescribe_grav)
           (bx.loVect(), bx.hiVect(), dx,
            BL_TO_FORTRAN(grav_vector[mfi]),
            geom.ProbLo(), &addToExisting);
    }
}
#endif

void
Gravity::AddParticlesToRhs (int               level,
                            MultiFab&         Rhs,
                            int               ngrow)
{
    BL_PROFILE("Gravity::AddParticlesToRhs()");
    // Use the same multifab for all particle types
    MultiFab particle_mf(grids[level], dmap[level], 1, ngrow);

    for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
    {
        particle_mf.setVal(0.);
        Nyx::theActiveParticles()[i]->AssignDensitySingleLevel(particle_mf, level);
        MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
    }
}

void
Gravity::AddParticlesToRhs(int base_level, int finest_level, int ngrow, const Vector<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddParticlesToRhsML()");
    const int num_levels = finest_level - base_level + 1;
    for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
    {
        Vector<std::unique_ptr<MultiFab> > PartMF;
        Nyx::theActiveParticles()[i]->AssignDensity(PartMF, base_level, 1, finest_level, ngrow);
        for (int lev = 0; lev < num_levels; lev++)
        {
            if (PartMF[lev]->contains_nan())
            {
                std::cout << "Testing particle density of type " << i << " at level " << base_level+lev << std::endl;
                amrex::Abort("...PartMF has NaNs in Gravity::actual_multilevel_solve()");
            }
        }

        for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        {
            amrex::average_down(*PartMF[lev+1], *PartMF[lev],
                                 0, 1, parent->refRatio(lev+base_level));
        }

        for (int lev = 0; lev < num_levels; lev++)
        {
            MultiFab::Add(*Rhs_particles[lev], *PartMF[lev], 0, 0, 1, 0);
        }
    }

}

void
Gravity::AddVirtualParticlesToRhs (int               level,
                                   MultiFab&         Rhs,
                                   int               ngrow)
{
    BL_PROFILE("Gravity::AddVirtualParticlesToRhs()");
    if (level <  parent->finestLevel())
    {
        // If we have virtual particles, add their density to the single level solve
        MultiFab particle_mf(grids[level], dmap[level], 1, ngrow);

        for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
        {
            particle_mf.setVal(0.);
            Nyx::theVirtualParticles()[i]->AssignDensitySingleLevel(particle_mf, level, 1, 1);
            MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddVirtualParticlesToRhs(int finest_level, const Vector<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddVirtualParticlesToRhsML()");
    if (finest_level < parent->finestLevel())
    {
        // Should only need ghost cells for virtual particles if they're near
        // the simulation boundary and even then only maybe
        MultiFab VirtPartMF(grids[finest_level], dmap[finest_level], 1, 1);
        VirtPartMF.setVal(0.0);

        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            Nyx::theVirtualParticles()[i]->AssignDensitySingleLevel(VirtPartMF, finest_level, 1, 1);
            MultiFab::Add(*Rhs_particles[finest_level], VirtPartMF, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddGhostParticlesToRhs (int               level,
                                 MultiFab&         Rhs)
{
    BL_PROFILE("Gravity::AddGhostParticlesToRhs()");
    if (level > 0)
    {
        // If we have ghost particles, add their density to the single level solve
        MultiFab ghost_mf(grids[level], dmap[level], 1, 1);

        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            ghost_mf.setVal(0.);
            Nyx::theGhostParticles()[i]->AssignDensitySingleLevel(ghost_mf, level, 1, -1);
            MultiFab::Add(Rhs, ghost_mf, 0, 0, 1, 0);
        }
    }
}

void
Gravity::AddGhostParticlesToRhs(int level, const Vector<MultiFab*>& Rhs_particles)
{
    BL_PROFILE("Gravity::AddGhostParticlesToRhsML()");
    if (level > 0)
    {
        // We require one ghost cell in GhostPartMF because that's how we handle
        // particles near fine-fine boundaries.  However we don't add any ghost
        // cells from GhostPartMF to the RHS.
        MultiFab GhostPartMF(grids[level], dmap[level], 1, 1);
        GhostPartMF.setVal(0.0);

        // Get the Ghost particle mass function. Note that Ghost particles should
        // only affect the coarsest level so we use a single level solve.  We pass in
        // -1 for the particle_lvl_offset because that makes the particles the size
        // of the coarse, not fine, dx.
        for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        {
            Nyx::theGhostParticles()[i]->AssignDensitySingleLevel(GhostPartMF, level, 1, -1);
            MultiFab::Add(*Rhs_particles[0], GhostPartMF, 0, 0, 1, 0);
        }
    }
}

void
Gravity::CorrectRhsUsingOffset(int level, MultiFab& Rhs)
{
    BL_PROFILE("Gravity::CorrectRhsUsingOffset()");
    if (verbose)
        amrex::Print() << " ... subtracting average density from RHS in solve ... "
                       << mass_offset << '\n';

    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
        Rhs[mfi].plus(-mass_offset);
    // This checks if mass has been conserved--in particular if
    // virtual particles are correctly representing finer particles.
    if (level == 0)
    {
        Nyx* nyx_level = dynamic_cast<Nyx*>(&(parent->getLevel(0)));
        Real sum = nyx_level->vol_weight_sum(Rhs,false);

        const Real eps = 1.e-10 * std::abs(mass_offset);
        if (std::abs(sum) > eps)
        {
            amrex::Print() << " ... current avg differs from mass_offset by " << sum << " " << '\n';
            amrex::Print() << " ... Gravity : single_level solve for phi -- total mass has changed!"
                          << '\n';;
        }

        if (verbose)
            amrex::Print() << " ... subtracting " << sum << " to ensure solvability " << '\n';

        Rhs.plus(-sum, 0, 1, 0);
    }
}

#ifdef USEHPGMG
void
Gravity::solve_with_HPGMG(int level,
                          MultiFab& soln,
                          const Vector<MultiFab*>& grad_phi,
                          MultiFab& rhs,
                          Real tol,
                          Real abs_tol)
{
  BL_PROFILE("Gravity::solve_with_HPGMG()");
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Box& domain = geom.Domain();

  BndryData bd(grids[level], dmap[level], 1, geom);
  set_boundary(bd, rhs, dx);
  const int n_cell = domain.length(0);

  // We'll get the max grid (box) size from the soln MultiFab already provided
  // to us. Just pick the first valid Box we can find and measure its size.
  // HPGMG requires the Boxes to be cubes, so if they're not then a sanity
  // check in MultiFab::CreateHPGMGLevel() will catch it and quit.
  MFIter mfi(soln);
  while (!mfi.isValid()) ++mfi;
  const Box& bx = mfi.validbox();
  const int max_grid_size = bx.length(0);

  const int my_rank = ParallelDescriptor::MyProc();
  const int num_ranks = ParallelDescriptor::NProcs();

  Laplacian lap_operator(bd, dx[0]);

  const Real a = 0.0; // coefficient in front of alpha in the Helmholtz operator
  // The canonical Helmholtz operator is a alpha u - b div (beta grad(u)) = f.
  // The self-gravity equation that we solve in Nyx is Lap(u) = f. So we need
  // either the betas or b to be -1. The other GMG solvers in amrex set the
  // b*beta to -1.
  const Real b = -1.0; // coefficient in front of beta in the Helmholtz operator

  const auto& ba = rhs.boxArray();
  const auto& dm = rhs.DistributionMap();
  MultiFab alpha(ba, dm, 1, 0);
  MultiFab beta_cc(ba, dm, 1, 1);
  alpha.setVal(0.0);
  beta_cc.setVal(1.0);

  const int domain_boundary_condition = BC_PERIODIC;

  const int minCoarseDim = 2; // avoid problems with black box calculation of D^{-1} for poisson with periodic BC's on a 1^3 grid

  level_type level_h;
  mg_type MG_h;

  int numVectors = 12;

  CreateHPGMGLevel(&level_h, rhs, n_cell, max_grid_size, my_rank, num_ranks, domain_boundary_condition, numVectors, *dx);
  SetupHPGMGCoefficients(a, b, alpha, beta_cc, &level_h);
  ConvertToHPGMGLevel(rhs, n_cell, max_grid_size, &level_h, VECTOR_F);
  ConvertToHPGMGLevel(soln, n_cell, max_grid_size, &level_h, VECTOR_U);

  if (level_h.boundary_condition.type == BC_PERIODIC)
  {
    Real average_value_of_f = mean (&level_h, VECTOR_F);
    if (average_value_of_f != 0.0)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        std::cerr << "WARNING: Periodic boundary conditions, but f does not sum to zero... mean(f)=" 
                << average_value_of_f << std::endl;
        std::cerr << "Subtracting mean(f) from RHS ..." << std::endl;
      }
      shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rebuild_operator(&level_h,NULL,a,b);    // i.e. calculate Dinv and lambda_max
  MGBuild(&MG_h,&level_h,a,b,minCoarseDim,ParallelDescriptor::Communicator()); // build the Multigrid Hierarchy
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  amrex::Print() << std::endl << std::endl << "===== STARTING SOLVE =====" << std::endl << std::flush;

  MGResetTimers (&MG_h);
  //zero_vector (MG_h.levels[0], VECTOR_U);
  FMGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, tol);

  // Now convert solution from HPGMG back to rhs MultiFab.
  ConvertFromHPGMGLevel(soln, &level_h, VECTOR_U);

  MGDestroy(&MG_h); // destroys all but the finest grid
  destroy_level(&level_h); // destroys the finest grid

  lap_operator.compFlux(*grad_phi[0],*grad_phi[1],*grad_phi[2],soln);

  // We have to multiply by -1 here because the compFlux routine returns grad(phi)
  grad_phi[0]->mult(-1.0);
  grad_phi[1]->mult(-1.0);
  grad_phi[2]->mult(-1.0);
}
#endif

void
Gravity::solve_for_phi_with_mlmg (int level, MultiFab& Rhs, MultiFab& phi,
                                  const Vector<MultiFab*>& grad_phi, Real time)
{
    BL_PROFILE("Gravity::solve_for_phi_with_mlmg");
    const MultiFab* crse_bcdata = nullptr;
    MultiFab CPhi;
    if (level > 0) {
        get_crse_phi(level, CPhi, time);
        crse_bcdata = &CPhi;
    }
    Real rel_eps = sl_tol;
    Real abs_eps = 0.;
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phi_aa;
    grad_phi_aa.push_back({AMREX_D_DECL(grad_phi[0], grad_phi[1], grad_phi[2])});
    level_solver_resnorm[level] =
        solve_with_MLMG(level, level, {&phi}, {&Rhs}, grad_phi_aa, crse_bcdata, rel_eps, abs_eps);
}

void
Gravity::solve_for_delta_phi_with_mlmg (int crse_level, int fine_level, MultiFab& CrseRhs,
                                        const Vector<MultiFab*>& delta_phi,
                                        const Vector<Vector<MultiFab*> >& grad_delta_phi)
{
    BL_PROFILE("Gravity::solve_for_delta_phi_with_mlmg");

    if (verbose)
    {
        amrex::Print() << "... solving for delta_phi at crse_level = " << crse_level << '\n';
        amrex::Print() << "...                    up to fine_level = " << fine_level << '\n';
    }
    
    const int num_levels = fine_level - crse_level + 1;

    BL_ASSERT(grad_delta_phi.size() == num_levels);
    BL_ASSERT(delta_phi.size() == num_levels);

    Vector<MultiFab> rhs(num_levels);
    Vector<const MultiFab*> rhsp(num_levels);

    for (int lev = 0; lev < num_levels; ++lev) {
        delta_phi[lev]->setVal(0.0);
        if (lev == 0) {
            rhsp[lev] = &CrseRhs;
        } else {
            rhs[lev].define(grids[lev+crse_level], dmap[lev+crse_level], 1, 0);
            rhs[lev].setVal(0.0);
            rhsp[lev] = &rhs[lev];
        }
    }

    Real rel_eps = delta_tol;
    // fine_level is not included.
    Real abs_eps = *(std::max_element(level_solver_resnorm.begin() + crse_level,
                                      level_solver_resnorm.begin() + fine_level));
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad;
    for (const auto& x : grad_delta_phi) {
        grad.push_back({AMREX_D_DECL(x[0],x[1],x[2])});
    }
    solve_with_MLMG(crse_level, fine_level, delta_phi, rhsp, grad, nullptr, rel_eps, abs_eps);
}

Real
Gravity::solve_with_MLMG (int crse_level, int fine_level,
                          const Vector<MultiFab*>& phi,
                          const Vector<const MultiFab*>& rhs,
                          const Vector<std::array<MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                          const MultiFab* const crse_bcdata,
                          Real rel_eps, Real abs_eps)
{
    BL_PROFILE("Gravity::solve_with_MLMG");

    const int nlevs = fine_level - crse_level + 1;

    Vector<Geometry> gmv;
    Vector<BoxArray> bav;
    Vector<DistributionMapping> dmv;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        gmv.push_back(parent->Geom(ilev+crse_level));
        bav.push_back(rhs[ilev]->boxArray());
        dmv.push_back(rhs[ilev]->DistributionMap());
    }

    LPInfo info;
    info.setAgglomeration(mlmg_agglomeration);
    info.setConsolidation(mlmg_consolidation);

    MLPoisson mlpoisson(gmv, bav, dmv, info);

    // BC
    mlpoisson.setDomainBC(mlmg_lobc, mlmg_hibc);

    if (mlpoisson.needsCoarseDataForBC())
    {
        mlpoisson.setCoarseFineBC(crse_bcdata, parent->refRatio(crse_level-1)[0]);
    }
    
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        mlpoisson.setLevelBC(ilev, phi[ilev]);
    }

    MLMG mlmg(mlpoisson);
    mlmg.setVerbose(verbose);
    if (crse_level == 0) {
	mlmg.setMaxFmgIter(mlmg_max_fmg_iter);
    } else {
	mlmg.setMaxFmgIter(0); // Vcycle
    }

    Real final_resnorm = mlmg.solve(phi, rhs, rel_eps, abs_eps);

    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad_phi_tmp;
    for (const auto& x: grad_phi) {
        grad_phi_tmp.push_back({AMREX_D_DECL(x[0],x[1],x[2])});
    }
    mlmg.getGradSolution(grad_phi_tmp);

    return final_resnorm;
}

void
Gravity::set_boundary(BndryData& bd, MultiFab& rhs, const Real* dx)
{
  BL_PROFILE("Gravity::set_boundary()");
  for (int n=0; n<BL_SPACEDIM; ++n) {
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index();

      // Our default will be that the face of this grid is either touching another grid
      //  across an interior boundary or a periodic boundary.
      {
        // Define the type of boundary conditions to be Dirichlet (even for periodic)
        bd.setBoundCond(Orientation(n, Orientation::low) ,i,0,LO_DIRICHLET);
        bd.setBoundCond(Orientation(n, Orientation::high),i,0,LO_DIRICHLET);

        // Set the boundary conditions to the cell centers outside the domain
        bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
        bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
      }
    }
  }
}


