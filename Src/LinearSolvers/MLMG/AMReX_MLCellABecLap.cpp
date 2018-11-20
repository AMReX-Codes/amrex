
#include <AMReX_MLCellABecLap.H>

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

namespace amrex {

MLCellABecLap::MLCellABecLap ()
{
}

MLCellABecLap::~MLCellABecLap () {}

void
MLCellABecLap::define (const Vector<Geometry>& a_geom,
                       const Vector<BoxArray>& a_grids,
                       const Vector<DistributionMapping>& a_dmap,
                       const LPInfo& a_info,
                       const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLCellABecLap::prepareForSolve ()
{
    MLCellLinOp::prepareForSolve();
}

void
MLCellABecLap::update ()
{
    if (MLCellLinOp::needsUpdate()) MLCellLinOp::update();
}

void
MLCellABecLap::getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                          const Vector<MultiFab*>& a_sol,
                          Location a_loc) const
{
    BL_PROFILE("MLMG::getFluxes()");

    const Real betainv = 1.0 / getBScalar();
    const int nlevs = NAMRLevels();
    for (int alev = 0; alev < nlevs; ++alev) {
        compFlux(alev, a_flux[alev], *a_sol[alev], a_loc);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            unapplyMetricTerm(alev, 0, *a_flux[alev][idim]);
            if (betainv != 1.0) {
                a_flux[alev][idim]->mult(betainv);
            }
        }
    }
}

#ifdef AMREX_USE_HYPRE
std::unique_ptr<Hypre>
MLCellABecLap::makeHypre (Hypre::Interface hypre_interface) const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    MPI_Comm comm = BottomCommunicator();

    auto hypre_solver = amrex::makeHypre(ba, dm, geom, comm, hypre_interface);

    hypre_solver->setScalars(getAScalar(), getBScalar());

    const int mglev = NMGLevels(0)-1;
    auto ac = getACoeffs(0, mglev);
    if (ac)
    {
        hypre_solver->setACoeffs(*ac);
    }
    else
    {
        MultiFab alpha(ba,dm,1,0,MFInfo(),factory);
        alpha.setVal(0.0);
        hypre_solver->setACoeffs(alpha);
    }

    auto bc = getBCoeffs(0, mglev);
    if (bc[0])
    {
        hypre_solver->setBCoeffs(bc);
    }
    else
    {
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            beta[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                              dm, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }
        hypre_solver->setBCoeffs(amrex::GetArrOfConstPtrs(beta));
    }

    return hypre_solver;
}
#endif

#ifdef AMREX_USE_PETSC
std::unique_ptr<PETScABecLap>
MLCellABecLap::makePETSc () const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    MPI_Comm comm = BottomCommunicator();
    
    auto petsc_solver = makePetsc(ba, dm, geom, comm);

    petsc_solver->setScalars(getAScalar(), getBScalar());

    const int mglev = NMGLevels(0)-1;
    auto ac = getACoeffs(0, mglev);
    if (ac)
    {
        petsc_solver->setACoeffs(*ac);
    }
    else
    {
        MultiFab alpha(ba,dm,1,0,MFInfo(),factory);
        alpha.setVal(0.0);
        petsc_solver->setACoeffs(alpha);
    }

    auto bc = getBCoeffs(0, mglev);
    if (bc[0])
    {
        petsc_solver->setBCoeffs(bc);
    }
    else
    {
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            beta[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                              dm, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }
        petsc_solver->setBCoeffs(amrex::GetArrOfConstPtrs(beta));
    }
    return petsc_solver;
}
#endif

}
