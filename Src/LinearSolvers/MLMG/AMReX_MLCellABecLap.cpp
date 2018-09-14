
#include <AMReX_MLCellABecLap.H>

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

}
