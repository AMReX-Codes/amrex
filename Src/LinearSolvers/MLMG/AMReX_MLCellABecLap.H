#ifndef AMREX_ML_CELL_ABECLAP_H_
#define AMREX_ML_CELL_ABECLAP_H_

#include <AMReX_MLCellLinOp.H>

namespace amrex {

class MLCellABecLap
    : public MLCellLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    MLCellABecLap ();
    virtual ~MLCellABecLap ();

    MLCellABecLap (const MLCellABecLap&) = delete;
    MLCellABecLap (MLCellABecLap&&) = delete;
    MLCellABecLap& operator= (const MLCellABecLap&) = delete;
    MLCellABecLap& operator= (MLCellABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual bool needsUpdate () const override {
        return MLCellLinOp::needsUpdate();
    }
    virtual void update () override;

    virtual void prepareForSolve () override;

    virtual void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                            const Vector<MultiFab*>& a_sol,
                            Location a_loc) const final override;

    virtual Real getAScalar () const = 0;
    virtual Real getBScalar () const = 0;
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const = 0;
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const = 0;

#ifdef AMREX_USE_HYPRE
    virtual std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override;
#endif

#ifdef AMREX_USE_PETSC
    virtual std::unique_ptr<PETScABecLap> makePETSc () const override;
#endif
};

}

#endif
