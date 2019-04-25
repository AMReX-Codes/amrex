#ifndef AMREX_YAFLUXREGISTER_H_
#define AMREX_YAFLUXREGISTER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <array>

namespace amrex {

/**
  YAFluxRegister is yet another Flux Register for refluxing.

  At the beginning of a coarse step, `reset()` is called.  In MFIter
  for the coarse level advance, `CrseAdd` is called with coarse flux.
  The flux is not scaled.  In MFIter for the fine level advance,
  `FineAdd` is called.  After the fine level finished its time steps,
  `Reflux` is called to update the coarse cells next to the
  coarse/fine boundary.
*/

class YAFluxRegister
{
public:

    YAFluxRegister () = default;

    YAFluxRegister (const BoxArray& fba, const BoxArray& cba,
                    const DistributionMapping& fdm, const DistributionMapping& cdm,
                    const Geometry& fgeom, const Geometry& cgeom,
                    const IntVect& ref_ratio, int fine_lev, int nvar);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 const IntVect& ref_ratio, int fine_lev, int nvar);

    void reset ();

    void CrseAdd (const MFIter& mfi,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt) noexcept;

    void FineAdd (const MFIter& mfi,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt) noexcept;

    void Reflux (MultiFab& state, int dc = 0);

    bool CrseHasWork (const MFIter& mfi) const noexcept {
        return m_crse_fab_flag[mfi.LocalIndex()] != crse_cell;
    }

    bool FineHasWork (const MFIter& mfi) const noexcept {
        return !(m_cfp_fab[mfi.LocalIndex()].empty());
    }

    enum CellType : int {
        // must be same as in AMReX_YAFluxRegiser_K.H
        crse_cell = 0, crse_fine_boundary_cell, fine_cell
    };

protected:

    MultiFab m_crse_data;
    iMultiFab m_crse_flag;
    Vector<int> m_crse_fab_flag;

    MultiFab m_cfpatch;                   //!< This is built on crse/fine patches
    MultiFab m_cfp_mask;
    Vector<Vector<FArrayBox*> > m_cfp_fab;  //!< The size of this is (# of local fine grids (# of crse/fine patches for that grid))
    Vector<int> m_cfp_localindex;

    Geometry m_fine_geom;
    Geometry m_crse_geom;

    IntVect m_ratio;
    int m_fine_level;
    int m_ncomp;
};

}

#endif
