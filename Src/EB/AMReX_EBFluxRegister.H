#ifndef AMREX_EBFLUXREGISTER_H_
#define AMREX_EBFLUXREGISTER_H_

#include <AMReX_YAFluxRegister.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

/**
  EBFluxRegister is used for refluxing, re-redistribution,
  re-refluxing, and re-re-redistribution.  See `Tutorials/EB/CNS` for
  an example of using this class.

  Reflux for EB is somewhat similar to non-EB.  At the beginning of a
  coarse step, `reset()` is called.  In MFIter for the coarse level
  advance, `CrseAdd` is called with coarse flux.  There are two
  versions of `CrseAdd`, one for regular fab/tile and the other for
  fab/tile containing cutcells.  For the cutcell version, the flux is
  is at centroid, not face center.  The flux is not scaled.  In MFIter
  for the fine level advance, `FineAdd` is called.  There are also two
  version of `FineAdd`, one for regular and the other for fab/tile
  containing cutcells.  The cutcell version also takes `dm` for
  re-redistribution explained below.  After the fine level finished
  its time steps, `Reflux` is called to update the coarse cells next
  to the coarse/fine boundary.  Note that re-redistribution is also
  performed in `Reflux`.

  Re-redistribution is unfortunately more complicated.  The coarse
  level needs to accumulate the *density* (e.g., g/cm^3 for mass
  density) *loss* of coarse/fine boundary cells (i.e., coarse cells on
  the coarse side of the coarse/fine boundary) due to redistribution
  from/to coarse cells covered by the fine level.  If a coarse/fine
  boundary cell redistributes mass to a covered cell, this counts as a
  positive loss.  Here covered means covered be the fine level, not
  EB.  If a covered cell redistributes mass to a boundary cell, this
  counts as a negative loss for that boundary cell.  The accumulation
  as the coarse level is performed on `FArrayBox` returned by
  `EBFluxRegister::getCrseData(MFIter const&)`.  EBFluxRegister also
  has `getCrseFlag(MFIte cons&)` function returning `IArrayBox` that
  contains flags for genuine coarse cell, coarse/fine boundary cells
  and covered cells.

  The fine level in re-redistribution needs to accumulate the *mass*
  (e.g., g for density) *gain* of ghost cells due to redistribution
  from/to valid cells.  The application code can use a local
  `FArrayBox` to store this and then EBFLuxRegister::FineAdd is called
  to add the part in ghost cells (excluding ghost cells covered by
  valid cells of other grids) to EBFluxRegister's internal data.
*/

class EBFluxRegister
    : public YAFluxRegister
{
public:

    EBFluxRegister () = default;

    EBFluxRegister (const BoxArray& fba, const BoxArray& cba,
                    const DistributionMapping& fdm, const DistributionMapping& cdm,
                    const Geometry& fgeom, const Geometry& cgeom,
                    const IntVect& ref_ratio, int fine_lev, int nvar);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 const IntVect& ref_ratio, int fine_lev, int nvar);

    using YAFluxRegister::CrseAdd;
    void CrseAdd (const MFIter& mfi,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt,
                  const FArrayBox& volfrac,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac);

    using YAFluxRegister::FineAdd;
    void FineAdd (const MFIter& mfi,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt,
                  const FArrayBox& volfrac,
                  const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac,
                  const FArrayBox& dm);

    void Reflux (MultiFab& crse_state, const amrex::MultiFab& crse_vfrac,
                 MultiFab& fine_state, const amrex::MultiFab& fine_vfrac);

    FArrayBox* getCrseData (const MFIter& mfi) {
        return &(m_crse_data[mfi]);
    }

    const IArrayBox* getCrseFlag (const MFIter& mfi) const {
        return &(m_crse_flag[mfi]);
    }

private:

    iMultiFab m_cfp_inside_mask;

    void defineExtra (const BoxArray& fba, const DistributionMapping& fdm);
};

}

#endif
