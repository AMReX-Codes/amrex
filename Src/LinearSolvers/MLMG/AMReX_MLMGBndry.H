#ifndef AMREX_MLMGBNDRY_H_
#define AMREX_MLMGBNDRY_H_
#include <AMReX_Config.H>

#include <AMReX_InterpBndryData.H>

namespace amrex {

/**
 * \file AMReX_MLMGBndry.H
 *
 * Boundary data helper for MLMG that wraps `InterpBndryDataT` with general
 * BC tagging utilities.
 */

template <typename MF>
/// \brief Boundary helper for MLMG that manages coarse/fine and physical BC metadata.
class MLMGBndryT
    : public InterpBndryDataT<MF>
{
public:

    using BCTuple = Array<BoundCond,2*BL_SPACEDIM>;
    using RealTuple = typename BndryDataT<MF>::RealTuple;

    /**
     * \brief Allocate BC metadata for \p _grids/\p _geom.
     *
     * \param _grids Grid array whose ghost faces require BCs.
     * \param _dmap  Distribution mapping for the grids.
     * \param _ncomp Number of components handled.
     * \param _geom  Geometry describing the domain.
     */
    MLMGBndryT (const BoxArray& _grids, const DistributionMapping& _dmap,
                int _ncomp, const Geometry& _geom);

    ~MLMGBndryT () = default;

    MLMGBndryT (MLMGBndryT<MF>&& rhs) = delete;
    MLMGBndryT (const MLMGBndryT<MF>& rhs) = delete;
    MLMGBndryT<MF>& operator= (const MLMGBndryT<MF>& rhs) = delete;
    MLMGBndryT<MF>& operator= (MLMGBndryT<MF>&& rhs) = delete;

    /**
     * \brief Populate BC tags at coarse/fine and physical boundaries
     * for an integer refinement ratio.
     *
     * \param lo  Per-level low-side BC types.
     * \param hi  Per-level high-side BC types.
     * \param ratio Refinement ratio between coarse/fine.
     * \param a_loc Distance from the coarse/fine boundary face to the
     *              coarse/fine boundary data.
     * \param a_crse_fine_bc_type Type imposed at coarse/fine interfaces.
     */
    void setLOBndryConds (const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& lo,
                          const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& hi,
                          int ratio, const RealVect& a_loc,
                          LinOpBCType a_crse_fine_bc_type = LinOpBCType::Dirichlet);

    /**
     * \brief Populate BC tags at coarse/fine and physical boundaries
     * for anisotropic refinement ratios.
     *
     * \param lo  Per-level low-side BC types.
     * \param hi  Per-level high-side BC types.
     * \param ratio Directional refinement ratios.
     * \param a_loc Distance from the coarse/fine boundary face to the
     *              coarse/fine boundary data.
     * \param a_crse_fine_bc_type Type imposed at coarse/fine interfaces.
     */
    void setLOBndryConds (const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& lo,
                          const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& hi,
                          IntVect const& ratio, const RealVect& a_loc,
                          LinOpBCType a_crse_fine_bc_type = LinOpBCType::Dirichlet);

    /**
     * \brief Helper that sets up BC tuples for a single box.
     *
     * \param bloc             Output BC offsets (distance to boundary).
     * \param bctag            Output BC tags.
     * \param bx               Target box.
     * \param domain           Domain bounding box.
     * \param lo               Low-side BC types.
     * \param hi               High-side BC types.
     * \param dx               Cell size array.
     * \param ratio            Refinement ratio relative to the interior.
     * \param interior_bloc    Distance from the coarse/fine boundary face
     *                          to the boundary data.
     * \param domain_bloc_lo   Low-side offsets for domain boundaries.
     * \param domain_bloc_hi   High-side offsets for domain boundaries.
     * \param is_periodic      Periodicity mask.
     * \param a_crse_fine_bc_type Type imposed at coarse/fine interfaces.
     */
    static void setBoxBC (RealTuple& bloc, BCTuple& bctag, const Box& bx,
                          const Box& domain,
                          const Array<LinOpBCType,AMREX_SPACEDIM>& lo,
                          const Array<LinOpBCType,AMREX_SPACEDIM>& hi,
                          const Real* dx, IntVect const& ratio,
                          const RealVect& interior_bloc,
                          const Array<Real,AMREX_SPACEDIM>& domain_bloc_lo,
                          const Array<Real,AMREX_SPACEDIM>& domain_bloc_hi,
                          const GpuArray<int,AMREX_SPACEDIM>& is_periodic,
                          LinOpBCType a_crse_fine_bc_type);
};

template <typename MF>
MLMGBndryT<MF>::MLMGBndryT (const BoxArray& _grids, const DistributionMapping& _dmap,
                            int _ncomp, const Geometry& _geom)
    : InterpBndryDataT<MF>(_grids,_dmap,_ncomp,_geom)
{}

template <typename MF>
void
MLMGBndryT<MF>::setLOBndryConds (const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& lo,
                                 const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& hi,
                                 int ratio, const RealVect& a_loc,
                                 LinOpBCType a_crse_fine_bc_type)
{
    setLOBndryConds(lo, hi, IntVect(ratio), a_loc, a_crse_fine_bc_type);
}

template <typename MF>
void
MLMGBndryT<MF>::setLOBndryConds (const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& lo,
                                 const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& hi,
                                 IntVect const& ratio, const RealVect& a_loc,
                                 LinOpBCType a_crse_fine_bc_type)
{
    const BoxArray& ba     = this->boxes();
    const Real*     dx     = this->geom.CellSize();
    const Box&      domain = this->geom.Domain();
    const GpuArray<int,AMREX_SPACEDIM>& is_periodic = this->geom.isPeriodicArray();

    AMREX_ASSERT(a_crse_fine_bc_type == LinOpBCType::Dirichlet ||
                 a_crse_fine_bc_type == LinOpBCType::Neumann);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(this->bndry[Orientation(0,Orientation::low)]);
         fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        const Box&                 grd   = ba[i];
        RealTuple&                 bloc  = this->bcloc[fsi];
        Vector< Vector<BoundCond> >& bctag = this->bcond[fsi];

        for (int icomp = 0; icomp < this->nComp(); ++icomp) {
            BCTuple bct;
            setBoxBC(bloc, bct, grd, domain, lo[icomp], hi[icomp], dx, ratio, a_loc,
                     {{AMREX_D_DECL(Real(0.),Real(0.),Real(0.))}},
                     {{AMREX_D_DECL(Real(0.),Real(0.),Real(0.))}},
                     is_periodic, a_crse_fine_bc_type);
            for (int idim = 0; idim < 2*AMREX_SPACEDIM; ++idim) {
                bctag[idim][icomp] = bct[idim];
            }
        }
    }
}

template <typename MF>
void
MLMGBndryT<MF>::setBoxBC (RealTuple& bloc, BCTuple& bctag,
                          const Box& bx, const Box& domain,
                          const Array<LinOpBCType,AMREX_SPACEDIM>& lo,
                          const Array<LinOpBCType,AMREX_SPACEDIM>& hi,
                          const Real* dx, IntVect const& ratio,
                          const RealVect& interior_bloc,
                          const Array<Real,AMREX_SPACEDIM>& domain_bloc_lo,
                          const Array<Real,AMREX_SPACEDIM>& domain_bloc_hi,
                          const GpuArray<int,AMREX_SPACEDIM>& is_periodic,
                          LinOpBCType a_crse_fine_bc_type)
{
    using T = typename MF::value_type;

    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face = fi();
        const int         dir  = face.coordDir();

        if (domain[face] == bx[face] && !is_periodic[dir])
        {
            // All physical bc values are located on face.
            bloc[face] = static_cast<T>
                (face.isLow() ? domain_bloc_lo[dir] : domain_bloc_hi[dir]);
            const auto linop_bc  = face.isLow() ? lo[dir] : hi[dir];
            if (linop_bc == LinOpBCType::Dirichlet) {
                bctag[face] = AMREX_LO_DIRICHLET;
            } else if (linop_bc == LinOpBCType::Neumann) {
                bctag[face] = AMREX_LO_NEUMANN;
            } else if (linop_bc == LinOpBCType::reflect_odd) {
                bctag[face] = AMREX_LO_REFLECT_ODD;
            } else {
                amrex::Abort("MLMGBndry::setBoxBC: Unknown LinOpBCType");
            }
        }
        else
        {
            // Internal bndry.
            bctag[face] = static_cast<int>(a_crse_fine_bc_type);
            bloc[face]  = static_cast<T>(ratio[dir] > 0 ? Real(0.5)*static_cast<Real>(ratio[dir])*dx[dir]
                                                   : interior_bloc[dir]);
            // If this is next to another same level box, BC type is wrong
            // and bloc is wrong.  But it doesn't matter, because we also
            // have mask.  It is used only if mask says it is next to coarse
            // cells.
        }
    }
}

using MLMGBndry = MLMGBndryT<MultiFab>;

}

#endif
