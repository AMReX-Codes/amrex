#ifndef AMREX_INTERP_FACE_REGISTER_H_
#define AMREX_INTERP_FACE_REGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex {

/**
 *  \brief InterpFaceRegister is a coarse/fine boundary register for
 *  interpolation of face data at the coarse/fine boundadry.
 */

class InterpFaceRegister
{
public:

    InterpFaceRegister () = default;

    /**
     * \brief InterpFaceRegister Constructor.
     *
     * \param fba  The fine level BoxArray
     * \param fdm  The fine level DistributionMapping
     * \param fgeom  The fine level Geometry
     * \param ref_ratio  The refinement ratio
     */
    InterpFaceRegister (BoxArray const& fba, DistributionMapping const& fdm,
                        Geometry const& fgeom, IntVect const& ref_ratio);

    /**
     * \brief Defines an InterpFaceRegister objecct.
     *
     * \param fba  The fine level BoxArray
     * \param fdm  The fine level DistributionMapping
     * \param fgeom  The fine level Geometry
     * \param ref_ratio  The refinement ratio
     */
    void define (BoxArray const& fba, DistributionMapping const& fdm,
                 Geometry const& fgeom, IntVect const& ref_ratio);

    /**
     * \brief Returns a coarse/fine boundary mask for a given face.
     *
     * This returns an iMultiFab at the coarse/fine boundary of a given
     * face.  The data are only defined on one face of each box in the fine
     * level BoxArray passed to the InterpFaceRegister constructor.  It has
     * two possible values: 1 for coarse/fine boundary and 0 for fine/fine
     * boundary including non-periodic physical domain face.
     *
     * \param face  The face
     */
    iMultiFab const& mask (Orientation face) const;

    /**
     * \brief Interpolates from coarse to fine data at coarse/fine
     * boundaries.
     *
     * \param fine  Array of pointers to the fine data.
     * \param crse  Array of const pointers to the coarse data.
     * \param scomp  Starting component
     * \param ncomp  Number of components
     */
    void interp (Array<MultiFab*, AMREX_SPACEDIM> const& fine,
                 Array<MultiFab const*, AMREX_SPACEDIM> const& crse,
                 int scomp, int ncomp);

private:
    BoxArray m_fine_ba;
    DistributionMapping m_fine_dm;
    Geometry m_fine_geom;
    IntVect m_ref_ratio;

    Geometry m_crse_geom;

    Array<BoxArray, 2*AMREX_SPACEDIM> m_fine_face_ba;
    Array<BoxArray, 2*AMREX_SPACEDIM> m_crse_face_ba;
    Array<iMultiFab,2*AMREX_SPACEDIM> m_face_mask; // crse/fine: 1, fine/fine & fine/physbc: 0
};

}

#endif
