#ifndef AMREX_INTERPOLATER_H_
#define AMREX_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_InterpBase.H>

namespace amrex {

class Geometry;
class FArrayBox;
class IArrayBox;

/**
* \brief Virtual base class for interpolaters.
*
* Specifies interpolater interface for coarse-to-fine interpolation in space.
*/

class Interpolater
    : public InterpBase
{
public:
    /**
    * \brief Coarse to fine interpolation in space.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const & bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) = 0;

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param fine_values
    * \param fine_known
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    virtual void interp_face (const FArrayBox& /*crse*/,
                              const int        /*crse_comp*/,
                              FArrayBox&       /*fine*/,
                              const int        /*fine_comp*/,
                              const int        /*ncomp*/,
                              const Box&       /*fine_region*/,
                              const IntVect&   /*ratio*/,
                              const IArrayBox& /*solve_mask*/,
                              const Geometry&  /*crse_geom*/,
                              const Geometry&  /*fine_geom*/,
                              Vector<BCRec> const & /*bcr*/,
                              const int        /*bccomp*/,
                              RunOn            /*runon*/)
    { amrex::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead."); }

    /**
    * \brief Coarse to fine interpolation in space.
    *        Takes an Array of FArrayBox*s
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param fine_values
    * \param fine_known
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& /*crse*/,
                             const int         /*crse_comp*/,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& /*fine*/,
                             const int         /*fine_comp*/,
                             const int         /*ncomp*/,
                             const Box&        /*fine_region*/,
                             const IntVect&    /*ratio*/,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& /*solve_mask*/,
                             const Geometry&   /*crse_geom*/,
                             const Geometry&   /*fine_geom*/,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                             const int         /*actual_comp*/,
                             const int         /*actual_state*/,
                             const RunOn       /*runon*/)

    { amrex::Abort("The Array<FArrayBox*, AMREX_SPACEDIM> version of this Interpolater is not implemented or does not apply. Call 'interp' instead."); }

    /**
    * \brief Re-visit the interpolation to protect against under- or overshoots.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param fine_state
    * \param state_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    */
    virtual void protect (const FArrayBox& /*crse*/,
                          int              /*crse_comp*/,
                          FArrayBox&       /*fine*/,
                          int              /*fine_comp*/,
                          FArrayBox&       /*fine_state*/,
                          int              /*state_comp*/,
                          int              /*ncomp*/,
                          const Box&       /*fine_region*/,
                          const IntVect&   /*ratio*/,
                          const Geometry&  /*crse_geom*/,
                          const Geometry&  /*fine_geom*/,
                          Vector<BCRec>&   /*bcr*/,
                          RunOn            /*runon*/) {}
};


/**
* \brief Bilinear interpolation on node centered data.
*
* Bilinear interpolation on node centered data.
*/

class NodeBilinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};


/**
* \brief Bilinear interpolation on cell centered data.
*
* Bilinear interpolation on cell centered data.
*/

class CellBilinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};


/**
* \brief Linear conservative interpolation on cell centered data
*
*
* Linear conservative interpolation on cell centered data
* I.e, conservative interpolation with a limiting scheme that
* preserves the value of any linear combination  of the
* fab components; //!< e.g.,
*
* if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
* sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
* in all fine cells if,jf covering coarse cell ic,jc.
*/

class CellConservativeLinear
    :
    public Interpolater
{
public:

    /**
    * \brief The constructor.
    *
    * \param do_linear_limiting_
    */
    explicit CellConservativeLinear (bool do_linear_limiting_= true);

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;
    /**
    * \brief Coarse to fine interpolation in space.
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              /*actual_comp*/,
                         int              /*actual_state*/,
                         RunOn            runon) override;

protected:

    bool do_linear_limiting;
};


/**
* \brief Lin. cons. interp. on cc data with protection against under/over-shoots.
*
* Linear conservative interpolation on cell centered data
* but with protection against undershoots or overshoots.
*/

class CellConservativeProtected
    :
    public CellConservativeLinear
{
public:

    /**
    * \brief The constructor.
    */
    CellConservativeProtected ();

    /**
    * \brief Re-visit the interpolation to protect against under- or overshoots.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param fine_state
    * \param state_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    */
    void protect (const FArrayBox& crse,
                          int              crse_comp,
                          FArrayBox&       fine,
                          int              fine_comp,
                          FArrayBox&       fine_state,
                          int              state_comp,
                          int              ncomp,
                          const Box&       fine_region,
                          const IntVect&   ratio,
                          const Geometry&  crse_geom,
                          const Geometry&  fine_geom,
                          Vector<BCRec>&   bcr,
                          RunOn            runon) override;
};


/**
* \brief Quadratic interpolation on cell centered data.
*
* Quadratic interpolation on cell centered data.
*/

class CellQuadratic
    :
    public Interpolater
{
public:

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};


/**
* \brief Piecewise Constant interpolation on cell centered data.
*/

class PCInterp
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const&  bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};


/**
* \brief Conservative quartic interpolation on cell averaged data.
*
* An order 4 polynomial is used to fit the data.  For each cell involved
* in constructing the polynomial, the average of the polynomial inside that
* cell is equal to the cell averaged value of the original data.
*/

class CellConservativeQuartic
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const&  bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};

/**
* \brief Divergence-free interpolation on face centered data.
*
* Divergence-free interpolation on face centered data.
*/

class FaceDivFree
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    * \param runon
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param solve_mask
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    * \param runon
    */
    void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                             int         crse_comp,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                             int         fine_comp,
                             int         ncomp,
                             const Box&        fine_region,
                             const IntVect&    ratio,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                             const Geometry&   crse_geom,
                             const Geometry&   fine_geom,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& bcr,
                             int         actual_comp,
                             int         actual_state,
                             RunOn       runon) override;
};


/**
* \brief Bilinear interpolation on face data.
*
* Bilinear interpolation on data.
*/

class FaceLinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param solve_mask
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    void interp_face (const FArrayBox& crse,
                              int        crse_comp,
                              FArrayBox&       fine,
                              int        fine_comp,
                              int        ncomp,
                              const Box&       fine_region,
                              const IntVect&   ratio,
                              const IArrayBox& solve_mask,
                              const Geometry&  crse_geom,
                              const Geometry&  fine_geom,
                              Vector<BCRec> const & bcr,
                              int        bccomp,
                              RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                             int         crse_comp,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                             int         fine_comp,
                             int         ncomp,
                             const Box&        fine_region,
                             const IntVect&    ratio,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                             const Geometry&   /*crse_geom*/,
                             const Geometry&   /*fine_geom*/,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                             int         /*actual_comp*/,
                             int         /*actual_state*/,
                             RunOn       runon) override;


};

/**
* \brief Quartic interpolation on cell centered data.
*
* Quartic interpolation on cell centered data.
*/

class CellQuartic
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};

//! CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
extern AMREX_EXPORT PCInterp                  pc_interp;
extern AMREX_EXPORT NodeBilinear              node_bilinear_interp;
extern AMREX_EXPORT FaceDivFree               face_divfree_interp;
extern AMREX_EXPORT FaceLinear                face_linear_interp;
extern AMREX_EXPORT CellConservativeLinear    lincc_interp;
extern AMREX_EXPORT CellConservativeLinear    cell_cons_interp;
extern AMREX_EXPORT CellBilinear              cell_bilinear_interp;
extern AMREX_EXPORT CellConservativeProtected protected_interp;
extern AMREX_EXPORT CellConservativeQuartic   quartic_interp;
extern AMREX_EXPORT CellQuadratic             quadratic_interp;
extern AMREX_EXPORT CellQuartic               cell_quartic_interp;

}

#endif /*_INTERPOLATER_H_*/
