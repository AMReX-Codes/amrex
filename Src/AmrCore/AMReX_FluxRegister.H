#ifndef AMREX_FLUXREGISTER_H_
#define AMREX_FLUXREGISTER_H_

#include <AMReX_BndryRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_Array.H>

namespace amrex {


/**
* \brief Flux Register.
*
* Stores and manipulates fluxes at coarse-fine interfaces.
*/

class FluxRegister
    :
    public BndryRegister
{
public:

    /**
    * \brief The default constructor.
    */
    FluxRegister();

    /**
    * \brief The constructor. This version allows setting the DistributionMapping.
    *
    * \param fine_boxes
    * \param dm
    * \param ref_ratio
    * \param fine_lev
    * \param nvar
    */
    FluxRegister (const BoxArray&            fine_boxes,
                  const DistributionMapping& dm,
                  const IntVect&             ref_ratio,
                  int                        fine_lev,
                  int                        nvar);

    /**
    * \brief The destructor.
    */
    virtual ~FluxRegister ();

    FluxRegister (FluxRegister&& rhs) noexcept = default;

    FluxRegister (const FluxRegister& rhs) = delete;
    FluxRegister& operator= (const FluxRegister& rhs) = delete;
    FluxRegister& operator= (FluxRegister&& rhs) = delete;

    //! An enum that says whether to add or copy src data to members.
    enum FrOp {COPY = 0, ADD = 1};

    /**
    * \brief Initialize after using default constructor.
    * This version allows setting the DistributionMapping.
    *
    * \param fine_boxes
    * \param dm
    * \param ref_ratio
    * \param fine_lev
    * \param nvar
    */
    void define (const BoxArray&            fine_boxes,
                 const DistributionMapping& dm,
                 const IntVect&             ref_ratio,
                 int                        fine_lev,
                 int                        nvar);

    void clear ();


    /**
    * \brief Returns the refinement ratio.
    */
    const IntVect& refRatio () const noexcept;

    /**
    * \brief Returns the level number of the fine level.
    */
    int fineLevel () const noexcept;

    /**
    * \brief Returns the level number of the coarse level (fineLevel()-1).
    */
    int crseLevel () const noexcept;

    /**
    * \brief The number of components.
    */
    int nComp () const noexcept;

    /**
    * \brief The coarsened boxes.
    */
    const BoxArray& coarsenedBoxes () const noexcept;

    /**
    * \brief Returns the sum of the registers.
    *
    * \param comp
    */
    Real SumReg (int comp) const;

    /**
    * \brief Initialize flux correction with coarse data.
    *
    * \param mflx
    * \param area
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    * \param op
    */
    void CrseInit (const MultiFab& mflx,
                   const MultiFab& area,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0,
                   FrOp            op = FluxRegister::COPY);

    /**
    * \brief Initialize flux correction with coarse data.
    *
    * \param mflx
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    * \param op
    */
    void CrseInit (const MultiFab& mflx,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0,
                   FrOp            op = FluxRegister::COPY);

    /**
    * \brief Add coarse fluxes to the flux register.
    * This is different from CrseInit with FluxRegister::ADD.
    * This is used for cases in which the grids covered by fine
    * do not have fluxes computed (e.g., FLASH).
    *
    * \param mflx
    * \param area
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    * \param geom
    */
    void CrseAdd (const MultiFab& mflx,
                  const MultiFab& area,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult,
                  const Geometry& geom);

    /**
    * \brief /in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
    *
    * \param mflx
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    * \param geom
    */
    void CrseAdd (const MultiFab& mflx,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult,
                  const Geometry& geom);

    /**
    * \brief Increment flux correction with fine data.
    *
    * /in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
    *
    * \param mflx
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    */
    void FineAdd (const MultiFab& mflx,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);

    /**
    * \brief Increment flux correction with fine data.
    *
    * \param mflx
    * \param area
    * \param dir
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    */
    void FineAdd (const MultiFab& mflx,
                  const MultiFab& area,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);

    /**
    * \brief Increment flux correction with fine data.
    *
    * in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
    *
    * \param flux
    * \param dir
    * \param boxno
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    */
    void FineAdd (const FArrayBox& flux,
                  int              dir,
                  int              boxno,
                  int              srccomp,
                  int              destcomp,
                  int              numcomp,
                  Real             mult) noexcept;

    /**
    * \brief Increment flux correction with fine data.
    *
    * \param flux
    * \param area
    * \param dir
    * \param boxno
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param mult
    */
    void FineAdd (const FArrayBox& flux,
                  const FArrayBox& area,
                  int              dir,
                  int              boxno,
                  int              srccomp,
                  int              destcomp,
                  int              numcomp,
                  Real             mult) noexcept;

    /**
    * \brief Set flux correction data for a fine box (given by boxno) to a given value.
    * This routine used by FLASH does NOT run on gpu for safety.
    *
    * \param dir
    * \param boxno
    * \param destcomp
    * \param numcomp
    * \param val
    */
    void FineSetVal (int              dir,
                     int              boxno,
                     int              destcomp,
                     int              numcomp,
                     Real             val) noexcept;

    /**
    * \brief Apply flux correction.  Note that this takes the coarse Geometry.
    *
    * \param mf
    * \param volume
    * \param scale
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param crse_geom
    */
    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 int             dir,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

    /**
    * \brief Constant volume version of Reflux().  Note that this takes the coarse Geometry.
    *
    * \param mf
    * \param scale
    * \param srccomp
    * \param destcomp
    * \param numcomp
    * \param crse_geom
    */
    void Reflux (MultiFab&       mf,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

    void Reflux (MultiFab&       mf,
                 int             dir,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

    void OverwriteFlux (Array<MultiFab*,AMREX_SPACEDIM> const& crse_fluxes,
                        Real scale, int srccomp, int destcomp, int numcomp,
                        const Geometry& crse_geom);


    /**
    * \brief Set internal borders to zero
    *
    * \param crse_geom
    */
    void ClearInternalBorders (const Geometry& crse_geom);

    /**
    * \brief Write (used for writing to checkpoint)
    *
    * \param name
    * \param os
    */
    void write (const std::string& name, std::ostream& os) const;

    /**
    * \brief Read (used for reading from checkpoint)
    *
    * \param name
    * \param is
    */
    void read (const std::string& name, std::istream& is);

// public for cuda
public:
    void Reflux (MultiFab& mf, const MultiFab& volume, Orientation face,
                 Real scale, int scomp, int dcomp, int nc, const Geometry& geom);

private:

    //! Refinement ratio
    IntVect ratio;

    //! Current level + 1.
    int fine_level;

    //! Number of state components.
    int ncomp;
};

}

#endif /*_FLUXREGISTER_H_*/
