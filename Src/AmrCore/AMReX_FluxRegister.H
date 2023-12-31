#ifndef AMREX_FLUXREGISTER_H_
#define AMREX_FLUXREGISTER_H_
#include <AMReX_Config.H>

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
    ~FluxRegister () = default;

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
    * \brief Initialize flux correction with coarse data (DG).
    *
    * \param SurfaceFlux
    * \param iDimX
    * \param nFields
    * \param nFields
    * \param nDOFX_X1
    * \param nDOFX_X2
    * \param nDOFX_X3
    * \param WeightsX_X1
    * \param WeightsX_X2
    * \param WeightsX_X3
    * \param op
    */
    void CrseInit_DG
           ( const MultiFab & SurfaceFlux,
             int              iDimX,
             int              nFields,
             int              nDOFX_X1,
             int              nDOFX_X2,
             int              nDOFX_X3,
             Real           * WeightsX_X1,
             Real           * WeightsX_X2,
             Real           * WeightsX_X3,
             FrOp            op = FluxRegister::COPY );

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
    * \brief /in this version the area is assumed to multiplied into the flux (if not, use scale to fix)
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
    * /in this version the area is assumed to multiplied into the flux (if not, use scale to fix)
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
    * \brief Increment flux correction with fine data (DG).
    *
    * \param SurfaceFluxes
    * \param iDimX
    * \param nFields
    * \param FaceRatio
    * \param nDOFX_X1
    * \param nDOFX_X2
    * \param nDOFX_X3
    * \param nFineX_X1
    * \param nFineX_X2
    * \param nFineX_X3
    * \param WeightsX_X1
    * \param WeightsX_X2
    * \param WeightsX_X3
    */
    void FineAdd_DG
           ( const MultiFab & SurfaceFluxes,
             int              iDimX,
             int              nFields,
             Real             FaceRatio,
             int              nDOFX_X1,
             int              nDOFX_X2,
             int              nDOFX_X3,
             int              nFineX_X1,
             int              nFineX_X2,
             int              nFineX_X3,
             Real           * WeightsX_X1,
             Real           * WeightsX_X2,
             Real           * WeightsX_X3,
             void           * vpLX_X1_Refined,
             void           * vpLX_X2_Refined,
             void           * vpLX_X3_Refined );

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
    * in this version the area is assumed to multiplied into the flux (if not, use scale to fix)
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
                  Real             mult,
                  RunOn            runon) noexcept;

    /**
    * \brief Increment flux correction with fine data for a given box (DG).
    *
    * \param SurfaceFluxes
    * \param iDimX
    * \param nFields
    * \param FaceRatio
    * \param BoxNumber
    * \param BoxNumber
    * \param nDOFX_X1
    * \param nDOFX_X2
    * \param nDOFX_X3
    * \param nFineX_X1
    * \param nFineX_X2
    * \param nFineX_X3
    * \param WeightsX_X1
    * \param WeightsX_X2
    * \param WeightsX_X3
    * \param runon
    */
    void FineAdd_DG
           ( const FArrayBox    & SurfaceFluxes,
             int                  iDimX,
             int                  nFields,
             Real                 FaceRatio,
             int                  nDOFX_X1,
             int                  nDOFX_X2,
             int                  nDOFX_X3,
             Real               * WeightsX_X1,
             Real               * WeightsX_X2,
             Real               * WeightsX_X3,
             Array4<Real const>   LX_X1_Refined,
             Array4<Real const>   LX_X2_Refined,
             Array4<Real const>   LX_X3_Refined,
             int                  BoxNumber,
             RunOn                runon ) noexcept;

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
                  Real             mult,
                  RunOn            runon) noexcept;

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
                     Real             val,
                     RunOn            runon) noexcept;

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
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 int             dir,
                 Real            scale,
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    /**
    * \brief Apply flux correction (DG).
    * Note that this takes the coarse Geometry.
    * MF_G is a multifab that contains geometry fields, in particular,
    * the square root of the determinant of the spatial three-metric.
    * MF_dU is a multifab that contains the fluid increment, dU.
    * crse_geom is the geometry on the coarse level.
    * nDOFX is the number of degrees of freedom per fluid field, per element.
    * nDOFX_X1 is the number of degrees of freedom per fluid field, per element,
    * per X1 face.
    * nDOFX_X2 is the number of degrees of freedom per fluid field, per element,
    * per X2 face.
    * nDOFX_X3 is the number of degrees of freedom per fluid field, per element,
    * per X3 face.
    * nFields is the number of distinct fluid fields (density, velocity, etc.).
    * iGF_SqrtGm is the index of MF_G that contains the square root of the
    * spatial three-metric determinant.
    * NodeNumberTableX_X1 is a 1D array that maps degrees of freedom on an
    * element to degrees of freedom on an interface in the X1-direction.
    * NodeNumberTableX_X2 is a 1D array that maps degrees of freedom on an
    * element to degrees of freedom on an interface in the X2-direction.
    * NodeNumberTableX_X3 is a 1D array that maps degrees of freedom on an
    * element to degrees of freedom on an interface in the X3-direction.
    * WeightsX_q is a 1D array that contains the Gauss--Legendre quadrature
    * weights for an element
    * LX_X1_Up is a 1D array containing the Lagrange polynomials interpolated
    * to the upper side of an element in the X1-direction.
    * LX_X1_Dn is a 1D array containing the Lagrange polynomials interpolated
    * to the lower side of an element in the X1-direction.
    * LX_X2_Up is a 1D array containing the Lagrange polynomials interpolated
    * to the upper side of an element in the X2-direction.
    * LX_X2_Dn is a 1D array containing the Lagrange polynomials interpolated
    * to the lower side of an element in the X2-direction.
    * LX_X3_Up is a 1D array containing the Lagrange polynomials interpolated
    * to the upper side of an element in the X3-direction.
    * LX_X3_Dn is a 1D array containing the Lagrange polynomials interpolated
    * to the lower side of an element in the X3-direction.
    * dX1 is the mesh width of the coarse level in the X1 direction
    * dX2 is the mesh width of the coarse level in the X2 direction
    * dX3 is the mesh width of the coarse level in the X3 direction
    *
    * \param MF_G
    * \param MF_dU
    * \param crse_geom
    * \param nDOFX
    * \param nDOFX_X1
    * \param nDOFX_X2
    * \param nDOFX_X3
    * \param nFields
    * \param iGF_SqrtGm
    * \param NodeNumberTableX_X1
    * \param NodeNumberTableX_X2
    * \param NodeNumberTableX_X3
    * \param WeightsX_q
    * \param LX_X1_Up
    * \param LX_X1_Dn
    * \param LX_X2_Up
    * \param LX_X2_Dn
    * \param LX_X3_Up
    * \param LX_X3_Dn
    * \param dX1
    * \param dX2
    * \param dX3
    */
    void Reflux_DG( MultiFab&       MF_G,
                    MultiFab&       MF_dU,
                    const Geometry& crse_geom,
                    int             nDOFX,
                    int             nDOFX_X1,
                    int             nDOFX_X2,
                    int             nDOFX_X3,
                    int             nFields,
                    int             iGF_SqrtGm,
                    Array4<int>     NodeNumberTableX_X1,
                    Array4<int>     NodeNumberTableX_X2,
                    Array4<int>     NodeNumberTableX_X3,
                    Array4<Real>    WeightsX_q,
                    Array4<Real>    LX_X1_Up,
                    Array4<Real>    LX_X1_Dn,
                    Array4<Real>    LX_X2_Up,
                    Array4<Real>    LX_X2_Dn,
                    Array4<Real>    LX_X3_Up,
                    Array4<Real>    LX_X3_Dn,
                    Real            dX1,
                    Real            dX2,
                    Real            dX3 );

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
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    void Reflux (MultiFab&       mf,
                 int             dir,
                 Real            scale,
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    /**
     * \brief Overwrite the coarse flux at the coarse/fine interface (and
     * the interface only) with the fine flux stored in the FluxRegister.
     *
     * \param crse_fluxes MultiFab pointers to coarse fluxes.
     * \param scale       scaling factor by which the fine flux is multiplied.
     * \param destcomp    starting component in coarse flux MultiFab
     * \param numcomp     number of components
     * \param crse_geom   coarse Geometry
     */
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

    void Reflux (MultiFab& mf, const MultiFab& volume, Orientation face,
                 Real scale, int scomp, int dcomp, int nc, const Geometry& geom);
    void Reflux_DG( MultiFab&       MF_G,
                    MultiFab&       MF_dU,
                    const Geometry& crse_geom,
                    int             nDOFX,
                    int             nDOFX_X1,
                    int             nDOFX_X2,
                    int             nDOFX_X3,
                    int             nFields,
                    int             iGF_SqrtGm,
                    Array4<int>     NodeNumberTableX_X1,
                    Array4<int>     NodeNumberTableX_X2,
                    Array4<int>     NodeNumberTableX_X3,
                    Array4<Real>    WeightsX_q,
                    Array4<Real>    LX_X1_Up,
                    Array4<Real>    LX_X1_Dn,
                    Array4<Real>    LX_X2_Up,
                    Array4<Real>    LX_X2_Dn,
                    Array4<Real>    LX_X3_Up,
                    Array4<Real>    LX_X3_Dn,
                    Real            dX1,
                    Real            dX2,
                    Real            dX3,
                    Orientation     face );

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
