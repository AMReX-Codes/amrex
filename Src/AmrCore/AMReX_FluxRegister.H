#ifndef AMREX_FLUXREGISTER_H_
#define AMREX_FLUXREGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_BndryRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_Array.H>

/**
 * \file AMReX_FluxRegister.H
 * \brief Classic flux register used for refluxing between coarse/fine levels.
 */

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
    * \brief Constructor.
    *
    * \param fine_boxes Fine-level BoxArray.
    * \param dm        DistributionMapping describing ownership of fine boxes.
    * \param ref_ratio Refinement ratio between coarse and fine.
    * \param fine_lev  Fine AMR level index.
    * \param nvar      Number of components carried in the registers.
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
    * \brief Initialize after default construction.
    *
    * \param fine_boxes Fine-level BoxArray.
    * \param dm        DistributionMapping describing ownership.
    * \param ref_ratio Refinement ratio between coarse and fine.
    * \param fine_lev  Fine AMR level index.
    * \param nvar      Number of components.
    */
    void define (const BoxArray&            fine_boxes,
                 const DistributionMapping& dm,
                 const IntVect&             ref_ratio,
                 int                        fine_lev,
                 int                        nvar);

    //! Release all internal data structures and reset to the default-constructed state.
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
    * \param mflx    Coarse flux MultiFab.
    * \param area    Area MultiFab on the coarse faces.
    * \param dir     Direction index of the flux.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component stored in the register.
    * \param numcomp Number of components to copy.
    * \param mult    Optional scale factor applied to the flux.
    * \param op      Whether to COPY or ADD the flux into the register.
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
    * \brief Initialize flux correction with coarse data (area already applied).
    *
    * \param mflx    Coarse flux MultiFab with area weighting applied.
    * \param dir     Direction index.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Optional scale factor.
    * \param op      Whether to COPY or ADD the flux into the register.
    */
    void CrseInit (const MultiFab& mflx,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0,
                   FrOp            op = FluxRegister::COPY);

    /**
    * \brief Add coarse fluxes to the register without overwriting existing data.
    *
    * Unlike CrseInit(ADD), this path is used when fine-grid fluxes are unavailable
    * on regions covered by fine levels (e.g., the FLASH code).
    *
    * \param mflx    Coarse flux MultiFab.
    * \param area    Coarse face areas.
    * \param dir     Direction index.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor applied to the flux.
    * \param geom    Coarse Geometry (needed for periodic wrap-around).
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
    * \brief Like CrseAdd above but assumes areas are already multiplied into \p mflx.
    *
    * \param mflx    Coarse flux with area weighting.
    * \param dir     Direction index.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor applied to the flux.
    * \param geom    Coarse Geometry (needed for periodic wrap-around).
    */
    void CrseAdd (const MultiFab& mflx,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult,
                  const Geometry& geom);

    /**
    * \brief Increment flux correction with fine data (areas already applied).
    *
    * \param mflx    Fine flux MultiFab.
    * \param dir     Direction index.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor (typically +1).
    */
    void FineAdd (const MultiFab& mflx,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);

    /**
    * \brief Increment flux correction with fine data, supplying explicit areas.
    *
    * \param mflx    Fine flux MultiFab.
    * \param area    Fine face areas.
    * \param dir     Direction index.
    * \param srccomp Starting component in \p mflx.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor.
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
    * This overload assumes areas have already been multiplied into \p flux.
    *
    * \param flux    Fine flux FAB.
    * \param dir     Direction index.
    * \param boxno   Index of the owning FAB in the register.
    * \param srccomp Starting component in \p flux.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor.
    * \param runon   Execution context (host/device) when copying into the register.
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
    * \brief Increment flux correction with fine flux and area.
    *
    * \param flux    Fine flux FAB.
    * \param area    Area FAB.
    * \param dir     Direction index.
    * \param boxno   Index of the owning FAB in the register.
    * \param srccomp Starting component in \p flux.
    * \param destcomp Starting component in the register.
    * \param numcomp Number of components.
    * \param mult    Scale factor.
    * \param runon   Execution context (host/device) when copying into the register.
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
    * \brief Set flux correction data for a fine box (given by \p boxno) to \p val.
    * This routine used by FLASH does NOT run on GPU for safety.
    *
    * \param dir      Direction index.
    * \param boxno    Index of the owning FAB in the register.
    * \param destcomp Starting component to set.
    * \param numcomp  Number of components to set.
    * \param val      Scalar value to assign.
    * \param runon    Execution context (host/device).
    */
    void FineSetVal (int              dir,
                     int              boxno,
                     int              destcomp,
                     int              numcomp,
                     Real             val,
                     RunOn            runon) noexcept;

    /**
    * \brief Apply flux correction and divide by coarse cell volumes.
    *
    * \param mf        Coarse MultiFab to correct.
    * \param volume    Coarse cell volumes.
    * \param scale     Scaling factor applied to the correction.
    * \param scomp     Starting component within the register.
    * \param dcomp     Starting component within \p mf.
    * \param nc        Number of components to correct.
    * \param crse_geom Coarse Geometry (used for periodic wrap and face orientation).
    */
    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 Real            scale,
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    //! Overload that only corrects faces in direction \p dir using explicit cell volumes.
    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 int             dir,
                 Real            scale,
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    /**
    * \brief Constant-volume version of Reflux() (no explicit volume MultiFab).
    *
    * \param mf        Coarse MultiFab to correct.
    * \param scale     Scaling factor applied to the correction.
    * \param scomp     Starting component within the register.
    * \param dcomp     Starting component within \p mf.
    * \param nc        Number of components to correct.
    * \param crse_geom Coarse Geometry.
    */
    void Reflux (MultiFab&       mf,
                 Real            scale,
                 int             scomp,
                 int             dcomp,
                 int             nc,
                 const Geometry& crse_geom);

    //! Overload that only corrects faces in direction \p dir assuming unit volumes.
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
     * \param srccomp     starting component in fine flux in the FluxRegister
     * \param destcomp    starting component in coarse flux MultiFab
     * \param numcomp     number of components
     * \param crse_geom   coarse Geometry
     */
    void OverwriteFlux (Array<MultiFab*,AMREX_SPACEDIM> const& crse_fluxes,
                        Real scale, int srccomp, int destcomp, int numcomp,
                        const Geometry& crse_geom);


    /**
    * \brief Set internal coarse/fine border fluxes to zero before reuse.
    *
    * \param crse_geom Coarse Geometry describing the domain (needed for iteration).
    */
    void ClearInternalBorders (const Geometry& crse_geom);

    /**
    * \brief Write data to checkpoint files.
    *
    * \param name Base filename prefix.
    * \param os   Stream storing metadata (e.g., BoxArray layout).
    */
    void write (const std::string& name, std::ostream& os) const;

    /**
    * \brief Read data from checkpoint files.
    *
    * \param name Base filename prefix.
    * \param is   Stream providing metadata.
    */
    void read (const std::string& name, std::istream& is);

// public for cuda

    /**
     * \brief Low-level helper used on GPUs: apply correction for a single face orientation.
     *
     * \param mf      Coarse MultiFab to correct.
     * \param volume  Coarse cell volumes.
     * \param face    Orientation being updated.
     * \param scale   Scaling factor applied to the accumulated flux.
     * \param scomp   Starting component within the register.
     * \param dcomp   Starting component within \p mf.
     * \param nc      Number of components to correct.
     * \param geom    Coarse Geometry (used for periodic wrap and face iteration).
     */
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
