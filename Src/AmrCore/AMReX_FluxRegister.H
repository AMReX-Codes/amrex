#ifndef AMREX_FLUXREGISTER_H_
#define AMREX_FLUXREGISTER_H_

#include <AMReX_BndryRegister.H>
#include <AMReX_Geometry.H>

namespace amrex {

//
// Flux Register.
//
// Stores and manipulates fluxes at coarse-fine interfaces.
//

class FluxRegister
    :
    public BndryRegister
{
public:
    //
    // The default constructor.
    //
    FluxRegister();
    //
    // The constructor. This version allows setting the DistributionMapping.
    //
    FluxRegister (const BoxArray&            fine_boxes, 
                  const DistributionMapping& dm,
                  const IntVect&             ref_ratio,
                  int                        fine_lev,
                  int                        nvar);
    //
    // The destructor.
    //
    virtual ~FluxRegister ();

    FluxRegister (FluxRegister&& rhs) noexcept = default;

    FluxRegister (const FluxRegister& rhs) = delete;
    FluxRegister& operator= (const FluxRegister& rhs) = delete;
    FluxRegister& operator= (FluxRegister&& rhs) = delete;

    //
    // An enum that says whether to add or copy src data to members.
    //
    enum FrOp {COPY = 0, ADD = 1};
    //
    // Initialize after using default constructor.
    // This version allows setting the DistributionMapping.
    //
    void define (const BoxArray&            fine_boxes, 
                 const DistributionMapping& dm,
                 const IntVect&             ref_ratio,
                 int                        fine_lev,
                 int                        nvar);

    void clear ();

    //
    // Returns the refinement ratio.
    //
    const IntVect& refRatio () const;
    //
    // Returns the level number of the fine level.
    //
    int fineLevel () const;
    //
    // Returns the level number of the coarse level (fineLevel()-1).
    //
    int crseLevel () const;
    //
    // The number of components.
    //
    int nComp () const;
    //
    // The coarsened boxes.
    //
    const BoxArray& coarsenedBoxes () const;
    //
    // Returns the sum of the registers.
    //
    Real SumReg (int comp) const;
    //
    // Initialize flux correction with coarse data.
    //
    void CrseInit (const MultiFab& mflx,
                   const MultiFab& area,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0,
                   FrOp            op = FluxRegister::COPY);
    //
    // Initialize flux correction with coarse data.
    //
    void CrseInit (const MultiFab& mflx,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0,
                   FrOp            op = FluxRegister::COPY);
    //
    //  Add coarse fluxes to the flux register.
    //  This is different from CrseInit with FluxRegister::ADD.
    //  This is used for cases in which the grids covered by fine do not have fluxes computed.
    //
    void CrseAdd (const MultiFab& mflx,
                  const MultiFab& area,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult,
                  const Geometry& geom);
  ///in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
    void CrseAdd (const MultiFab& mflx,
                  int             dir,
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult,
                  const Geometry& geom);
    //
    // Increment flux correction with fine data.
    //
    ///in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
  void FineAdd (const MultiFab& mflx,
                  int             dir, 
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);
    //
    // Increment flux correction with fine data.
    //
    void FineAdd (const MultiFab& mflx,
                  const MultiFab& area,
                  int             dir, 
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);
    //
    // Increment flux correction with fine data.
    //
  ///in this version the area is assumed to muliplied into the flux (if not, use scale to fix)
    void FineAdd (const FArrayBox& flux,
                  int              dir,
                  int              boxno,
                  int              srccomp,
                  int              destcomp,
                  int              numcomp,
                  Real             mult);
    //
    // Increment flux correction with fine data.
    //
    void FineAdd (const FArrayBox& flux,
                  const FArrayBox& area,
                  int              dir,
                  int              boxno,
                  int              srccomp,
                  int              destcomp,
                  int              numcomp,
                  Real             mult);
    //
    // Apply flux correction.  Note that this takes the coarse Geometry.
    //
    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);
    //
    // Constant volume version of Reflux().  Note that this takes the coarse Geometry.
    //
    void Reflux (MultiFab&       mf,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

    // Set internal borders to zero
    void ClearInternalBorders (const Geometry& crse_geom);
    //
    // Write (used for writing to checkpoint)
    //
    void write (const std::string& name, std::ostream& os) const;
    //
    // Read (used for reading from checkpoint)
    //
    void read (const std::string& name, std::istream& is);

private:
    //
    // Helper member function.
    //
    void increment (const FArrayBox& fab, int dir);
    //
    // Refinement ratio
    //
    IntVect ratio;
    //
    // Current level + 1.
    //
    int fine_level;
    //
    // Number of state components.
    //
    int ncomp;
};

}

#endif /*_FLUXREGISTER_H_*/
