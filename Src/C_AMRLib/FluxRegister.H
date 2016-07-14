#ifndef _FLUXREGISTER_H_
#define _FLUXREGISTER_H_

#include <BndryRegister.H>

class Geometry;

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
    // The constructor. Note that this takes the fine BoxArray.
    //
    FluxRegister (const BoxArray& fine_boxes, 
                  const IntVect&  ref_ratio,
                  int             fine_lev,
                  int             nvar);
    //
    // The constructor. This version allows setting the DistributionMapping.
    //
    FluxRegister (const BoxArray&            fine_boxes, 
                  const IntVect&             ref_ratio,
                  int                        fine_lev,
                  int                        nvar,
                  const DistributionMapping& dm);
    //
    // The destructor.
    //
    virtual ~FluxRegister ();
    //
    // An enum that says whether to add or copy src data to members.
    //
    enum FrOp {COPY = 0, ADD = 1};
    //
    // Initialize after using default constructor.
    // Note that this takes the fine BoxArray.
    //
    void define (const BoxArray& fine_boxes, 
                 const IntVect&  ref_ratio,
                 int             fine_lev,
                 int             nvar);
    //
    // Initialize after using default constructor.
    // This version allows setting the DistributionMapping.
    //
    void define (const BoxArray&            fine_boxes, 
                 const IntVect&             ref_ratio,
                 int                        fine_lev,
                 int                        nvar,
                 const DistributionMapping& dm);
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
    // Increment flux correction with fine data.
    //
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

    //
    // Write (used for writing to checkpoint)
    //
    void write (const std::string& name, std::ostream& os) const;
    //
    // Read (used for reading from checkpoint)
    //
    void read (const std::string& name, std::istream& is);

    //
    // AddProcsToComp
    //
    virtual void AddProcsToComp(int ioProcNumSCS, int ioProcNumAll,
                                int scsMyId, MPI_Comm scsComm);

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

#endif /*_FLUXREGISTER_H_*/
