
#ifndef _FLUXREGISTER_H_
#define _FLUXREGISTER_H_

//
// $Id: FluxRegister.H,v 1.11 2000-10-02 20:48:42 lijewski Exp $
//

#include <BndryRegister.H>

class Geometry;

//
//@Man:
//@Memo: Flux Register.
/*@Doc:

  Stores and manipulates fluxes at coarse-fine interfaces.
*/

class FluxRegister
    :
    public BndryRegister
{
public:
    //
    //@ManDoc: The default constructor.
    //
    FluxRegister();

    /*@ManDoc: The constructor.
               Note that this takes the fine BoxArray.
    */
    FluxRegister (const BoxArray& fine_boxes, 
                  const IntVect&  ref_ratio,
                  int             fine_lev,
                  int             nvar);
    //
    //@ManDoc: The destructor.
    //
    virtual ~FluxRegister ();

    /*@ManDoc: Initialize after using default constructor.
               Note that this takes the fine BoxArray.
    */
    virtual void define (const BoxArray& fine_boxes, 
                         const IntVect&  ref_ratio,
                         int             fine_lev,
                         int             nvar);
    //
    //@ManDoc: Returns the refinement ratio.
    //
    const IntVect& refRatio () const { return ratio; }
    //
    //@ManDoc: Returns the level number of the fine level.
    //
    int fineLevel () const { return fine_level; }
    //
    //@ManDoc: Returns the level number of the coarse level (fineLevel()-1).
    //
    int crseLevel () const { return fine_level-1; }
    //
    //@ManDoc: The number of components.
    //
    int nComp () const { return ncomp; }
    //
    //@ManDoc: The coarsened boxes.
    //
    const BoxArray& coarsenedBoxes () const { return grids; }
    //
    //@ManDoc: Returns the sum of the registers.
    //
    Real SumReg (int comp) const;
    //
    //@ManDoc: Copy from flux register to FAB `flx'.
    //
    void copyTo (FArrayBox& flx,
                 int        dir,
                 int        src_comp,
                 int        dest_comp,
                 int        num_comp);
    //
    //@ManDoc: flx = flx + mult*register.
    //
    void addTo (FArrayBox& flx,
                int        dir,
                Real       mult,
                int        src_comp,
                int        dest_comp,
                int        num_comp);
    //
    //@ManDoc: mf = fm + mult*register.
    //
    void addTo (MultiFab& mf,
                int       dir,
                Real      mult,
                int       src_comp,
                int       dest_comp,
                int       num_comp);
    //
    //@ManDoc: Initialize flux correction with coarse data.
    //
    void CrseInit (const MultiFab& mflx,
                   const MultiFab& area,
                   int             dir,
                   int             srccomp,
                   int             destcomp,
                   int             numcomp,
                   Real            mult = -1.0);
    //
    //@ManDoc: Initialize flux correction with coarse data. 
    //
    void CrseInit (const FArrayBox& flux,
                   const Box&       subbox,
                   int              dir,
                   int              srccomp,
                   int              destcomp,
                   int              numcomp,
                   Real             mult = -1.0);
    //
    // This isn't documented  --  it's a parallel feature :-(
    //
    // It MUST be called after all CrseInit()s outside of
    // any implicit or explicit parallel loop.
    //
    void CrseInitFinish ();
    //
    //@ManDoc: Increment flux correction with fine data.
    //
    void FineAdd (const MultiFab& mflx,
                  int             dir, 
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);
    //
    //@ManDoc: Increment flux correction with fine data.
    //
    void FineAdd (const MultiFab& mflx,
                  const MultiFab& area,
                  int             dir, 
                  int             srccomp,
                  int             destcomp,
                  int             numcomp,
                  Real            mult);
    //
    //@ManDoc: Increment flux correction with fine data.
    //
    void FineAdd (const FArrayBox& flux,
                  int              dir,
                  int              boxno,
                  int              srccomp,
                  int              destcomp,
                  int              numcomp,
                  Real             mult);
    //
    //@ManDoc: Increment flux correction with fine data.
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
    //@ManDoc: Apply flux correction.  Note that this takes the coarse Geometry.
    //
    void Reflux (MultiFab&       mf,
                 const MultiFab& volume,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);
    //
    //@ManDoc: Constant volume version of Reflux().  Note that this takes the coarse Geometry.
    //
    void Reflux (MultiFab&       mf,
                 Real            scale,
                 int             srccomp,
                 int             destcomp,
                 int             numcomp,
                 const Geometry& crse_geom);

protected:
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
    //
    // Helper member function.
    //
    void increment (const FArrayBox& fab, int dir);
};

#endif /*_FLUXREGISTER_H_*/
