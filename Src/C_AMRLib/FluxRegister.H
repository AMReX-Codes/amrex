//BL_COPYRIGHT_NOTICE

#ifndef _FLUXREGISTER_H_
#define _FLUXREGISTER_H_

//
// $Id: FluxRegister.H,v 1.3 1997-11-26 19:18:49 lijewski Exp $
//

#include <BndryRegister.H>

class Geometry;

/// Stores and manipulates fluxes at coarse-fine interfaces
class FluxRegister  : public BndryRegister
{
 /** This class is a collection of storage, access, and utility
     functions to deal with fluxes at the interface between
     coarse and fine grids
 */
public:
    //@ManDoc: constructor
    FluxRegister();
    //@ManDoc: constructor
    FluxRegister(const BoxArray& fine_boxes, 
		 const IntVect & ref_ratio, int fine_lev, int nvar);

//@ManDoc: destructor
    virtual ~FluxRegister();

///
//@ManDoc: data member initialization
    virtual void define(const BoxArray& fine_boxes, 
		        const IntVect & ref_ratio, int fine_lev, int nvar);

      /// access functions
///
    const IntVect & refRatio() const {
	return ratio;
    }
///
    int fineLevel() const {
	return fine_level;
    }
///
    int crseLevel() const {
	return fine_level-1;
    }
///
    int nComp() const {
	return ncomp;
    }
///
    const BoxArray& coarsenedBoxes() const {
	return grids;
    }
///

      /// utility functions

    Real SumReg(int comp) const;

      //@ManDoc: copy from flux register to fab
    void copyTo(FArrayBox& flx, int dir,
	      int src_comp, int dest_comp, int num_comp);

      //@ManDoc: copy from flux register to multifab
    void copyTo(MultiFab& mf, int dir,
	      int src_comp, int dest_comp, int num_comp);

      //@ManDoc: flx = flx + mult*register
    void addTo(FArrayBox& flx, int dir, Real mult,
	       int src_comp, int dest_comp, int num_comp);

      //@ManDoc: mf = fm + mult*register
    void addTo(MultiFab& mf, int dir, Real mult,
	       int src_comp, int dest_comp, int num_comp);

      //@ManDoc: flx = flx + mult*register/area
    void scaleAddTo(FArrayBox& flx, const FArrayBox& area,
                    int dir, Real mult,
                    int src_comp, int dest_comp, int num_comp);

      //@ManDoc: mf = fm + mult*register/area
    void scaleAddTo(MultiFab& mf, const MultiFab& area,
                    int dir, Real mult,
                    int src_comp, int dest_comp, int num_comp);

    ///coarse-fine fixup

///
    void CrseInit(const MultiFab& mflx, int dir,
                  int srccomp, int destcomp, int numcomp,
                  Real mult = -1.0);

///
    void CrseInit(const MultiFab& mflx, const MultiFab& area, int dir,
                  int srccomp, int destcomp, int numcomp,
                  Real mult = -1.0);

///
    void CrseInit(const FArrayBox& flux, const Box& subbox, int dir,
                  int srccomp, int destcomp, int numcomp,
                  Real mult = -1.0);

    void CrseInitFinish();
///
    void CrseInit(const FArrayBox& flux, const FArrayBox& area,
                  const Box& subbox, int dir,
                  int srccomp, int destcomp, int numcomp,
                  Real mult = -1.0);

///
    void FineAdd(const MultiFab& mflx, int dir, 
		 int srccomp, int destcomp, int numcomp, Real mult);

///
    void FineAdd(const MultiFab& mflx, const MultiFab& area, int dir, 
		 int srccomp, int destcomp, int numcomp, Real mult);

///
    void FineAdd(const FArrayBox& flux, int dir, int boxno,
		 int srccomp, int destcomp, int numcomp, Real mult);

///
    void FineAdd(const FArrayBox& flux, const FArrayBox& area,
                 int dir, int boxno,
		 int srccomp, int destcomp, int numcomp, Real mult);

///
    //@ManDoc: reflux
    void Reflux(MultiFab& mf, const MultiFab& volume, Real scale,
                int srccomp, int destcomp, int numcomp, const Geometry&geom);

      //@ManDoc: constant volume version of reflux
    void Reflux(MultiFab& mf, Real scale,
                int srccomp, int destcomp, int numcomp, const Geometry&geom);

//@ManDoc: i-o function
    void print(ostream& os);

protected:
      
//@ManDoc:refinement ratio
    IntVect ratio;
//@ManDoc: current level + 1
    int fine_level;
//@ManDoc: number of state components
    int ncomp;
///

      /// helper member functions
///
    void increment(const FArrayBox& fab, int dir);
};

#endif
