
#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H

#include <winstd.H>

#include <stdint.h>

#include <BLassert.H>
#include <FArrayBox.H>
#include <FabArray.H>


//
// Forward declaration.
//
typedef FabArrayId MultiFabId;

class MultiFabCopyDescriptor;

class MultiFab;

namespace BoxLib
{
    void InterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
		       BoxList*                returnUnfilledBoxes,
		       Array<FillBoxId>&       returnedFillBoxIds,
		       const Box&              subbox,
		       MultiFabId              faid1,
		       MultiFabId              faid2,
		       Real                    t1,
		       Real                    t2,
		       Real                    t,
		       int                     src_comp,
		       int                     dest_comp,
		       int                     num_comp,
		       bool                    extrap);

    void InterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
			const Array<FillBoxId>& fillBoxIds,
			MultiFabId              faid1,
			MultiFabId              faid2,
			FArrayBox&              dest,
			Real                    t1,
			Real                    t2,
			Real                    t,
			int                     src_comp,
			int                     dest_comp,
			int                     num_comp,
			bool                    extrap);
}
//
// A Collection of FArrayBoxes
//
// The MultiFab class is publically derived from the
// FabArray<Real,FArrayBox> class.  It is a collection (stored as an array) of
// FArrayBoxes useful for storing floating point data on a domain defined by
// a union of rectangular regions embedded in a uniform index space.  The
// MultiFab class extends the function of the underlying FabArray class just
// as the FArrayBox class extends the funtion of BaseFab<Real>.  Additional
// member functions are defined for I/O and simple arithmetic operations on
// these aggregate objects.
//
// This class does NOT provide a copy constructor or assignment operator.
//
class MultiFab
    :
    public FabArray<FArrayBox>
{
public:
    //
    // Constructs an empty MultiFab.  Data can be defined at a later
    // time using the define member functions inherited
    // from FabArray.
    //
    MultiFab ();
    //
    // Constructs a MultiFab with a valid region defined by bxs and
    // a region of definition defined by the grow factor ngrow.
    // If mem_mode is defined to be Fab_allocate then FArrayBoxes are
    // allocated for each Box in the BoxArray.  The size of the
    // FArrayBox is given by the Box grown by ngrow and the number of
    // components is given by ncomp.  If mem_mode is defined to be
    // Fab_noallocate then no FArrayBoxes are allocated at this time
    // but can be defined later.
    //
    MultiFab (const BoxArray& bs,
              int             ncomp,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate,
	      const IntVect&  nodal = IntVect::TheZeroVector());

    MultiFab (const BoxArray&            bs,
              int                        ncomp,
              int                        ngrow,
              const DistributionMapping& dm,
              FabAlloc                   mem_mode = Fab_allocate,
	      const IntVect&             nodal = IntVect::TheZeroVector());

    MultiFab (const BoxArray& bs,
              int             ncomp,
              int             ngrow,
	      ParallelDescriptor::Color color);

    virtual ~MultiFab() BL_OVERRIDE;


    void operator= (const Real& r);
    //
    // Initialize FArrayBoxes
    //
    void initVal ();
    //
    // Returns the minimum value contained in component comp of the
    // MultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the minimum.  The default is to
    // search only the valid regions of the FArrayBoxes.
    //
    Real min (int comp,
              int nghost = 0,
	      bool local = false) const;
    //
    // Identical to the previous min() function, but confines its
    // search to intersection of Box b and the MultiFab.
    //
    Real min (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    //
    // Returns the maximum value contained in component comp of the
    // MultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the maximum.  The default is to
    // search only the valid regions of the FArrayBoxes.
    //
    Real max (int comp,
              int nghost = 0,
	      bool local = false) const;
    //
    // Identical to the previous max() function, but confines its
    // search to intersection of Box b and the MultiFab.
    //
    Real max (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the MultiFab.
    //
    Real norm0 (int comp = 0, int nghost = 0, bool local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the MultiFab in the intersection of the BoxArray
    // and the valid regions of the MultiFab.
    //
    Real norm0 (int comp, const BoxArray& ba, int nghost = 0, bool local = false) const;
    //
    // Returns the maximum *absolute* values contained in 
    // each component of "comps" of the MultiFab.  No ghost cells are used.
    //
    Array<Real> norm0 (const Array<int>& comps, int nghost = 0, bool local = false) const;
    //
    // Returns the L1 norm of component "comp" over the MultiFab.
    // ngrow ghost cells are used.
    //
    Real norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    //
    // Returns the L1 norm of each component of "comps" over the MultiFab.
    // ngrow ghost cells are used.
    //
    Array<Real> norm1 (const Array<int>& comps, int ngrow = 0, bool local = false) const;
    //
    // Returns the L2 norm of component "comp" over the MultiFab.
    // No ghost cells are used.
    //
    Real norm2 (int comp = 0) const;
    //
    // Returns the L2 norm of each component of "comps" over the MultiFab.
    // No ghost cells are used.
    //
    Array<Real> norm2 (const Array<int>& comps) const;
    //
    // Returns the sum of component "comp" over the MultiFab -- no ghost cells are included.
    //
    Real sum (int comp = 0, bool local = false) const;
    //
    // Adds the scalar value val to the value of each cell in the
    // specified subregion of the MultiFab.  The subregion consists
    // of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each FArrayBox in the subregion that should
    // be modified.  
    //
    void plus (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of plus(), with the
    // restriction that the subregion is further constrained to
    // the intersection with Box region.
    //
    void plus (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the MultiFab.  The value
    // of nghost specifies the number of cells in the boundary
    // region that should be modified.
    //
    void plus (Real val,
               int  nghost);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the MultiFab, that also
    // intersects the Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in
    // the subregion that should be modified.
    //
    void plus (Real       val,
               const Box& region,
               int        nghost);
    //
    // Scales the value of each cell in the specified subregion of the
    // MultiFab by the scalar val (a[i] <- a[i]*val). The subregion
    // consists of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each FArrayBox in the subregion that should
    // be modified.  
    //
    void mult (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of mult(), with the
    // restriction that the subregion is further constrained to the
    // intersection with Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in
    // the subregion that should be modified.
    //
    void mult (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the MultiFab by the scalar val (a[i] <- a[i]*val).
    // The value of nghost specifies the number of cells in the
    // boundary region that should be modified.
    //
    void mult (Real val,
               int  nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the MultiFab by the scalar val (a[i] <- a[i]*val),
    // that also intersects the Box region.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // FArrayBox in the subregion that should be modified.
    //
    void mult (Real       val,
               const Box& region,
               int        nghost = 0);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab with its reciprocal multiplied by the value of
    // numerator. The subregion consists of the num_comp components
    // starting at component comp.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in the
    // subregion that should be modified.
    //
    void invert (Real numerator,
                 int  comp,
                 int  num_comp,
                 int  nghost = 0);
    //
    // Identical to the previous version of invert(), with the
    // restriction that the subregion is further constrained to the
    // intersection with Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in the
    // subregion that should be modified.
    //
    void invert (Real       numerator,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab with its reciprocal multiplied by the value of
    // numerator.  The value of nghost specifies the number of cells
    // in the boundary region that should be modified.
    //
    void invert (Real numerator,
                 int  nghost);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab, that also intersects the Box region, with its
    // reciprocal multiplied by the value of numerator.  The value
    // of nghost specifies the number of cells in the boundary region
    // of each FArrayBox in the subregion that should be modified.
    //
    void invert (Real       numerator,
                 const Box& region,
                 int        nghost);
    //
    // Negates the value of each cell in the specified subregion of
    // the MultiFab.  The subregion consists of the num_comp
    // components starting at component comp.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // FArrayBox in the subregion that should be modified.  
    //
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);
    //
    // Identical to the previous version of negate(), with the
    // restriction that the subregion is further constrained to
    // the intersection with Box region.
    //
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Negates the value of each cell in the valid region of
    // the MultiFab.  The value of nghost specifies the number of
    // cells in the boundary region that should be modified.  
    //
    void negate (int nghost = 0);
    //
    // Negates the value of each cell in the valid region of
    // the MultiFab that also intersects the Box region.  The value
    // of nghost specifies the number of cells in the boundary region
    // that should be modified.  
    //
    void negate (const Box& region,
                 int        nghost = 0);

    IntVect minIndex (int comp,
                      int nghost = 0) const;

    IntVect maxIndex (int comp,
                      int nghost = 0) const;
    //
    // This function adds the values of the cells in mf to the corresponding
    // cells of this MultiFab.  mf is required to have the same BoxArray or
    // "valid region" as this MultiFab.  The addition is done only to num_comp
    // components, starting with component number strt_comp.  The parameter
    // nghost specifies the number of boundary cells that will be modified.
    // If nghost == 0, only the valid region of each FArrayBox will be
    // modified.
    //
    void plus (const MultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);
    //
    //
    // This function subtracts the values of the cells in mf from the
    // corresponding cells of this MultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this MultiFab.  The subtraction is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each FArrayBox will be modified.
    //
    void minus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);
    //
    //
    // This function divides the values of the cells in mf from the
    // corresponding cells of this MultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this MultiFab.  The division is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each FArrayBox will be modified.  Note, nothing is done to protect 
    // against divide by zero.
    //
    void divide (const MultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);
    //
    // Add src to dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Add (MultiFab&       dst,
		     const MultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);
    //
    // Copy from src to dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    // The copy is local.  The parallel copy function is in FabArray.H
    //
    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);
    //
    // Subtract src from dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Subtract (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Multiply dst by src including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Multiply (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Divide dst by src including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Divide (MultiFab&       dst,
			const MultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);
    //
    // dst += a*src
    //
    static void Saxpy (MultiFab&       dst,
		       Real            a, 
		       const MultiFab& src,
		       int             srccomp,
		       int             dstcomp,
		       int             numcomp,
		       int             nghost);

    void define (const BoxArray& bxs,
		 int             nvar,
		 int             ngrow,
		 FabAlloc        mem_mode,
		 const IntVect&  nodal = IntVect::TheZeroVector(),
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());

    void define (const BoxArray&            bxs,
		 int                        nvar,
		 int                        ngrow,
		 const DistributionMapping& dm,
		 FabAlloc                   mem_mode,
		 const IntVect&             nodal = IntVect::TheZeroVector());

    //
    // Are there any NaNs in the MF?
    // This may return false, even if the MF contains NaNs, if the machine
    // doesn't support the appropriate NaN testing functions.
    //
    // This version checks all components and grow cells.
    //
    bool contains_nan (bool local=false) const;

    bool contains_nan (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    //
    // Are there any Infs in the MF?
    // This may return false, even if the MF contains Infs, if the machine
    // doesn't support the appropriate Inf testing functions.
    //
    // This version checks all components and grow cells.
    //
    bool contains_inf (bool local=false) const;

    //
    // This tests on whether the MultiFab is fully nodal.
    //
    bool is_nodal () const;

    bool contains_inf (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    //
    // Sum ghost cells that are covered by valid cells into the valid cells.
    //
    void SumBoundary (const Periodicity& period = Periodicity::NonPeriodic());

    void SumBoundary (int scomp, int ncomp, const Periodicity& period = Periodicity::NonPeriodic());

    static void Initialize ();
    static void Finalize ();

    virtual void AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                                 int scsMyId, MPI_Comm scsComm);

private:
    //
    // These are disabled.
    //
    MultiFab (const MultiFab& rhs);
    MultiFab& operator= (const MultiFab& rhs);

    //
    // Some useful typedefs.
    //
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef FabArrayBase::MapOfCopyComTagContainers MapOfCopyComTagContainers;
};

class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<FArrayBox>
{
  public:

    MultiFabCopyDescriptor ();

    ~MultiFabCopyDescriptor ();

    MultiFabId RegisterMultiFab (MultiFab* mf) { return RegisterFabArray(mf); }

  private:
    //
    // These are disallowed.
    //
    MultiFabCopyDescriptor (const MultiFabCopyDescriptor&);
    MultiFabCopyDescriptor& operator= (const MultiFabCopyDescriptor&);
};

#endif /*BL_MULTIFAB_H*/
