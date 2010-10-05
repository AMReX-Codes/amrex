#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H
//
// $Id: MultiFab.H,v 1.55 2010-10-05 23:00:00 lijewski Exp $
//
#include <winstd.H>

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
    void linInterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
                          BoxList*                returnUnfilledBoxes,
                          Array<FillBoxId>&       returnedFillBoxIds,
                          const Box&              subbox,
                          const MultiFabId&       faid1,
                          const MultiFabId&       faid2,
                          Real                    t1,
                          Real                    t2,
                          Real                    t,
                          int                     src_comp,
                          int                     dest_comp,
                          int                     num_comp,
                          bool                    extrap);

    void linInterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
                           const Array<FillBoxId>& fillBoxIds,
                           const MultiFabId&       faid1,
                           const MultiFabId&       faid2,
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
//@Man:
//@Memo: A Collection of FArrayBoxes
/*@Doc:

  The MultiFab class is publically derived from the
  FabArray<Real,FArrayBox> class.  It is a collection (stored as an array) of
  FArrayBoxes useful for storing floating point data on a domain defined by
  a union of rectangular regions embedded in a uniform index space.  The
  MultiFab class extends the function of the underlying FabArray class just
  as the FArrayBox class extends the funtion of BaseFab<Real>.  Additional
  member functions are defined for I/O and simple arithmetic operations on
  these aggregate objects.

  This class does NOT provide a copy constructor or assignment operator.
*/
class MultiFab
    :
    public FabArray<FArrayBox>
{
public:

    /*@ManDoc: Constructs an empty MultiFab.  Data can be defined at a later
               time using the `define' member functions inherited
               from FabArray.
    */
    MultiFab ();

    /*@ManDoc: Constructs a MultiFab with a valid region defined by bxs and
               a region of definition defined by the grow factor ngrow.
               If mem\_mode is defined to be Fab\_allocate then FArrayBoxes are
               allocated for each Box in the BoxArray.  The size of the
               FArrayBox is given by the Box grown by ngrow and the number of
               components is given by ncomp.  If mem\_mode is defined to be
               Fab\_noallocate then no FArrayBoxes are allocated at this time
               but can be defined later.
    */
    MultiFab (const BoxArray& bs,
              int             ncomp,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);

    void operator= (const Real& r);

    /*@ManDoc: Returns the minimum value contained in component comp of the
               MultiFab.  The parameter nghost determines the number of
               boundary cells to search for the minimum.  The default is to
               search only the valid regions of the FArrayBoxes.
    */
    Real min (int comp,
              int nghost = 0) const;

    /*@ManDoc: Identical to the previous min() function, but confines its
               search to intersection of Box b and the MultiFab.
    */
    Real min (const Box& b,
              int        comp,
              int        nghost = 0) const;

    /*@ManDoc: Returns the maximum value contained in component comp of the
               MultiFab.  The parameter nghost determines the number of
               boundary cells to search for the maximum.  The default is to
               search only the valid regions of the FArrayBoxes.
    */
    Real max (int comp,
              int nghost = 0) const;

    /*@ManDoc: Returns the maximum *absolute* value contained in 
               component comp of the MultiFab.  No ghost cells are used.
    */
    Real norm0 (int comp = 0) const;

    /*@ManDoc: Returns the L2 norm of component "comp" over the MultiFab.
               No ghost cells are used.
    */
    Real norm2 (int comp = 0) const;

    /*@ManDoc: Identical to the previous max() function, but confines its
               search to intersection of Box b and the MultiFab.
    */
    Real max (const Box& b,
              int        comp,
              int        nghost = 0) const;

    /*@ManDoc: Adds the scalar value val to the value of each cell in the
               specified subregion of the MultiFab.  The subregion consists
               of the num\_comp components starting at component comp.
               The value of nghost specifies the number of cells in the
               boundary region of each FArrayBox in the subregion that should
               be modified.  
    */
    void plus (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);

    /*@ManDoc: Identical to the previous version of plus(), with the
               restriction that the subregion is further constrained to
               the intersection with Box region.
    */
    void plus (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);

    /*@ManDoc: Adds the scalar value val to the value of each cell in the
               valid region of each component of the MultiFab.  The value
               of nghost specifies the number of cells in the boundary
               region that should be modified.
    */
    void plus (Real val,
               int  nghost);

    /*@ManDoc: Adds the scalar value val to the value of each cell in the
               valid region of each component of the MultiFab, that also
               intersects the Box region.  The value of nghost specifies the
               number of cells in the boundary region of each FArrayBox in
               the subregion that should be modified.

    */
    void plus (Real       val,
               const Box& region,
               int        nghost);

    /*@ManDoc: Scales the value of each cell in the specified subregion of the
               MultiFab by the scalar val (a[i] <- a[i]*val). The subregion
               consists of the num\_comp components starting at component comp.
               The value of nghost specifies the number of cells in the
               boundary region of each FArrayBox in the subregion that should
               be modified.  
    */
    void mult (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);

    /*@ManDoc: Identical to the previous version of mult(), with the
               restriction that the subregion is further constrained to the
               intersection with Box region.  The value of nghost specifies the
               number of cells in the boundary region of each FArrayBox in
               the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);

    /*@ManDoc: Scales the value of each cell in the valid region of each
               component of the MultiFab by the scalar val (a[i] <- a[i]*val).
               The value of nghost specifies the number of cells in the
               boundary region that should be modified.
    */
    void mult (Real val,
               int  nghost = 0);

    /*@ManDoc: Scales the value of each cell in the valid region of each
               component of the MultiFab by the scalar val (a[i] <- a[i]*val),
               that also intersects the Box region.  The value of nghost
               specifies the number of cells in the boundary region of each
               FArrayBox in the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        nghost = 0);

    /*@ManDoc: Replaces the value of each cell in the specified subregion of
               the MultiFab with its reciprocal multiplied by the value of
               numerator. The subregion consists of the num\_comp components
               starting at component comp.  The value of nghost specifies the
               number of cells in the boundary region of each FArrayBox in the
               subregion that should be modified.
    */
    void invert (Real numerator,
                 int  comp,
                 int  num_comp,
                 int  nghost = 0);

    /*@ManDoc: Identical to the previous version of invert(), with the
               restriction that the subregion is further constrained to the
               intersection with Box region.  The value of nghost specifies the
               number of cells in the boundary region of each FArrayBox in the
               subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Replaces the value of each cell in the specified subregion of
               the MultiFab with its reciprocal multiplied by the value of
               numerator.  The value of nghost specifies the number of cells
               in the boundary region that should be modified.
    */
    void invert (Real numerator,
                 int  nghost);

    /*@ManDoc: Replaces the value of each cell in the specified subregion of
               the MultiFab, that also intersects the Box region, with its
               reciprocal multiplied by the value of numerator.  The value
               of nghost specifies the number of cells in the boundary region
               of each FArrayBox in the subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        nghost);

    /*@ManDoc: Negates the value of each cell in the specified subregion of
               the MultiFab.  The subregion consists of the num\_comp
               components starting at component comp.  The value of nghost
               specifies the number of cells in the boundary region of each
               FArrayBox in the subregion that should be modified.  
    */
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);

    /*@ManDoc: Identical to the previous version of negate(), with the
               restriction that the subregion is further constrained to
               the intersection with Box region.
    */
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Negates the value of each cell in the valid region of
               the MultiFab.  The value of nghost specifies the number of
               cells in the boundary region that should be modified.  
    */
    void negate (int nghost = 0);

    /*@ManDoc: Negates the value of each cell in the valid region of
               the MultiFab that also intersects the Box region.  The value
               of nghost specifies the number of cells in the boundary region
               that should be modified.  
    */
    void negate (const Box& region,
                 int        nghost = 0);

    /*@ManDoc:
      This function adds the values of the cells in mf to the corresponding
      cells of this MultiFab.  mf is required to have the same BoxArray or
      "valid region" as this MultiFab.  The addition is done only to num\_comp
      components, starting with component number strt\_comp.  The parameter
      nghost specifies the number of boundary cells that will be modified.
      If nghost == 0, only the valid region of each FArrayBox will be
      modified.
    */
    void plus (const MultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);

    /*@ManDoc:
      This function subtracts the values of the cells in mf from the
      corresponding cells of this MultiFab.  mf is required to have the
      same BoxArray or "valid region" as this MultiFab.  The subtraction is
      done only to num\_comp components, starting with component number
      strt\_comp.  The parameter nghost specifies the number of boundary
      cells that will be modified.  If nghost == 0, only the valid region of
      each FArrayBox will be modified.
    */
    void minus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);

    /*@ManDoc:
      This function divides the values of the cells in mf from the
      corresponding cells of this MultiFab.  mf is required to have the
      same BoxArray or "valid region" as this MultiFab.  The division is
      done only to num\_comp components, starting with component number
      strt\_comp.  The parameter nghost specifies the number of boundary
      cells that will be modified.  If nghost == 0, only the valid region of
      each FArrayBox will be modified.  Note, nothing is done to protect 
      against divide by zero.
    */
    void divide (const MultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);

    /*@ManDoc: Add 'src' to `dst' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Add (MultiFab&       dst,
		     const MultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);

    /*@ManDoc: Copy from `src' to `dst' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    /*@ManDoc: Subtract 'src' from `dst' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Subtract (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    /*@ManDoc: Multiply 'dst' by `src' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Multiply (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    /*@ManDoc: Divide 'dst' by `src' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Divide (MultiFab&       dst,
			const MultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);

    void define (const BoxArray& bxs,
                 int             nvar,
                 int             ngrow,
                 FabAlloc        mem_mode);

    void define (const BoxArray&            bxs,
                 int                        nvar,
                 int                        ngrow,
		 const DistributionMapping& dm,
                 FabAlloc                   mem_mode);
    //
    // Are there any NaNs in the MF?
    // This may return false, even if the MF contains NaNs, if the machine
    // doesn't support the appropriate NaN testing functions.
    //
    bool contains_nan () const;
    //
    // Are there any Infs in the MF?
    // This may return false, even if the MF contains Infs, if the machine
    // doesn't support the appropriate Inf testing functions.
    //
    bool contains_inf () const;

    const FArrayBox& operator[] (int K) const;

    FArrayBox& operator[] (int K);

    const FArrayBox& operator[] (const MFIter& mfi) const { return operator[](mfi.index()); }

    const FArrayBox& get (const MFIter& mfi) const { return operator[](mfi.index()); }

    FArrayBox& operator[] (const MFIter& mfi) { return operator[](mfi.index()); }

    FArrayBox& get (const MFIter& mfi) { return operator[](mfi.index()); }
    //
    // Set ghost cells that are covered by valid cells to the value of valid cells.
    //
    void FillBoundary (bool local = false);

    void FillBoundary (int scomp, int ncomp, bool local = false);
    //
    // Sum ghost cells that are covered by valid cells into the valid cells.
    //
    void SumBoundary ();

    void SumBoundary (int scomp, int ncomp);

    static void Initialize ();

protected:

    static bool check_for_nan;
    static bool check_for_inf;

private:
    //
    // These are disabled.
    //
    MultiFab (const MultiFab& rhs);
    MultiFab& operator= (const MultiFab& rhs);
};

class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<FArrayBox>
{
  public:

    MultiFabCopyDescriptor ()
        :
        FabArrayCopyDescriptor<FArrayBox>() {}

    MultiFabId RegisterMultiFab (MultiFab* mf)
    {
        return RegisterFabArray(mf);
    }

  private:
    //
    // These are disallowed.
    //
    MultiFabCopyDescriptor (const MultiFabCopyDescriptor&);
    MultiFabCopyDescriptor& operator= (const MultiFabCopyDescriptor&);
};

#endif /*BL_MULTIFAB_H*/
