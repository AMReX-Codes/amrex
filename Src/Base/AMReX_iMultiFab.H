
#ifndef BL_IMULTIFAB_H
#define BL_IMULTIFAB_H

#include <AMReX_BLassert.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>

namespace amrex {

//
// A Collection of IArrayBoxes
//
// The iMultiFab class is publically derived from the
// FabArray<int,IArrayBox> class.  It is a collection (stored as an array) of
// IArrayBoxes useful for storing integer data on a domain defined by
// a union of rectangular regions embedded in a uniform index space.  The
// iMultiFab class extends the function of the underlying FabArray class just
// as the IArrayBox class extends the funtion of BaseFab<int>.  Additional
// member functions are defined for I/O and simple arithmetic operations on
// these aggregate objects.
//
// This class does NOT provide a copy constructor or assignment operator.
//
class iMultiFab
    :
    public FabArray<IArrayBox>
{
public:
    //
    // Constructs an empty iMultiFab.  Data can be defined at a later
    // time using the define member functions inherited
    // from FabArray.
    //
    iMultiFab ();
    //
    // Constructs a iMultiFab with a valid region defined by bxs and
    // a region of definition defined by the grow factor ngrow.
    //
    iMultiFab (const BoxArray&            bs,
	       const DistributionMapping& dm,
	       int                        ncomp,
	       int                        ngrow,
#ifdef AMREX_STRICT_MODE
               const MFInfo&              info,
               const FabFactory<IArrayBox>& factory);
#else
               const MFInfo&              info = MFInfo(),
               const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>());
#endif

    /**
     * \brief Make an alias iMultiFab. maketype must be
     * amrex::make_alias.  scomp is the starting component of the
     * alias and ncomp is the number of components in the new aliasing
     * iMultiFab.
     */
    iMultiFab (const iMultiFab& rhs, MakeType maketype, int scomp, int ncomp);

    virtual ~iMultiFab () override = default;
    
    iMultiFab (iMultiFab&& rhs) noexcept = default;
    iMultiFab& operator= (iMultiFab&& rhs) noexcept = delete;

    iMultiFab (const iMultiFab& rhs) = delete;
    iMultiFab& operator= (const iMultiFab& rhs) = delete;

    void operator= (const int& r);
    //
    // Returns the minimum value contained in component comp of the
    // iMultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the minimum.  The default is to
    // search only the valid regions of the IArrayBoxes.
    //
    int min (int comp,
             int nghost = 0,
	     bool local = false) const;
    //
    // Identical to the previous min() function, but confines its
    // search to intersection of Box b and the iMultiFab.
    //
    int min (const Box& b,
             int        comp,
             int        nghost = 0,
	     bool       local = false) const;
    //
    // Returns the maximum value contained in component comp of the
    // iMultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the maximum.  The default is to
    // search only the valid regions of the IArrayBoxes.
    //
    int max (int comp,
	     int nghost = 0,
	     bool local = false) const;
    //
    // Identical to the previous max() function, but confines its
    // search to intersection of Box b and the iMultiFab.
    //
    int max (const Box& b,
             int        comp,
             int        nghost = 0,
	     bool       local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the iMultiFab.
    //
    int norm0 (int comp = 0, int nghost = 0, bool local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the iMultiFab in the intersection of the BoxArray
    // and the valid regions of the iMultiFab.
    //
    int norm0 (int comp, const BoxArray& ba, int nghost = 0, bool local = false) const;

    //
    // Returns the L1 norm of component "comp" over the iMultiFab.
    // ngrow ghost cells are used.
    //
    int norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    
    //
    // Returns the L2 norm of component "comp" over the iMultiFab.
    // No ghost cells are used.
    //
    int norm2 (int comp = 0) const;
    //
    // Adds the scalar value val to the value of each cell in the
    // specified subregion of the iMultiFab.  The subregion consists
    // of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each IArrayBox in the subregion that should
    // be modified.  
    //
    void plus (int val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of plus(), with the
    // restriction that the subregion is further constrained to
    // the intersection with Box region.
    //
    void plus (int       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the iMultiFab.  The value
    // of nghost specifies the number of cells in the boundary
    // region that should be modified.
    //
    void plus (int val,
               int  nghost);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the iMultiFab, that also
    // intersects the Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each IArrayBox in
    // the subregion that should be modified.
    //
    void plus (int       val,
               const Box& region,
               int        nghost);
    //
    // Scales the value of each cell in the specified subregion of the
    // iMultiFab by the scalar val (a[i] <- a[i]*val). The subregion
    // consists of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each IArrayBox in the subregion that should
    // be modified.  
    //
    void mult (int val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of mult(), with the
    // restriction that the subregion is further constrained to the
    // intersection with Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each IArrayBox in
    // the subregion that should be modified.
    //
    void mult (int       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the iMultiFab by the scalar val (a[i] <- a[i]*val).
    // The value of nghost specifies the number of cells in the
    // boundary region that should be modified.
    //
    void mult (int val,
               int  nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the iMultiFab by the scalar val (a[i] <- a[i]*val),
    // that also intersects the Box region.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // IArrayBox in the subregion that should be modified.
    //
    void mult (int       val,
               const Box& region,
               int        nghost = 0);
    //
    // Negates the value of each cell in the specified subregion of
    // the iMultiFab.  The subregion consists of the num_comp
    // components starting at component comp.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // IArrayBox in the subregion that should be modified.  
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
    // the iMultiFab.  The value of nghost specifies the number of
    // cells in the boundary region that should be modified.  
    //
    void negate (int nghost = 0);
    //
    // Negates the value of each cell in the valid region of
    // the iMultiFab that also intersects the Box region.  The value
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
    // cells of this iMultiFab.  mf is required to have the same BoxArray or
    // "valid region" as this iMultiFab.  The addition is done only to num_comp
    // components, starting with component number strt_comp.  The parameter
    // nghost specifies the number of boundary cells that will be modified.
    // If nghost == 0, only the valid region of each IArrayBox will be
    // modified.
    //
    void plus (const iMultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);
    //
    //
    // This function subtracts the values of the cells in mf from the
    // corresponding cells of this iMultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this iMultiFab.  The subtraction is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each IArrayBox will be modified.
    //
    void minus (const iMultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);
    //
    //
    // This function divides the values of the cells in mf from the
    // corresponding cells of this iMultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this iMultiFab.  The division is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each IArrayBox will be modified.  Note, nothing is done to protect 
    // against divide by zero.
    //
    void divide (const iMultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);
    //
    // Add src to dst including nghost ghost cells.
    // The two iMultiFabs MUST have the same underlying BoxArray.
    //
    static void Add (iMultiFab&       dst,
		     const iMultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);
    //
    // Copy from src to dst including nghost ghost cells.
    // The two iMultiFabs MUST have the same underlying BoxArray.
    //
    static void Copy (iMultiFab&       dst,
                      const iMultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);
    //
    // Subtract src from dst including nghost ghost cells.
    // The two iMultiFabs MUST have the same underlying BoxArray.
    //
    static void Subtract (iMultiFab&       dst,
			  const iMultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Multiply dst by src including nghost ghost cells.
    // The two iMultiFabs MUST have the same underlying BoxArray.
    //
    static void Multiply (iMultiFab&       dst,
			  const iMultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Divide dst by src including nghost ghost cells.
    // The two iMultiFabs MUST have the same underlying BoxArray.
    //
    static void Divide (iMultiFab&       dst,
			const iMultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);

    virtual void define (const BoxArray&            bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 int                        ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<IArrayBox>& factory) override;
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>()) override;
#endif

    const IArrayBox& operator[] (int K) const;

    IArrayBox& operator[] (int K);

    const IArrayBox& operator[] (const MFIter& mfi) const { 
	return this->FabArray<IArrayBox>::get(mfi); } 

    const IArrayBox& get (const MFIter& mfi) const { return operator[](mfi); }

    IArrayBox& operator[] (const MFIter& mfi) { 
	return this->FabArray<IArrayBox>::get(mfi); } 

    IArrayBox& get (const MFIter& mfi) { return operator[](mfi); }

    static void Initialize ();
    static void Finalize ();
};

}

#endif /*BL_IMULTIFAB_H*/
