//BL_COPYRIGHT_NOTICE

#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H

//
// $Id: MultiFab.H,v 1.29 1999-08-02 19:22:17 sstanley Exp $
//

#include <BLassert.H>
#include <FArrayBox.H>
#include <FabArray.H>

//
// Forward declaration.
//
typedef FabArrayId MultiFabId;

class MultiFabCopyDescriptor;

class MultiFab;
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
                                 bool                    extrap = false);

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
                                  bool                    extrap = false);


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
    public FabArray<Real, FArrayBox>
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
    //
    //@ManDoc: The (virtual) destructor.
    //
    virtual ~MultiFab();

    /*@ManDoc: This is a diagnostic function.  It writes (in ASCII) the state
               of the cell in the MultiFab with index location given by pnt
               to the output stream strm.
    */
    void probe (ostream& strm,
                IntVect& pt);

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

    friend void linInterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
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

    friend void linInterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
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

    /*@ManDoc: Copy on intersection within MultiFab.  Data is copied from
               valid regions to intersecting regions of definition.  The
               purpose is to fill in the boundary regions of each FAB in
               the MultiFab.
    */
    void FillBoundary ();

    /*@ManDoc: Same as FillBoundary(), but only copies
               num\_comp components starting at src\_comp.
    */
    void FillBoundary (int src_comp,
                       int num_comp);
    //
    //@ManDoc: The cached descriptor used by FillPeriodicBoundary().
    //
    MultiFabCopyDescriptor& theFPBmfcd (int   src_comp,
                                        int   num_comp,
                                        bool  no_overlap,
                                        bool  do_corners);
    //
    //@ManDoc: Flush the cache of self-intersection info used by FillBoundary.
    //
    static void FlushSICache ();
    //
    //@ManDoc: The size of the cache of self-intersection info.
    //
    static int SICacheSize ();

    /*@ManDoc: Copy from `src' to `dst' including `nghost' ghost cells.
               The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);
private:
    //
    // These are disabled.
    //
    MultiFab (const MultiFab& rhs);
    MultiFab& operator= (const MultiFab& rhs);
    //
    // Build and register the MFCD.  Guarantee the Register() MultiFabId == 0.
    //
    MultiFabCopyDescriptor& theFBmfcd (int src_comp,
                                       int num_comp);
    //
    // Some additional cached data.
    //
    MultiFabCopyDescriptor* m_FB_mfcd;    /* Hack for FillBoundary */
    int                     m_FB_scomp;   /* Hack for FillBoundary */
    int                     m_FB_ncomp;   /* Hack for FillBoundary */
    MultiFabCopyDescriptor* m_FPB_mfcd;   /* Hack for FillPeriodicBoundary */
    int                     m_FPB_scomp;  /* Hack for FillPeriodicBoundary */
    int                     m_FPB_ncomp;  /* Hack for FillPeriodicBoundary */
    bool                    m_FPB_noovlp; /* Hack for FillPeriodicBoundary */
    bool                    m_FPB_docrns; /* Hack for FillPeriodicBoundary */
};

class ConstMultiFabIterator;

//
//@Man:
//@Memo: An iterator over the FABs in a MultiFab.
/*@Doc:

  This class encapsulates an iterator over the FABs in a MultiFab such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the MultiFab that are actually owned by that CPU.
*/
class MultiFabIterator
    :
    public FabArrayIterator<Real,FArrayBox>
{
  public:
    //
    //@ManDoc: The constructor.
    //
    MultiFabIterator (MultiFab& multifab)
        :
        FabArrayIterator<Real, FArrayBox>(multifab) {}

  private:
    //
    // These are disallowed.
    //
    MultiFabIterator ();
    MultiFabIterator (const MultiFabIterator&);
    MultiFabIterator& operator= (const MultiFabIterator&);
};

//
//@Man:
//@Memo: An iterator over a dependent MultiFab.
/*@Doc:

  This class is used to iterator over a MultiFab (the dependent MultiFab) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a MultiFabIterator.  This allows us to iterate
  over a MultiFab and access the associated FAB(s) in one or more dependent
  MultiFabs.
*/

class DependentMultiFabIterator
    :
    public DependentFabArrayIterator<Real,FArrayBox>
{
  public:
    //
    //@ManDoc: Constructor over a mutable dependent MultiFab.
    //
    DependentMultiFabIterator (MultiFabIterator& controllermfiter,
                               MultiFab&         dependentmultifab)
        :
        DependentFabArrayIterator<Real,FArrayBox>(controllermfiter,
                                                  dependentmultifab) {}
    //
    //@ManDoc: Constructor over an immutable dependent MultiFab.
    //
    DependentMultiFabIterator (MultiFabIterator& controllermfiter,
                               const MultiFab&   dependentmultifab)
        :
        DependentFabArrayIterator<Real,FArrayBox>(controllermfiter,
                                                  dependentmultifab) {}
  private:
    //
    // These are disallowed.
    //
    DependentMultiFabIterator ();
    DependentMultiFabIterator (const DependentMultiFabIterator&);
    DependentMultiFabIterator& operator=(const DependentMultiFabIterator&);
};

//
//@Man:
//@Memo: A constant iterator over the FABs in a MultiFab.
/*@Doc:

  This class encapsulates an iterator over the FABs in a MultiFab such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the MultiFab that are actually owned by that CPU.  It is
  a constant iterator; you can access the FABs but you cannot change them.
*/

class ConstMultiFabIterator
    :
    public ConstFabArrayIterator<Real,FArrayBox>
{
  public:
    //
    //@ManDoc: The constructor.
    //
    ConstMultiFabIterator (const MultiFab& multifab)
        :
        ConstFabArrayIterator<Real,FArrayBox>(multifab) {}

  private:
    //
    // These are disallowed.
    //
    ConstMultiFabIterator ();
    ConstMultiFabIterator (const ConstMultiFabIterator&);
    ConstMultiFabIterator& operator= (const ConstMultiFabIterator&);
};

//
//@Man:
//@Memo: A constant iterator over a dependent MultiFab.
/*@Doc:

  This class is used to iterator over a MultiFab (the dependent MultiFab) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a ConstMultiFabIterator.  This allows us to
  iterate over a constant MultiFab and access the associated FAB(s) in one or
  more dependent MultiFabs.
*/

class ConstDependentMultiFabIterator
    :
    public ConstDependentFabArrayIterator<Real,FArrayBox>
{
  public:
    //
    //@ManDoc: The constructor.
    //
    ConstDependentMultiFabIterator (ConstMultiFabIterator& controllermfiter,
                                    const MultiFab&        dependentmultifab)
        :
        ConstDependentFabArrayIterator<Real,FArrayBox>(controllermfiter,
                                                       dependentmultifab) {}
  private:
    //
    // These are disallowed.
    //
    ConstDependentMultiFabIterator ();
    ConstDependentMultiFabIterator (const ConstDependentMultiFabIterator &);
    ConstDependentMultiFabIterator& operator= (const ConstDependentMultiFabIterator &);
};

class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<Real,FArrayBox>
{
  public:

    MultiFabCopyDescriptor ()
        :
        FabArrayCopyDescriptor<Real, FArrayBox>() {}

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

//
// Inlines.
//

inline
MultiFabCopyDescriptor&
MultiFab::theFBmfcd (int src_comp,
                     int num_comp)
{
    BL_ASSERT(num_comp > 0);
    BL_ASSERT(src_comp >= 0);

    if (!(m_FB_scomp == src_comp && m_FB_ncomp == num_comp))
    {
        m_FB_scomp = src_comp;
        m_FB_ncomp = num_comp;

        if (m_FB_mfcd == 0)
        {
            m_FB_mfcd = new MultiFabCopyDescriptor;
        }
        else
        {
            m_FB_mfcd->clear();
        }

        MultiFabId mfid = m_FB_mfcd->RegisterMultiFab(this);

        BL_ASSERT(mfid == MultiFabId(0));
    }

    return *m_FB_mfcd;
}

inline
MultiFabCopyDescriptor&
MultiFab::theFPBmfcd (int   src_comp,
                      int   num_comp,
                      bool  no_overlap,
                      bool  do_corners)
{
    BL_ASSERT(src_comp >= 0);
    BL_ASSERT(num_comp > 0);

    if (!(m_FPB_scomp  == src_comp   &&
          m_FPB_ncomp  == num_comp   &&
          m_FPB_noovlp == no_overlap &&
          m_FPB_docrns == do_corners))
    {
        m_FPB_scomp  = src_comp;
        m_FPB_ncomp  = num_comp;
        m_FPB_noovlp = no_overlap;
        m_FPB_docrns = do_corners;

        if (m_FPB_mfcd == 0)
        {
            m_FPB_mfcd = new MultiFabCopyDescriptor;
        }
        else
        {
            m_FPB_mfcd->clear();
        }

        MultiFabId mfid = m_FPB_mfcd->RegisterMultiFab(this);

        BL_ASSERT(mfid == MultiFabId(0));
    }

    return *m_FPB_mfcd;
}

inline
void
MultiFab::FillBoundary ()
{
    FillBoundary(0, n_comp);
}

inline
void
MultiFab::plus (Real val,
                int  nghost)
{
    plus(val,0,n_comp,nghost);
}

inline
void
MultiFab::plus (Real       val,
                const Box& region,
                int        nghost)
{
    plus(val,region,0,n_comp,nghost);
}

inline
void
MultiFab::mult (Real val,
                int  nghost)
{
    mult(val,0,n_comp,nghost);
}

inline
void
MultiFab::mult (Real       val,
                const Box& region,
                int        nghost)
{
    mult(val,region,0,n_comp,nghost);
}

inline
void
MultiFab::invert (Real numerator,
                  int  nghost)
{
    invert(numerator,0,n_comp,nghost);
}

inline
void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        nghost)
{
    invert(numerator,region,0,n_comp,nghost);
}

inline
void
MultiFab::negate (int nghost)
{
    negate(0,n_comp,nghost);
}

inline
void
MultiFab::negate (const Box& region,
                  int        nghost)
{
    negate(region,0,n_comp,nghost);
}

#endif /*BL_MULTIFAB_H*/
