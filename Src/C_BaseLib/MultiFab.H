//BL_COPYRIGHT_NOTICE

#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H

//
// $Id: MultiFab.H,v 1.19 1998-06-15 23:54:24 lijewski Exp $
//

#include <Assert.H>
#include <FArrayBox.H>
#include <FabArray.H>

#ifdef BL_USE_NEW_HFILES
#include <list>
using std::list;
#else
#include <list.h>
#endif

typedef FabArrayId MultiFabId;

class MultiFabCopyDescriptor;

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
                                 bool                    extrap = false);

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
                                  bool                    extrap = false);

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
    //@ManDoc: The cached descriptor used by FillPeriodicBoundary(). May be NULL.
    //
    MultiFabCopyDescriptor* getFPBmfcd (int  src_comp,
                                        int  num_comp,
                                        bool no_ovlp);
    //
    //@ManDoc: Set the cached copy descriptor used by FillPeriodicBoundary().
    //
    void setFPBmfcd (int                     src_comp,
                     int                     num_comp,
                     bool                    no_ovlp,
                     MultiFabCopyDescriptor* mfcd);
    //
    //@ManDoc: Flush the cache of SI information.
    //
    static void FlushSICache ();
    //
    //@ManDoc: The size of the SI cache.
    //
    static int SICacheSize ();

private:
    //
    // Holds single grid intersection record.
    //
    struct SIRec
    {
	SIRec ()
            :
            m_i(-1),
            m_j(-1) {}

	SIRec (int        i,
               int        j,
               const Box& bx)
	    :
            m_i(i),
            m_j(j),
            m_bx(bx)
        {
            assert(i >= 0);
            assert(j >= 0);
        }

        int       m_i;
	int       m_j;
	Box       m_bx;
	FillBoxId m_fbid;
    };
    //
    // Used in caching self-intersection info for FillBoundary().
    //
    struct SI
    {
        SI ()
            :
            m_scomp(-1),
            m_ncomp(-1),
            m_ngrow(-1) {}

        SI (const BoxArray& ba,
            int             scomp,
            int             ncomp,
            int             ngrow)
            :
            m_ba(ba),
            m_scomp(scomp),
            m_ncomp(ncomp),
            m_ngrow(ngrow)
        {
            assert(ncomp > 0);
            assert(scomp >= 0);
            assert(ngrow >= 0);
        }

        ~SI () {}

        bool operator== (const SI& rhs) const
        {
            return
                m_scomp == rhs.m_scomp &&
                m_ncomp == rhs.m_ncomp &&
                m_ngrow == rhs.m_ngrow &&
                m_ba == rhs.m_ba;
        }

        vector<SIRec> m_sirec;
        BoxArray      m_ba;
        int           m_scomp;
        int           m_ncomp;
        int           m_ngrow;
    };
    //
    // A useful typedef.
    //
    typedef list<SI> SIList;
    //
    // These are disabled.
    //
    MultiFab (const MultiFab& rhs);
    MultiFab& operator= (const MultiFab& rhs);
    //
    // FillBoundary() Helper functions.
    //
    vector<SIRec>& theFBsirec () const;

    vector<SIRec>& buildFBsirec (const SI& si) const;

    void buildFBmfcd (int src_comp,
                      int num_comp);

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
    //
    // Cache of grid self-intersection info used by FillBoundary().
    //
    static SIList m_SICache;
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
    assert(src_comp >= 0);
    assert(num_comp > 0);

    if (!(m_FB_scomp == src_comp && m_FB_ncomp == num_comp))
    {
        buildFBmfcd(src_comp, num_comp);
    }

    return *m_FB_mfcd;
}

inline
MultiFabCopyDescriptor*
MultiFab::getFPBmfcd (int  src_comp,
                      int  num_comp,
                      bool no_ovlp)
{
    assert(src_comp >= 0);
    assert(num_comp > 0);

    if (!(m_FPB_scomp == src_comp &&
          m_FPB_ncomp == num_comp &&
          m_FPB_noovlp == no_ovlp))
    {
        //
        // Null it out so it'll be rebuilt.
        //
        delete m_FPB_mfcd;

        m_FPB_mfcd   = 0;
        m_FPB_scomp  = -1;
        m_FPB_ncomp  = -1;
        m_FPB_noovlp = false;
    }

    return m_FPB_mfcd;
}

inline
void
MultiFab::setFPBmfcd (int                     src_comp,
                      int                     num_comp,
                      bool                    no_ovlp,
                      MultiFabCopyDescriptor* mfcd)
{
    assert(num_comp > 0);
    assert(!(mfcd == 0));
    assert(src_comp >= 0);
    
    delete m_FPB_mfcd;

    m_FPB_mfcd   = mfcd;
    m_FPB_scomp  = src_comp;
    m_FPB_ncomp  = num_comp;
    m_FPB_noovlp = no_ovlp;
}

inline
int
MultiFab::SICacheSize ()
{
    return m_SICache.size();
}

inline
vector<MultiFab::SIRec>&
MultiFab::theFBsirec () const
{
    //
    // We assume `m_FB_scomp' & `m_FB_ncomp' have already been set.
    //
    assert(m_FB_ncomp > 0);
    assert(m_FB_scomp >= 0);

    const SI si(boxArray(), m_FB_scomp, m_FB_ncomp, nGrow());
    
    for (SIList::iterator it = m_SICache.begin(); it != m_SICache.end(); ++it)
    {
        if (*it == si)
        {
            return (*it).m_sirec;
        }
    }

    return buildFBsirec(si);
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
