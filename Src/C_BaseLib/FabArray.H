#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H
//
// $Id: FabArray.H,v 1.193 2008-09-30 17:17:04 lijewski Exp $
//
#include <iostream>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <BLassert.H>
#include <PArray.H>
#include <Array.H>

#include <Box.H>
#include <BoxLib.H>
#include <BoxArray.H>
#include <BoxDomain.H> 
#include <FArrayBox.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <ccse-mpi.H>

#include <Profiler.H>

//
// Helper class
//

class FillBoxId
{
  public:

    FillBoxId ()
        :
        m_fillBoxId(-1),
        m_fabIndex(-1)
        {}
    FillBoxId (int newid, const Box& fillbox)
        :
        m_fillBox(fillbox),
        m_fillBoxId(newid),
        m_fabIndex(-1)
        {}

    int Id () const              { return m_fillBoxId;    }
    int FabIndex () const        { return m_fabIndex;     }
    void FabIndex (int fabindex) { m_fabIndex = fabindex; }
    const Box& box () const      { return m_fillBox;      }

private:

    Box m_fillBox;
    int m_fillBoxId;
    int m_fabIndex;
};

//
// Used to cache some CommData stuff in CollectData().
//

class CommDataCache
{
public:

    CommDataCache ();

    void operator= (const Array<CommData>& rhs);

    bool isValid () const { return m_valid; }

    Array<CommData>& theCommData () { return m_commdata; }

    const Array<CommData>& theCommData () const { return m_commdata; }

private:

    Array<CommData> m_commdata;
    bool            m_valid;
};

class FabArrayBase
{
public:

    FabArrayBase ();
    FabArrayBase (const BoxArray& bx, int nvar, int ngrow);
    FabArrayBase (const BoxArray& bx, int nvar, int ngrow, const DistributionMapping& map);

    virtual ~FabArrayBase();
    //
    // Returns the grow factor that defines the region of definition.
    //
    int nGrow () const;
    //
    // Returns number of variables associated with each point (nvar).
    //
    int nComp () const;
    //
    // Returns a constant reference to the BoxArray that defines the
    // valid region associated with this FabArray.
    //
    const BoxArray& boxArray () const;
    //
    // Returns a constant reference to the Kth Box in the BoxArray.
    // That is, the valid region of the Kth grid.
    //
    virtual const Box& box (int K) const;
    //
    // Returns the Kth FABs Box in the FabArray.
    // That is, the region the Kth fab is actually defined on.
    //
    virtual Box fabbox (int K) const;
    //
    // Returns the number of FABs in the FabArray..
    //
    int size () const;
    //
    // Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const;
    //
    // Flush the cache of self-intersection info used by FillBoundary.
    //
    static void FlushSICache ();
    //
    // The size of the cache of self-intersection info.
    //
    static int SICacheSize ();
    //
    // Holds single grid intersection record.
    //
    struct SIRec
    {
        SIRec ()
            :
            m_i(-1),
            m_j(-1)
            {}

        SIRec (int        i,
               int        j,
               const Box& bx)
            :
            m_i(i),
            m_j(j),
            m_bx(bx)
            {
                BL_ASSERT(i >= 0);
                BL_ASSERT(j >= 0);
            }

        SIRec (const SIRec& rhs)
            :
            m_i(rhs.m_i),
            m_j(rhs.m_j),
            m_bx(rhs.m_bx),
            m_fbid(rhs.m_fbid) {}

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
        SI ();

        SI (const BoxArray&            ba,
            const DistributionMapping& dm,
            int                        ngrow);

        SI (const SI& rhs);

        ~SI ();

        bool operator== (const SI& rhs) const;
        bool operator!= (const SI& rhs) const;

        Array<int>          m_cache;    // Snds cached for CollectData().
        CommDataCache       m_commdata; // Yet another cache for CollectData().
        std::vector<SIRec>  m_sirec;
        BoxArray            m_ba;
        DistributionMapping m_dm;
        int                 m_ngrow;
        bool                m_reused;
    };
    //
    // Builds self-intersection records.
    //
    static SI& BuildFBsirec (const SI& si, const FabArrayBase& mf);
    //
    // Returns cached self-intersection records or builds them.
    //
    static SI& TheFBsirec (int scomp, int ncomp, const FabArrayBase& mf);

protected:
    //
    // The data ...
    //
    mutable BoxArray    boxarray;   /* So FabSets can modify'm */
    DistributionMapping distributionMap;
    int                 n_grow;
    int                 n_comp;
};

class MFIter
{
public:
    //
    //@ManDoc: Construct a MFIter.
    //
    explicit MFIter (const FabArrayBase& fabarray);
    //
    //@ManDoc: Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    //@ManDoc: Returns the Box of the FAB at which we currently point.
    //
    Box fabbox () const;
    //
    //@ManDoc: Increments iterator to the next FAB in the FabArray that we own.
    //
    void operator++ ();
    //
    //@ManDoc: Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid ();
    //
    //@ManDoc: The index into the underlying BoxArray of the current FAB.
    //
    int index () const;
    //
    //@ManDoc: Constant reference to FabArray over which we're iterating.
    //
    const FabArrayBase& theFabArrayBase () const;

protected:

    const FabArrayBase& fabArray;
    int                 currentIndex;
};

//
// A forward declaration.
//
template <class FAB> class FabArray;
template <class FAB> class FabArrayCopyDescriptor;

/*
  A Collection of Fortran Array-like Objects


  The FabArray<FAB> class implements a collection (stored as an array) of
  Fortran array-like objects.  The parameterized type FAB is intended to be
  any class derived from BaseFab<T>.  For example, FAB may be a BaseFab of
  integers, so we could write:

    FabArray<BaseFab<int> > int_fabs;

  Then int_fabs is a FabArray that can hold a collection of BaseFab<int>
  objects.

  FabArray is not just a general container class for Fortran arrays.  It is
  intended to hold "grid" data for use in finite difference calculations in
  which the data is defined on a union of (usually disjoint) rectangular
  regions embedded in a uniform index space.  This region, called the valid
  region, is represented by a BoxArray.  For the purposes of this discussion,
  the Kth Box in the BoxArray represents the interior region of the Kth grid.

  Since the intent is to be used with finite difference calculations a
  FabArray also includes the notion of a boundary region for each grid.  The
  boundary region is specified by the ngrow parameter which tells the FabArray
  to allocate each FAB to be ngrow cells larger in all directions than the
  underlying Box.  The larger region covered by the union of all the FABs is
  called the region of definition.  The underlying notion is that the valid
  region contains the grid interior data and the region of definition includes
  the interior region plus the boundary areas.

  Operations are available to copy data from the valid regions into these
  boundary areas where the two overlap.  The number of components, that is,
  the number of values that can be stored in each cell of a FAB, is either
  given as an argument to the constructor or is inherent in the definition of
  the underlying FAB.  Each FAB in the FabArray will have the same number of
  components.

  In summary, a FabArray is an array of FABs.  The Kth element contains a FAB
  that holds the data for the Kth grid, a Box that defines the valid region
  of the Kth grid.

  A typical use for a FabArray would be to hold the solution vector or
  right-hand-side when solving a linear system of equations on a union of
  rectangular grids.  The copy operations would be used to copy data from the
  valid regions of neighboring grids into the boundary regions after each
  relaxation step of the iterative method.  If a multigrid method is used, a
  FabArray could be used to hold the data at each level in the multigrid
  hierarchy.

  This class is a concrete class not a polymorphic one.

  This class does NOT provide a copy constructor or assignment operator.
*/

//
// An enumumeration that controls whether or not the memory for a FAB
// will actually be allocated on construction of a FabArray.
// Possible values are: Fab\_noallocate and Fab\_allocate.
//

enum FabAlloc { Fab_noallocate = 0, Fab_allocate };

template <class FAB>
class FabArray
    :
    public FabArrayBase
{
public:

    typedef typename FAB::value_type value_type;
    //
    //@ManDoc: Constructs an empty FabArray<FAB>.
    //
    FabArray ();

    /*@ManDoc: Construct a FabArray<FAB> with a valid region defined by bxs
               and a region of definition defined by the grow factor ngrow
               and the number of components nvar.
               If mem\_mode is defined to be Fab\_allocate then FABs are
               allocated for each Box in the BoxArray.  The size of the Kth
               FAB is given by bxs[K] grown by ngrow.  If mem\_mode is defined
               to be Fab\_noallocate, then no FABs are allocated at this time,
               but can be defined later.  The number of components in each
               FAB is not specified and is expected to be implicit in the
               definition of the FAB class.  That is, the FAB constructor will
               take only a Box argument.  Call this constructor number two.
    */
    FabArray (const BoxArray& bxs,
              int             nvar,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);
    //
    //@ManDoc: The (virtual) destructor -- deletes all FABs in the array.
    //
    virtual ~FabArray ();

    /*@ManDoc: Define this FabArray identically to that performed by
               the constructor having an analogous function signature.
               This is only valid if this FabArray was defined using
               the default constructor.
    */
    void define (const BoxArray& bxs,
                 int             nvar,
                 int             ngrow,
                 FabAlloc        mem_mode);

    void define (const BoxArray&            bxs,
                 int                        nvar,
                 int                        ngrow,
		 const DistributionMapping& dm,
                 FabAlloc                   mem_mode);

    /*@ManDoc: Returns true if the FabArray is well-defined.  That is,
               if FABs are allocated for each Box in the BoxArray and the
               sizes of the FABs and the number of components are consistent
               with the definition of the FabArray.
    */
    bool ok () const;

    /*@ManDoc: Returns a constant reference to the FAB associated with the
               Kth element.
    */
    const FAB& operator[] (const MFIter& mfi) const;

    const FAB& get (const MFIter& mfi) const;
    //
    //@ManDoc: Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (const MFIter& mfi);

    FAB& get (const MFIter& mfi);

    /*@ManDoc: Returns a constant reference to the FAB associated with the
               Kth element.
    */
    const FAB& operator[] (int K) const;
    //
    //@ManDoc: Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (int K);
    //
    //@ManDoc: Explicitly set the Kth FAB in the FabArray to point to elem.
    //
    void setFab (int K, FAB* elem);
    //
    //@ManDoc: Releases FAB memory in the FabArray.
    //
    void clear ();
    //
    //@ManDoc: Set all components in the valid region of each FAB to val.
    //
    void setVal (value_type val);
    void operator= (const value_type& val);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp to val.
               Also set the value of nghost boundary cells.
    */
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp, as well
               as nghost boundary cells, to val, provided they also intersect
               with the Box region.
    */
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells.
    */
    void setVal (value_type val,
                 int        nghost);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells, that also
               intersect the Box region.
    */
    void setVal (value_type val,
                 const Box& region,
                 int        nghost);
    //
    //@ManDoc: Set all values in the boundary region to val.
    //
    void setBndry (value_type val);

    /*@ManDoc: Set ncomp values in the boundary region, starting at
               start\_comp to val.
    */
    void setBndry (value_type val,
                   int        strt_comp,
                   int        ncomp);

    /*@ManDoc: This function copies data from fa to this FabArray.  Each FAB
               in fa is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the valid region of each FAB.
    */
    void copy (const FabArray<FAB>& fa);

    /*@ManDoc: This function copies data from src to this FabArray.  Each FAB
               in src is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the num\_comp components starting at src\_comp
               in the FabArray src, with the destination components in this
               FabArray starting at dest\_comp.
    */
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp);

    /*@ManDoc: Copies the values contained in the intersection of the
               valid region of this FabArray with the FAB dest into dest.
    */
    void copy (FAB& dest) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               valid region of this FabArray with the FAB dest and the Box
               subbox into that subregion of dest.
    */
    void copy (FAB&       dest,
               const Box& subbox) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               num\_comp component valid region of this FabArray, starting at
               component src\_comp, with the FAB dest into dest, starting at
               component dest\_comp in dest.
    */
    void copy (FAB& dest,
               int  src_comp,
               int  dest_comp,
               int  num_comp) const;

    /*@ManDoc: Copies the values contained in the intersection of the
               num\_comp component valid region of this FabArray, starting at
               component src\_comp, with the FAB dest and the Box subbox, into
               dest, starting at component dest\_comp in dest.
    */
    void copy (FAB&       dest,
               const Box& subbox,
               int        src_comp,
               int        dest_comp,
               int        num_comp) const;

    void shift (const IntVect& v);

    bool defined (int i) const;
    //
    // Copy on intersection within a FabArray.  Data is copied from
    // valid regions to intersecting regions of definition.  The
    // purpose is to fill in the boundary regions of each FAB in
    // the FabArray.
    //
    void FillBoundary ();
    //
    // Same as FillBoundary(), but only copies ncomp components starting at scomp.
    //
    void FillBoundary (int scomp, int ncomp);

protected:

    PArray<FAB> fabparray;

private:
    //
    // These are disallowed.
    //
    FabArray (const FabArray<FAB>&);
    FabArray<FAB>& operator= (const FabArray<FAB>&);
    //
    // This is used locally in all define functions.
    //
    void AllocFabs ();
};

class FabArrayId
{
public:

    explicit FabArrayId (int newid = -1)
        :
        fabArrayId(newid) {}

    int Id () const { return fabArrayId; }

    bool operator== (const FabArrayId& rhs) const
    {
        return fabArrayId == rhs.fabArrayId;
    }

private:

    int fabArrayId;
};

//
// This enum and the FabCopyDescriptor class should really be nested
// in FabArrayCopyDescriptor (not done for portability reasons).
//

enum FillType { FillLocally, FillRemotely, Unfillable };

template <class FAB>
struct FabCopyDescriptor
{
    FabCopyDescriptor ();

    ~FabCopyDescriptor ();

    FAB*     localFabSource;
    Box      subBox;
    int      myProc;
    int      copyFromProc;
    int      copyFromIndex;
    int      fillBoxId;
    int      srcComp;
    int      destComp;
    int      nComp;
    FillType fillType;
    bool     cacheDataAllocated;

private:
    //
    // Disallowed.
    //
    FabCopyDescriptor (const FabCopyDescriptor&);
    FabCopyDescriptor& operator= (const FabCopyDescriptor&);
};

template <class FAB>
FabCopyDescriptor<FAB>::FabCopyDescriptor ()
    :
    localFabSource(0),
    myProc(-1),
    copyFromProc(-1),
    copyFromIndex(-1),
    fillBoxId(-1),
    srcComp(-1),
    destComp(-1),
    nComp(-1),
    fillType(Unfillable),
    cacheDataAllocated(false)
{}

template <class FAB>
FabCopyDescriptor<FAB>::~FabCopyDescriptor ()
{
    if (cacheDataAllocated)
        delete localFabSource;
}

//
// This class orchestrates filling a destination fab of size destFabBox
// from fabarray on the local processor (myProc).
//

template <class FAB>
class FabArrayCopyDescriptor
{
  typedef std::multimap<int,FabCopyDescriptor<FAB>*> FCDMap;
  typedef typename FCDMap::value_type                FCDMapValueType;
  typedef typename FCDMap::iterator                  FCDMapIter;
  typedef typename FCDMap::const_iterator            FCDMapConstIter;

  public:

    FabArrayCopyDescriptor ();

    ~FabArrayCopyDescriptor ();

    FabArrayId RegisterFabArray(FabArray<FAB> *fabarray);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp);
    //
    // Add a box but only from FabArray[fabarrayindex].
    //
    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        fabarrayindex,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp,
                      bool       bUseValidBox = true);

    void CollectData (Array<int>*    snd_cache = 0,
                      CommDataCache* cd_cache  = 0);

    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab);

    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab,
                  const Box&       destBox);

    void PrintStats () const;

    bool DataAvailable () const { return dataAvailable; }

    void clear ();

    int nFabArrays () const { return fabArrays.size(); }

    int nFabComTags () const { return fabComTagList.size(); }

    int nFabCopyDescs () const { return fabCopyDescList.size(); }

  protected:
    //
    // Helper function for AddBox() routines.
    //
    void AddBoxDoIt (FabArrayId fabarrayid,
                     const Box& destFabBox,
                     BoxList*   returnedUnfilledBoxes,
                     int        faindex,
                     int        srccomp,
                     int        destcomp,
                     int        numcomp,
                     bool       bUseValidBox,
                     BoxDomain& unfilledBoxDomain);

    std::vector<FabArray<FAB>*> fabArrays;
    std::vector<FCDMap>           fabCopyDescList;
    std::vector<FabComTag>        fabComTagList;
    int                           nextFillBoxId;
    bool                          dataAvailable;

  private:
    //
    // These are disallowed.
    //
    FabArrayCopyDescriptor (const FabArrayCopyDescriptor<FAB>&);
    FabArrayCopyDescriptor<FAB>& operator= (const FabArrayCopyDescriptor<FAB> &);
};


inline
int
FabArrayBase::nGrow () const
{
    return n_grow;
}

inline
const BoxArray&
FabArrayBase::boxArray () const
{
    return boxarray;
}

inline
const Box&
FabArrayBase::box (int K) const
{
    return boxarray[K];
}

inline
int
FabArrayBase::size () const
{
    return boxarray.size();
}

inline
int
FabArrayBase::nComp () const
{
    return n_comp;
}

inline
const DistributionMapping&
FabArrayBase::DistributionMap () const
{
    return distributionMap;
}

inline
const FabArrayBase&
MFIter::theFabArrayBase () const
{
    return fabArray;
}

inline
int
MFIter::index () const
{
    return currentIndex;
}

template <class FAB>
inline
bool
FabArray<FAB>::defined (int i) const
{
    return fabparray.defined(i);
}

template <class FAB>
inline
const FAB&
FabArray<FAB>::operator[] (const MFIter& mfi) const
{
    return fabparray[mfi.index()];
}

template <class FAB>
inline
const FAB&
FabArray<FAB>::get (const MFIter& mfi) const
{
    return fabparray[mfi.index()];
}

template <class FAB>
inline
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    return fabparray[mfi.index()];
}

template <class FAB>
inline
FAB&
FabArray<FAB>::get (const MFIter& mfi)
{
    return fabparray[mfi.index()];
}

template <class FAB>
inline
const FAB&
FabArray<FAB>::operator[] (int K) const
{
    return fabparray[K];
}

template <class FAB>
inline
FAB&
FabArray<FAB>::operator[] (int K)
{
    return fabparray[K];
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    fabparray.clear();
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       int        nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type   val,
                         const Box& region,
                         int        nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

template <class FAB>
FabArray<FAB>::FabArray ()
    :
    fabparray(0,PArrayManage)
{}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc)
    :
    fabparray(0, PArrayManage)
{
    define(bxs,nvar,ngrow,alloc);
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{}

template <class FAB>
bool
FabArray<FAB>::ok () const
{
    long isok = true;

    for (MFIter fai(*this); fai.isValid() && isok; ++fai)
    {
        if (fabparray.defined(fai.index()))
        {
            if (get(fai).box() != BoxLib::grow(box(fai.index()),n_grow))
            {
                isok = false;
            }
        }
        else
        {
            isok = false;
        }
    }

    ParallelDescriptor::ReduceLongAnd(isok);

    return isok != 0;
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray& bxs,
                       int             nvar,
                       int             ngrow,
                       FabAlloc        alloc)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap.define(boxarray,ParallelDescriptor::NProcs());
    fabparray.resize(bxs.size());
    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       int                        nvar,
                       int                        ngrow,
                       const DistributionMapping& dm,
                       FabAlloc                   alloc)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);
    BL_ASSERT(dm.ProcessorMap().size() == bxs.size()+1);
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap = dm;
    fabparray.resize(bxs.size());
    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class FAB>
void
FabArray<FAB>::AllocFabs ()
{
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Box tmp = BoxLib::grow(fai.validbox(), n_grow);

        fabparray.set(fai.index(), new FAB(tmp, n_comp));
    }
}

template <class FAB>
void
FabArray<FAB>::setFab (int  boxno,
                       FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();

    BL_ASSERT(n_comp == elem->nComp());
    BL_ASSERT(boxarray.size() > 0);
    BL_ASSERT(elem->box() == BoxLib::grow(boxarray[boxno],n_grow));
    BL_ASSERT(!fabparray.defined(boxno));
    BL_ASSERT(distributionMap[boxno] == ParallelDescriptor::MyProc());

    fabparray.set(boxno,elem);
}

template <class FAB>
void
FabArray<FAB>::setBndry (value_type val)
{
    setBndry(val, 0, n_comp);
}

template <class FAB>
void
FabArray<FAB>::setBndry (value_type val,
                         int        strt_comp,
                         int        ncomp)
{
    if (n_grow > 0)
    {
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).setComplement(val, fai.validbox(), strt_comp, ncomp);
        }
    }
}

//
// Used by FabArray::copy().
//
struct CopyComTag
{
    CopyComTag ();
    CopyComTag (const CopyComTag& cct);

    Box box;
    int fabIndex;
    int srcIndex;
};

struct CPC
{
    CPC ();

    CPC (const BoxArray&            dstba,
         const BoxArray&            srcba,
         const DistributionMapping& dstdm,
         const DistributionMapping& srcdm);

    CPC (const CPC& rhs);

    ~CPC ();

    bool operator== (const CPC& rhs) const;
    bool operator!= (const CPC& rhs) const;

    static CPC& TheCPC (const CPC& cpc, bool& got_from_cache);

    static void FlushCache ();

    BoxArray            m_dstba;
    BoxArray            m_srcba;
    DistributionMapping m_dstdm;
    DistributionMapping m_srcdm;
    //
    // The cache of local and send/recv info per FabArray::copy().
    //
    std::vector<CopyComTag>          m_LocalTags;
    Array< std::vector<CopyComTag> > m_SndTags;
    Array< std::vector<CopyComTag> > m_RcvTags;
    //
    // Has this CPC ever been reused in the cache?
    //
    bool m_reused;
};

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::copy()");

    if (boxarray == src.boxarray && distributionMap == src.distributionMap)
    {
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).copy(src[fai],fai.validbox(),scomp,fai.validbox(),dcomp,ncomp);
        }
        return;
    }

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    CPC cpc(boxarray, src.boxarray, distributionMap, src.distributionMap);

    bool got_from_cache = false;

    CPC& thecpc = CPC::TheCPC(cpc, got_from_cache);

    if (got_from_cache)
    {
        //
        // All we need to do here are the local copy()s.
        //
        for (int j = 0; j < thecpc.m_LocalTags.size(); j++)
        {
            CopyComTag& tag = thecpc.m_LocalTags[j];

            BL_ASSERT(distributionMap[tag.fabIndex] == MyProc);
            BL_ASSERT(src.distributionMap[tag.srcIndex] == MyProc);
            
            fabparray[tag.fabIndex].copy(src[tag.srcIndex],tag.box,scomp,tag.box,dcomp,ncomp);
        }
    }
    else
    {
        CopyComTag tag;

        thecpc.m_SndTags.resize(NProcs);
        thecpc.m_RcvTags.resize(NProcs);

        for (int i = 0; i < size(); i++)
        {
            std::vector< std::pair<int,Box> > isects = src.boxarray.intersections(boxarray[i]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                const Box& bx  = isects[ii].second;
                const int  iii = isects[ii].first;

                if (distributionMap[i] == MyProc)
                {
                    if (src.distributionMap[iii] == MyProc)
                    {
                        fabparray[i].copy(src[iii],bx,scomp,bx,dcomp,ncomp);

                        tag.box      = bx;
                        tag.fabIndex = i;
                        tag.srcIndex = iii;

                        thecpc.m_LocalTags.push_back(tag);
                    }
                    else
                    {
                        tag.box      = bx;
                        tag.fabIndex = i;
                        thecpc.m_RcvTags[src.distributionMap[iii]].push_back(tag);
                    }
                }
                else
                {
                    if (src.distributionMap[iii] == MyProc)
                    {
                        tag.box      = bx;
                        tag.fabIndex = iii;
                        thecpc.m_SndTags[distributionMap[i]].push_back(tag);
                    }
                }
            }
        }
    }
    
    if (NProcs == 1) return;

#ifdef BL_USE_MPI

    static bool do_alltoallv     = false;
    static bool isdef_alltoallv  = false;

    if (!isdef_alltoallv)
    {
        isdef_alltoallv = true;
        ParmParse pp("fabarray");
        pp.query("do_alltoallv",do_alltoallv);
    }

    if (do_alltoallv)
    {
        size_t Nsend = 0, Nrecv = 0;

        Array<int> sendcount(NProcs,0), recvcount(NProcs,0);
        Array<int> senddspls(NProcs,0), recvdspls(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            int Ns = 0;
            for (int j = 0; j < thecpc.m_SndTags[i].size(); j++)
                Ns += thecpc.m_SndTags[i][j].box.numPts() * ncomp;
            sendcount[i] = Ns;
            Nsend       += Ns;

            int Nr = 0;
            for (int j = 0; j < thecpc.m_RcvTags[i].size(); j++)
                Nr += thecpc.m_RcvTags[i][j].box.numPts() * ncomp;
            recvcount[i] = Nr;
            Nrecv       += Nr;
        }

        for (int i = 1; i < NProcs; i++)
        {
            senddspls[i] = senddspls[i-1] + sendcount[i-1];
            recvdspls[i] = recvdspls[i-1] + recvcount[i-1];
        }

        BL_ASSERT((Nsend*sizeof(value_type)) < std::numeric_limits<size_t>::max());
        BL_ASSERT((Nrecv*sizeof(value_type)) < std::numeric_limits<size_t>::max());

        value_type* senddata = static_cast<value_type*>(BoxLib::The_Arena()->alloc(Nsend*sizeof(value_type)));
        value_type* recvdata = static_cast<value_type*>(BoxLib::The_Arena()->alloc(Nrecv*sizeof(value_type)));

        FAB fab;

        value_type* senddptr = senddata;

        for (int i = 0; i < NProcs; i++)
        {
            for (int j = 0; j < thecpc.m_SndTags[i].size(); j++)
            {
                const Box& bx = thecpc.m_SndTags[i][j].box;
                fab.resize(bx, ncomp);
                fab.copy(src[thecpc.m_SndTags[i][j].fabIndex],bx,scomp,bx,0,ncomp);
                const int N = bx.numPts() * ncomp;
                memcpy(senddptr, fab.dataPtr(), N*sizeof(value_type));
                senddptr += N;
            }
        }

        BL_ASSERT(senddata+Nsend == senddptr);

        BL_MPI_REQUIRE( MPI_Alltoallv(senddata,
                                      sendcount.dataPtr(),
                                      senddspls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<value_type>::type(),
                                      recvdata,
                                      recvcount.dataPtr(),
                                      recvdspls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<value_type>::type(),
                                      ParallelDescriptor::Communicator()) );

        BoxLib::The_Arena()->free(senddata);

        value_type* recvdptr = recvdata;

        for (int i = 0; i < NProcs; i++)
        {
            for (int j = 0; j < thecpc.m_RcvTags[i].size(); j++)
            {
                const Box& bx = thecpc.m_RcvTags[i][j].box;
                fab.resize(bx, ncomp);
                const int N = bx.numPts() * ncomp;
                memcpy(fab.dataPtr(), recvdptr, N*sizeof(value_type));
                fabparray[thecpc.m_RcvTags[i][j].fabIndex].copy(fab,bx,0,bx,dcomp,ncomp);
                recvdptr += N;
            }
        }

        BL_ASSERT(recvdata+Nrecv == recvdptr);

        BoxLib::The_Arena()->free(recvdata);
    }
    else
    {
        Array<value_type*> fab_data(NProcs);
        Array<int>         indx(NProcs);
        Array<MPI_Status>  status(NProcs);
        Array<MPI_Request> reqs(NProcs);

        const int seqno = ParallelDescriptor::SeqNum();

        int NWaits = 0;
        //
        // Post one receive for each chunk being sent by other CPUs.
        //
        for (int i = 0; i < NProcs; i++)
        {
            reqs[i] = MPI_REQUEST_NULL;

            if (!thecpc.m_RcvTags[i].empty())
            {
                NWaits++;
                size_t N = 0;
                for (unsigned int j = 0; j < thecpc.m_RcvTags[i].size(); j++)
                    N += thecpc.m_RcvTags[i][j].box.numPts() * ncomp;
                BL_ASSERT(N < std::numeric_limits<size_t>::max());
                fab_data[i] = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
                reqs[i] = ParallelDescriptor::Arecv(fab_data[i],N,i,seqno).req();
            }
        }

        FAB fab;
        //
        // Send the agglomerated FAB data.
        //
        for (int i = 0; i < NProcs; i++)
        {
            if (!thecpc.m_SndTags[i].empty())
            {
                size_t N = 0;
                for (unsigned int j = 0; j < thecpc.m_SndTags[i].size(); j++)
                    N += thecpc.m_SndTags[i][j].box.numPts() * ncomp;
                BL_ASSERT(N < std::numeric_limits<size_t>::max());
                value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
                value_type* dptr = data;
                for (unsigned int j = 0; j < thecpc.m_SndTags[i].size(); j++)
                {
                    const Box& bx = thecpc.m_SndTags[i][j].box;
                    fab.resize(bx, ncomp);
                    fab.copy(src[thecpc.m_SndTags[i][j].fabIndex],bx,scomp,bx,0,ncomp);
                    int count = bx.numPts() * ncomp;
                    memcpy(dptr, fab.dataPtr(), count*sizeof(value_type));
                    dptr += count;
                }
                BL_ASSERT(data+N == dptr);
                ParallelDescriptor::Send(data, N, i, seqno);
                BoxLib::The_Arena()->free(data);
            }
        }
        //
        // Now receive and unpack FAB data.
        //
        for (int completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(reqs, completed, indx, status);

            for (int k = 0; k < completed; k++)
            {
                value_type* dptr = fab_data[indx[k]];
                BL_ASSERT(!(dptr == 0));
                for (unsigned int j = 0; j < thecpc.m_RcvTags[indx[k]].size(); j++)
                {
                    const Box& bx = thecpc.m_RcvTags[indx[k]][j].box;
                    fab.resize(bx, ncomp);
                    int N = bx.numPts() * ncomp;
                    memcpy(fab.dataPtr(), dptr, N*sizeof(value_type));
                    fabparray[thecpc.m_RcvTags[indx[k]][j].fabIndex].copy(fab,bx,0,bx,dcomp,ncomp);
                    dptr += N;
                }
                BoxLib::The_Arena()->free(fab_data[indx[k]]);
            }
        }
    }

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src)
{
    copy(src,0,0,nComp());
}

//
// Copies to FABs, note that destination is first arg.
//

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest) const
{
    copy(dest, dest.box(), 0, 0, dest.nComp());
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
                     const Box& subbox) const
{
    copy(dest, subbox, 0, 0, dest.nComp());
}

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest,
                     int  scomp,
                     int  dcomp,
                     int  ncomp) const
{
    copy(dest, dest.box(), scomp, dcomp, ncomp);
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
                     const Box& subbox,
                     int        scomp,
                     int        dcomp,
                     int        ncomp) const
{
    BL_ASSERT(dcomp + ncomp <= dest.nComp());

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::copy(Fab&, ...)");

    if (ParallelDescriptor::NProcs() == 1)
    {
        for (int j = 0; j < size(); ++j)
        {
            if (boxarray[j].intersects(subbox))
            {
                Box destbox = boxarray[j] & subbox;

                dest.copy(fabparray[j],destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    FArrayBox ovlp;

    for (int i = 0; i < size(); i++)
    {
        if (subbox.intersects(boxarray[i]))
        {
            Box bx = subbox & boxarray[i];

            ovlp.resize(bx,ncomp);

            if (ParallelDescriptor::MyProc() == distributionMap[i])
            {
                ovlp.copy(fabparray[i],bx,scomp,bx,0,ncomp);
            }

            const int N = bx.numPts()*ncomp;

            ParallelDescriptor::Bcast(ovlp.dataPtr(),N,distributionMap[i]);

            dest.copy(ovlp,bx,0,bx,dcomp,ncomp);
        }
    }
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val)
{
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).setVal(val);
    }
}

template <class FAB>
inline
void
FabArray<FAB>::operator= (const value_type& val)
{
    setVal(val);
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).setVal(val,BoxLib::grow(fai.validbox(),nghost), comp, ncomp);
    }
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Box b = BoxLib::grow(fai.validbox(),nghost) & region;

        if (b.ok())
            get(fai).setVal(val, b, comp, ncomp);
    }
}


template <class FAB>
void
FabArray<FAB>::shift (const IntVect& v)
{
    for(int id(0); id < BL_SPACEDIM; ++id)
    {
      boxarray.shift(id, v[id]);
    }
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).shift(v);
    }
}


template <class FAB>
FabArrayCopyDescriptor<FAB>::FabArrayCopyDescriptor ()
    :
    nextFillBoxId(0),
    dataAvailable(false)
{}

template <class FAB>
FabArrayId
FabArrayCopyDescriptor<FAB>::RegisterFabArray(FabArray<FAB>* fabarray)
{
    BL_ASSERT(fabArrays.size() == fabCopyDescList.size());

    fabArrays.push_back(fabarray);  /* Bump size() by one */

    fabCopyDescList.resize(fabArrays.size(), FCDMap());

    return FabArrayId(fabArrays.size() - 1);
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::AddBoxDoIt (FabArrayId fabarrayid,
                                         const Box& destFabBox,
                                         BoxList*   returnedUnfilledBoxes,
                                         int        faindex,
                                         int        srccomp,
                                         int        destcomp,
                                         int        numcomp,
                                         bool       bUseValidBox,
                                         BoxDomain& unfilledBoxDomain)
{
    const int MyProc = ParallelDescriptor::MyProc();

    FabArray<FAB>* fabArray = fabArrays[fabarrayid.Id()];

    BL_ASSERT(faindex >= 0 && faindex < fabArray->size());

    Box intersect = destFabBox;

    if (bUseValidBox)
    {
        intersect &= fabArray->box(faindex);
    }
    else
    {
        intersect &= fabArray->fabbox(faindex);
    }

    if (intersect.ok())
    {
        FabCopyDescriptor<FAB>* fcd = new FabCopyDescriptor<FAB>;

        int remoteProc     = fabArray->DistributionMap()[faindex];
        fcd->fillBoxId     = nextFillBoxId;
        fcd->subBox        = intersect;
        fcd->myProc        = MyProc;
        fcd->copyFromProc  = remoteProc;
        fcd->copyFromIndex = faindex;
        fcd->srcComp       = srccomp;
        fcd->destComp      = destcomp;
        fcd->nComp         = numcomp;

        if (MyProc == remoteProc)
        {
            //
            // Data is local.
            //
            fcd->fillType       = FillLocally;
            fcd->localFabSource = &(*fabArray)[faindex];
        }
        else
        {
            //
            // Data is remote.
            //
            FabComTag fabComTag;

            dataAvailable               = false;
            fcd->fillType               = FillRemotely;
            fcd->localFabSource         = new FAB(intersect, numcomp);
            fcd->cacheDataAllocated     = true;
            fabComTag.fabArrayId        = fabarrayid.Id();
            fabComTag.fillBoxId         = nextFillBoxId;
            fabComTag.fabIndex          = faindex;
            fabComTag.procThatNeedsData = MyProc;
            fabComTag.procThatHasData   = remoteProc;
            fabComTag.box               = intersect;
            fabComTag.srcComp           = srccomp;
            fabComTag.destComp          = destcomp;
            fabComTag.nComp             = numcomp;
            //
            // Do not send the data yet.
            //
            fabComTagList.push_back(fabComTag);
        }

        fabCopyDescList[fabarrayid.Id()].insert(FCDMapValueType(fcd->fillBoxId,fcd));

        if (!(returnedUnfilledBoxes == 0))
        {
            unfilledBoxDomain.rmBox(intersect);
        }
    }
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (!(returnedUnfilledBoxes == 0))
    {
        unfilledBoxDomain.add(destFabBox);
    }

    std::vector< std::pair<int,Box> > isects = fabArrays[fabarrayid.Id()]->boxArray().intersections(destFabBox);

    for (int ii = 0; ii < isects.size(); ii++)
    {
        AddBoxDoIt(fabarrayid,
                   destFabBox,
                   returnedUnfilledBoxes,
                   isects[ii].first,
                   srccomp,
                   destcomp,
                   numcomp,
                   true,
                   unfilledBoxDomain);
    }

    if (!(returnedUnfilledBoxes == 0))
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        fabarrayindex,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp,
                                     bool       bUseValidBox)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (!(returnedUnfilledBoxes == 0))
    {
        unfilledBoxDomain.add(destFabBox);
    }

    AddBoxDoIt(fabarrayid,
               destFabBox,
               returnedUnfilledBoxes,
               fabarrayindex,
               srccomp,
               destcomp,
               numcomp,
               bUseValidBox,
               unfilledBoxDomain);

    if (!(returnedUnfilledBoxes == 0))
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes)
{
    return AddBox(fabarrayid,
                  destFabBox,
                  returnedUnfilledBoxes,
                  0,
                  0,
                  fabArrays[fabarrayid.Id()]->nComp(),
                  true);
}

template <class FAB>
FabArrayCopyDescriptor<FAB>::~FabArrayCopyDescriptor()
{
   clear();
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::clear ()
{
   for (unsigned int i = 0; i < fabCopyDescList.size(); ++i)
   {
       FCDMapIter fmi = fabCopyDescList[i].begin();

       for ( ; fmi != fabCopyDescList[i].end(); ++fmi)
       { 
           delete (*fmi).second;
       }
   }

   fabArrays.clear();
   fabCopyDescList.clear();
   fabComTagList.clear();

   nextFillBoxId = 0;
   dataAvailable = false;
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData (Array<int>*    snd_cache,
                                          CommDataCache* cd_cache)
{
    typedef typename FAB::value_type value_type;

    dataAvailable = true;

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1) return;

#if BL_USE_MPI

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::CollectData()");

    static bool do_alltoallv     = false;
    static bool isdef_alltoallv  = false;
    static bool do_not_use_cache = false;

    if (!isdef_alltoallv)
    {
        isdef_alltoallv = true;
        ParmParse pp("fabarray");
        pp.query("do_alltoallv",do_alltoallv);
        pp.query("do_not_use_cache",do_not_use_cache);
    }
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(CommData) == CommData::DIM*sizeof(int));

    const int MyProc = ParallelDescriptor::MyProc();

    Array<int> Snds(NProcs,0), Rcvs(NProcs,0), N_per_cpu(NProcs,0);
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (unsigned int i = 0; i < fabComTagList.size(); i++)
    {
        BL_ASSERT(fabComTagList[i].box.ok());
        BL_ASSERT(fabComTagList[i].procThatNeedsData == MyProc);
        BL_ASSERT(fabComTagList[i].procThatHasData   != MyProc);

        Rcvs[fabComTagList[i].procThatHasData]++;

        N_per_cpu[fabComTagList[i].procThatHasData] += fabComTagList[i].box.numPts()*fabComTagList[i].nComp;
    }
    BL_ASSERT(Rcvs[MyProc] == 0);

    if (snd_cache == 0 || snd_cache->size() == 0 || do_not_use_cache)
    {
        BL_MPI_REQUIRE( MPI_Alltoall(Rcvs.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<int>::type(),
                                     Snds.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<int>::type(),
                                     ParallelDescriptor::Communicator()) );
        BL_ASSERT(Snds[MyProc] == 0);

        if (snd_cache && !do_not_use_cache)
            *snd_cache = Snds;
    }
    else
    {
        Snds = *snd_cache;
    }

    Array<CommData> cd_others_need;

    if (cd_cache == 0 || !cd_cache->isValid() || do_not_use_cache)
    {
        int send_sz = 0;
        int recv_sz = 0;

        Array<int> sdispls(NProcs,0);
        Array<int> rdispls(NProcs,0);

        for (int i = 0; i < NProcs; i++) recv_sz += Snds[i];
        for (int i = 0; i < NProcs; i++) send_sz += Rcvs[i];

        for (int i = 1; i < NProcs; i++) rdispls[i] = rdispls[i-1] + Snds[i-1];
        for (int i = 1; i < NProcs; i++) sdispls[i] = sdispls[i-1] + Rcvs[i-1];

        cd_others_need.resize(recv_sz+1);     // +1 so there's always at least one element.

        Array<CommData> cd_i_need(send_sz+1); // ditto

        Array<int> index = sdispls;

        for (unsigned int i = 0; i < fabComTagList.size(); i++)
        {
            CommData data(0,
                          fabComTagList[i].fabIndex,
                          MyProc,
                          0,
                          fabComTagList[i].nComp,
                          fabComTagList[i].srcComp,
                          fabComTagList[i].fabArrayId,
                          fabComTagList[i].box);

            cd_i_need[index[fabComTagList[i].procThatHasData]++] = data;
        }
        //
        // Now increment displacements to indicate integers not CommData.
        //
        Array<int> scnts = Rcvs;
        Array<int> rcnts = Snds;

        for (int i = 0; i < NProcs; i++) scnts[i] *= CommData::DIM;
        for (int i = 0; i < NProcs; i++) rcnts[i] *= CommData::DIM;

        for (int i = 1; i < NProcs; i++) sdispls[i] *= CommData::DIM;
        for (int i = 1; i < NProcs; i++) rdispls[i] *= CommData::DIM;

        BL_MPI_REQUIRE( MPI_Alltoallv(cd_i_need.dataPtr(),
                                      scnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      cd_others_need.dataPtr(),
                                      rcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );

        if (cd_cache && !do_not_use_cache)
            *cd_cache = cd_others_need;
    }
    else
    {
        cd_others_need = cd_cache->theCommData();
    }

    if (do_alltoallv)
    {
        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), strt(NProcs,0);
	
        int send_sz = 0, recv_sz = 0;

        for (int i = 1; i < NProcs; i++)
            strt[i] = strt[i-1] + Snds[i-1];

        for (int i = 0; i < NProcs; i++)
        {
            recv_sz    += N_per_cpu[i];
            recvcnts[i] = N_per_cpu[i];

            size_t send_N = 0;
            for (int j = 0; j < Snds[i]; j++)
                send_N += cd_others_need[strt[i]+j].box().numPts() * cd_others_need[strt[i]+j].nComp();
            send_sz    += send_N;
            sendcnts[i] = send_N;

            if (i < NProcs-1)
            {
                rdispls[i+1] = rdispls[i] + recvcnts[i];
                sdispls[i+1] = sdispls[i] + sendcnts[i];
            }
        }

        FAB fab;

        BL_ASSERT((send_sz*sizeof(value_type)) < std::numeric_limits<size_t>::max());
        BL_ASSERT((recv_sz*sizeof(value_type)) < std::numeric_limits<size_t>::max());

        value_type* sendbuf = static_cast<value_type*>(BoxLib::The_Arena()->alloc(send_sz*sizeof(value_type)));
        value_type* recvbuf = static_cast<value_type*>(BoxLib::The_Arena()->alloc(recv_sz*sizeof(value_type)));

        for (int i = 0; i < NProcs; i++) 
        {
            value_type* dptr = &sendbuf[sdispls[i]];

            for (int j = 0; j < Snds[i]; j++)
            {
                const CommData& cd = cd_others_need[strt[i]+j];
                fab.resize(cd.box(), cd.nComp());
                const int N = fab.box().numPts()*fab.nComp();
                fab.copy((*fabArrays[cd.fabarrayid()])[cd.fabindex()],
                         fab.box(),
                         cd.srcComp(),
                         fab.box(),
                         0,
                         cd.nComp());
                memcpy(dptr, fab.dataPtr(), N*sizeof(value_type));
                dptr += N;
            }
        }

        cd_others_need.clear();

        BL_MPI_REQUIRE( MPI_Alltoallv(sendbuf,
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<value_type>::type(),
                                      recvbuf,
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<value_type>::type(),
                                      ParallelDescriptor::Communicator()) );

        BoxLib::The_Arena()->free(sendbuf);

        std::pair<FCDMapIter,FCDMapIter> match;

        for (int i = 0; i < NProcs; i++) 
        {
            if (Rcvs[i]) 
            {
                value_type* dptr = &recvbuf[rdispls[i]];

                for (unsigned int j = 0; j < fabComTagList.size(); j++)
                {
                    if (fabComTagList[j].procThatHasData == i)
                    {
                        const FabComTag& tag = fabComTagList[j];

                        match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                        for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
                        {
                            FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                            BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                            if (fcdp->subBox == tag.box)
                            {
                                value_type* dataPtr = fcdp->localFabSource->dataPtr();
                                BL_ASSERT(!(dataPtr == 0));
                                BL_ASSERT(fcdp->localFabSource->box() == tag.box);
                                const int N = tag.box.numPts()*tag.nComp;
                                memcpy(dataPtr, dptr, N*sizeof(value_type));
                                dptr += N;
                                break;
                            }
                        }
                    }
                }
            }
        }

        BoxLib::The_Arena()->free(recvbuf);
    }
    else
    {
        const int seqno = ParallelDescriptor::SeqNum();

        Array<value_type*> fab_data(NProcs);
        Array<MPI_Request> req_data(NProcs,MPI_REQUEST_NULL);
        Array<MPI_Request> req_cd(NProcs,MPI_REQUEST_NULL);
        Array<MPI_Status>  status(NProcs);
        //
        // Post one receive for data being requested by other CPUs.
        //
        for (int i = 0; i < NProcs; i++)
        {
            if (Rcvs[i] > 0)
            {
                BL_ASSERT(N_per_cpu[i] > 0);
                BL_ASSERT(N_per_cpu[i] < std::numeric_limits<int>::max());
                fab_data[i] = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N_per_cpu[i]*sizeof(value_type)));
                req_data[i] = ParallelDescriptor::Arecv(fab_data[i],N_per_cpu[i],i,seqno).req();
            }
        }
        //
        // Send the agglomerated FAB data.  Start Send()ing to upper neighbor.
        //
        FAB fab;

        int idx = 0;

        for (int k = 0, i = MyProc+1; k < NProcs; k++, i++)
        {
            i %= NProcs;

            if (Snds[i] > 0)
            {
                int strt = 0;
                for (int j = 0; j < i; j++)
                    strt += Snds[j];

                size_t N = 0;
                for (int j = 0; j < Snds[i]; j++)
                    N += cd_others_need[strt+j].box().numPts() * cd_others_need[strt+j].nComp();

                BL_ASSERT(N < std::numeric_limits<size_t>::max());

                value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
                value_type* dptr = data;

                for (int j = 0; j < Snds[i]; j++)
                {
                    const CommData& cd = cd_others_need[strt+j];

                    BL_ASSERT(cd.id() == 0);
                    BL_ASSERT(cd.fromproc() == i);
                    fab.resize(cd.box(), cd.nComp());
                    int count = fab.box().numPts() * fab.nComp();
                    fab.copy((*fabArrays[cd.fabarrayid()])[cd.fabindex()],
                             fab.box(),
                             cd.srcComp(),
                             fab.box(),
                             0,
                             cd.nComp());
                    memcpy(dptr, fab.dataPtr(), count*sizeof(value_type));
                    dptr += count;
                }
                BL_ASSERT(data+N == dptr);

                ParallelDescriptor::Send(data, N, i, seqno);

                BoxLib::The_Arena()->free(data);

                idx += Snds[i];
            }
        }
        //
        // Now receive and unpack FAB data.
        //
        std::pair<FCDMapIter,FCDMapIter> match;

        int NWaits = 0;
        for (int i = 0; i < NProcs; i++)
            if (req_data[i] != MPI_REQUEST_NULL)
                NWaits++;

        Array<int> indx(NProcs);

        for (int completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(req_data, completed, indx, status);

            for (int k = 0; k < completed; k++)
            {
                int          Processed = 0;
                value_type*  dptr      = fab_data[indx[k]];

                BL_ASSERT(!(dptr == 0));

                for (unsigned int j = 0; j < fabComTagList.size(); j++)
                {
                    if (fabComTagList[j].procThatHasData == indx[k])
                    {
                        const FabComTag& tag = fabComTagList[j];

                        match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                        for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
                        {
                            FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                            BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                            if (fcdp->subBox == tag.box)
                            {
                                value_type* dataPtr = fcdp->localFabSource->dataPtr();
                                BL_ASSERT(!(dataPtr == 0));
                                BL_ASSERT(fcdp->localFabSource->box() == tag.box);
                                int N = tag.box.numPts()*tag.nComp;
                                memcpy(dataPtr, dptr, N*sizeof(value_type));
                                dptr += N;
                                Processed++;
                                break;
                            }
                        }
                    }
                }

                BL_ASSERT(Processed == Rcvs[indx[k]]);

                BoxLib::The_Arena()->free(fab_data[indx[k]]);
            }
        }
    }

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                      const FillBoxId& fillboxid,
                                      FAB&             destFab)
{
    BL_ASSERT(dataAvailable);

    std::pair<FCDMapIter,FCDMapIter> match = fabCopyDescList[faid.Id()].equal_range(fillboxid.Id());

    for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
    {
        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

        BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

        destFab.copy(*fcdp->localFabSource,
                     fcdp->subBox,
                     fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                     fcdp->subBox,
                     fcdp->destComp,
                     fcdp->nComp);
    }
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                        const FillBoxId& fillboxid,
                                        FAB&             destFab,
                                        const Box&       destBox)
{
    BL_ASSERT(dataAvailable);

    FCDMapIter fmi = fabCopyDescList[faid.Id()].lower_bound(fillboxid.Id());

    BL_ASSERT(!(fmi == fabCopyDescList[faid.Id()].end()));

    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

    BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

    BL_ASSERT(fcdp->subBox.sameSize(destBox));

    destFab.copy(*fcdp->localFabSource,
                 fcdp->subBox,
                 fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                 destBox,
                 fcdp->destComp,
                 fcdp->nComp);

    BL_ASSERT(++fmi == fabCopyDescList[faid.Id()].upper_bound(fillboxid.Id()));
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::PrintStats () const
{
    const int MyProc = ParallelDescriptor::MyProc();

    std::cout << "----- "
         << MyProc
         << ":  Parallel stats for FabArrayCopyDescriptor:" << '\n';

    for (int fa = 0; fa < fabArrays.size(); ++fa)
    {
      std::cout << "fabArrays["
             << fa
             << "]->boxArray() = "
             << fabArrays[fa]->boxArray()
             << '\n';
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary ()
{
    FillBoundary(0, nComp());
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp,
                             int ncomp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::FillBoundary()");

    if ( n_grow <= 0 ) return;

    FabArrayCopyDescriptor<FAB> mfcd;
    FabArrayBase::SI&           si   = FabArrayBase::TheFBsirec(scomp,ncomp,*this);
    const FabArrayId            mfid = mfcd.RegisterFabArray(this);
    //
    // Add boxes we need to collect.
    //
    for (unsigned int i = 0; i < si.m_sirec.size(); i++)
    {
        si.m_sirec[i].m_fbid = mfcd.AddBox(mfid,
                                           si.m_sirec[i].m_bx,
                                           0,
                                           si.m_sirec[i].m_j,
                                           scomp,
                                           scomp,
                                           ncomp);
    }

    mfcd.CollectData(&si.m_cache,&si.m_commdata);

    for (unsigned int i = 0; i < si.m_sirec.size(); i++)
    {
        BL_ASSERT(DistributionMap()[si.m_sirec[i].m_i] == ParallelDescriptor::MyProc());
        //
        // Directly fill the FAB.
        //
        mfcd.FillFab(mfid,si.m_sirec[i].m_fbid,(*this)[si.m_sirec[i].m_i]);
    }
}

#endif /*BL_FABARRAY_H*/
