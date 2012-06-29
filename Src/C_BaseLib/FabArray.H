
#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

#include <iostream>
#include <cstring>
#include <limits>
#include <map>
#include <deque>
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
// This is meant to be a concrete class not a polymorphic one.
//

class FabArrayBase
{
public:

    FabArrayBase ();

    ~FabArrayBase();
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
    const Box& box (int K) const;
    //
    // Returns the Kth FABs Box in the FabArray.
    // That is, the region the Kth fab is actually defined on.
    //
    const Box fabbox (int K) const;
    //
    // Returns the number of FABs in the FabArray..
    //
    int size () const;
    //
    // Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const;
    //
    // Returns constant reference to indices in the FabArray that we own.
    //
    const Array<int>& IndexMap () const;
    //
    // Flush the cache of self-intersection info used by FillBoundary.
    //
    static void FlushSICache ();
    //
    // The size of the cache of self-intersection info.
    //
    static int SICacheSize ();
    //
    // Used to cache some CommData stuff in CollectData().
    //
    struct CommDataCache
    {
        CommDataCache ();

        void operator= (const Array<ParallelDescriptor::CommData>& rhs);

        bool isValid () const { return m_valid; }

        Array<ParallelDescriptor::CommData>& theCommData () { return m_commdata; }

        const Array<ParallelDescriptor::CommData>& theCommData () const { return m_commdata; }

        Array<ParallelDescriptor::CommData> m_commdata;
        bool                                m_valid;
    };
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
        //
        // A handy typedef.
        //
        typedef std::vector<SIRec> SIRecVector;

        Array<int>          m_cache;    // Snds cached for CollectData().
        CommDataCache       m_commdata; // Yet another cache for CollectData().
        SIRecVector         m_sirec;
        BoxArray            m_ba;
        DistributionMapping m_dm;
        int                 m_ngrow;
        bool                m_reused;
    };
    //
    // Used in copy()
    //
    struct CopyComTag
    {
        CopyComTag () {}
        CopyComTag (const CopyComTag& cct)
            :
            box(cct.box),
            fabIndex(cct.fabIndex),
            srcIndex(cct.srcIndex)
            {}

        Box box;
        int fabIndex;
        int srcIndex;
    };
    //
    // When copy()ing from one FabArray to another we can do a copy or an add.
    //
    enum CpOp { COPY = 0, ADD = 1 };

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
        int                 m_rcvs_numpts;
        //
        // Some useful typedefs.
        //
        typedef std::vector<CopyComTag> VectorOfCopyComTags;

        typedef std::map< int,std::vector<CopyComTag> > MapOfCopyComTagVectors;
        //
        // The cache of local and send/recv info per FabArray::copy().
        //
        VectorOfCopyComTags    m_LocTags;
        MapOfCopyComTagVectors m_SndTags;
        MapOfCopyComTagVectors m_RcvTags;
        //
        // Has this CPC ever been reused in the cache?
        //
        bool m_reused;
    };
    //
    // Used for collecting information used in communicating FABs.
    //
    struct FabComTag
    {
        int fromProc;
        int toProc;
        int fabIndex;
        int fineIndex;
        int srcComp;
        int destComp;
        int nComp;
        int face;
        int fabArrayId;
        int fillBoxId;
        int procThatNeedsData;
        int procThatHasData;
        Box box;

        FabComTag ();
    };
    //
    // Builds self-intersection records.
    //
    static SI& BuildFBsirec (const SI& si, const FabArrayBase& mf);
    //
    // Returns cached self-intersection records or builds them.
    //
    static SI& TheFBsirec (int scomp, int ncomp, const FabArrayBase& mf);
    //
    // Use MPI_Alltoallv() instead of MPI_IRecv()/MPI_Send().
    //
    // Turn on via ParmParse using "fabarray.do_alltoallv=1" in inputs file.
    //
    // Default is false.
    //
    static bool do_alltoallv;
    //
    // Use MPI_Asend() instead of MPI_Send() in CollectData() and copy().
    //
    // Turn on via ParmParse using "fabarray.do_async_sends=1" in inputs file.
    //
    // Default is false.
    //
    static bool do_async_sends;
    //
    // Do not use the CollectData() cache.
    //
    // Turn on via ParmParse using "fabarray.do_not_use_cache=1" in inputs file.
    //
    // Default is false.
    //
    static bool do_not_use_cache;
    //
    // Print out some stuff; default is false.
    //
    static bool verbose;
    //
    // Initialize from ParmParse with "fabarray" prefix.
    //
    static void Initialize ();
    static void Finalize ();

protected:
    //
    // The data ...
    //
    mutable BoxArray    boxarray;   /* So FabSets can modify'm */
    DistributionMapping distributionMap;
    Array<int>          indexMap;
    int                 n_grow;
    int                 n_comp;
};

class MFIter
{
public:
    //
    // Construct a MFIter.
    //
    explicit MFIter (const FabArrayBase& fabarray);
    //
    // Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    // Returns the Box of the FAB at which we currently point.
    //
    const Box fabbox () const;
    //
    // Increments iterator to the next FAB in the FabArray that we own.
    //
    void operator++ ();
    //
    // Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid ();
    //
    // The index into the underlying BoxArray of the current FAB.
    //
    int index () const;
    //
    // Constant reference to FabArray over which we're iterating.
    //
    const FabArrayBase& theFabArrayBase () const;

private:

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
// Possible values are: Fab_noallocate and Fab_allocate.
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
    // Constructs an empty FabArray<FAB>.
    //
    FabArray ();
    //
    // Construct a FabArray<FAB> with a valid region defined by bxs
    // and a region of definition defined by the grow factor ngrow
    // and the number of components nvar.
    // If mem_mode is defined to be Fab_allocate then FABs are
    // allocated for each Box in the BoxArray.  The size of the Kth
    // FAB is given by bxs[K] grown by ngrow.  If mem_mode is defined
    // to be Fab_noallocate, then no FABs are allocated at this time,
    // but can be defined later.  The number of components in each
    // FAB is not specified and is expected to be implicit in the
    // definition of the FAB class.  That is, the FAB constructor will
    // take only a Box argument.  Call this constructor number two.
    //
    FabArray (const BoxArray& bxs,
              int             nvar,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);
    //
    // The destructor -- deletes all FABs in the array.
    //
    ~FabArray ();
    //
    // Define this FabArray identically to that performed by
    // the constructor having an analogous function signature.
    // This is only valid if this FabArray was defined using
    // the default constructor.
    //
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
    // Returns true if the FabArray is well-defined.  That is,
    // if FABs are allocated for each Box in the BoxArray and the
    // sizes of the FABs and the number of components are consistent
    // with the definition of the FabArray.
    //
    bool ok () const;
    //
    // Returns a constant reference to the FAB associated with the Kth element.
    //
    const FAB& operator[] (const MFIter& mfi) const;

    const FAB& get (const MFIter& mfi) const;
    //
    // Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (const MFIter& mfi);

    FAB& get (const MFIter& mfi);
    //
    // Returns a constant reference to the FAB associated with the Kth element.
    //
    const FAB& operator[] (int K) const;

    const FAB& get (int K) const;
    //
    // Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (int K);

    FAB& get (int K);
    //
    // Explicitly set the Kth FAB in the FabArray to point to elem.
    //
    void setFab (int K, FAB* elem);
    //
    // Releases FAB memory in the FabArray.
    //
    void clear ();
    //
    // Set all components in the valid region of each FAB to val.
    //
    void setVal (value_type val);
    void operator= (const value_type& val);
    //
    // Set the value of num_comp components in the valid region of
    // each FAB in the FabArray, starting at component comp to val.
    // Also set the value of nghost boundary cells.
    //
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Set the value of num_comp components in the valid region of
    // each FAB in the FabArray, starting at component comp, as well
    // as nghost boundary cells, to val, provided they also intersect
    // with the Box region.
    //
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Set all components in the valid region of each FAB in the
    // FabArray to val, including nghost boundary cells.
    //
    void setVal (value_type val,
                 int        nghost);
    //
    // Set all components in the valid region of each FAB in the
    // FabArray to val, including nghost boundary cells, that also
    // intersect the Box region.
    //
    void setVal (value_type val,
                 const Box& region,
                 int        nghost);
    //
    // Set all values in the boundary region to val.
    //
    void setBndry (value_type val);
    //
    // Set ncomp values in the boundary region, starting at start_comp to val.
    //
    void setBndry (value_type val,
                   int        strt_comp,
                   int        ncomp);
    //
    // This function copies data from fa to this FabArray.  Each FAB
    // in fa is intersected with all FABs in this FabArray and a copy
    // is performed on the region of intersection.  The intersection
    // is restricted to the valid region of each FAB.
    //
    void copy (const FabArray<FAB>& fa, CpOp op = FabArrayBase::COPY);
    //
    // This function copies data from src to this FabArray.  Each FAB
    // in src is intersected with all FABs in this FabArray and a copy
    // is performed on the region of intersection.  The intersection
    // is restricted to the num_comp components starting at src_comp
    // in the FabArray src, with the destination components in this
    // FabArray starting at dest_comp.
    //
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
               CpOp                 op = FabArrayBase::COPY);
    //
    // Copies the values contained in the intersection of the
    // valid region of this FabArray with the FAB dest into dest.
    //
    void copy (FAB& dest) const;
    //
    // Copies the values contained in the intersection of the
    // valid region of this FabArray with the FAB dest and the Box
    // subbox into that subregion of dest.
    //
    void copy (FAB&       dest,
               const Box& subbox) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid region of this FabArray, starting at
    // component src_comp, with the FAB dest into dest, starting at
    // component dest_comp in dest.
    //
    void copy (FAB& dest,
               int  src_comp,
               int  dest_comp,
               int  num_comp) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid region of this FabArray, starting at
    // component src_comp, with the FAB dest and the Box subbox, into
    // dest, starting at component dest_comp in dest.
    //
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

    void defineDoit (const BoxArray&            bxs,
                     int                        nvar,
                     int                        ngrow,
                     FabAlloc                   mem_mode,
                     const DistributionMapping* dm);

    std::map<int,FAB*> m_fabs;

    typedef typename std::map<int,FAB*>::value_type ValueType;

    typedef typename std::map<int,FAB*>::iterator Iterator;

    typedef typename std::map<int,FAB*>::const_iterator Const_Iterator;

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

  typedef typename FCDMap::value_type     FCDMapValueType;
  typedef typename FCDMap::iterator       FCDMapIter;
  typedef typename FCDMap::const_iterator FCDMapConstIter;

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

    void CollectData (Array<int>*                  snd_cache = 0,
                      FabArrayBase::CommDataCache* cd_cache  = 0);

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

private:
    //
    // These are disallowed.
    //
    FabArrayCopyDescriptor (const FabArrayCopyDescriptor<FAB>&);
    FabArrayCopyDescriptor<FAB>& operator= (const FabArrayCopyDescriptor<FAB> &);
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
    //
    // The data.
    //
    std::vector<FabArray<FAB>*>         fabArrays;
    std::deque<FCDMap>                  fabCopyDescList;
    std::deque<FabArrayBase::FabComTag> fabComTagList;
    int                                 nextFillBoxId;
    bool                                dataAvailable;
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
const Array<int>&
FabArrayBase::IndexMap () const
{
    return indexMap;
}

inline
MFIter::MFIter (const FabArrayBase& fabarray)
    :
    fabArray(fabarray),
    currentIndex(0)
{}

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
    return fabArray.IndexMap()[currentIndex];
}

inline
void
MFIter::operator++ ()
{
    ++currentIndex;
}

inline
bool
MFIter::isValid ()
{
    return currentIndex < fabArray.IndexMap().size();
}

inline
const Box&
MFIter::validbox () const
{
    return fabArray.box(fabArray.IndexMap()[currentIndex]);
}

inline
const Box
MFIter::fabbox () const
{
    return fabArray.fabbox(fabArray.IndexMap()[currentIndex]);
}

template <class FAB>
bool
FabArray<FAB>::defined (int i) const
{
    bool result = false;

    Const_Iterator it = m_fabs.find(i);

    if (it != m_fabs.end()) result = (it->second != 0);

    return result;
}

template <class FAB>
const FAB&
FabArray<FAB>::operator[] (const MFIter& mfi) const
{
    Const_Iterator it = m_fabs.find(mfi.index());

    BL_ASSERT(it != m_fabs.end());

    return *(it->second);
}

template <class FAB>
inline
const FAB&
FabArray<FAB>::get (const MFIter& mfi) const
{
    return this->operator[](mfi);
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    Iterator it = m_fabs.find(mfi.index());

    BL_ASSERT(it != m_fabs.end());

    return *(it->second);
}

template <class FAB>
inline
FAB&
FabArray<FAB>::get (const MFIter& mfi)
{
    return this->operator[](mfi);
}

template <class FAB>
const FAB&
FabArray<FAB>::operator[] (int K) const
{
    Const_Iterator it = m_fabs.find(K);

    BL_ASSERT(it != m_fabs.end());

    return *(it->second);
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (int K)
{
    Iterator it = m_fabs.find(K);

    BL_ASSERT(it != m_fabs.end());

    return *(it->second);
}

template <class FAB>
inline
const FAB&
FabArray<FAB>::get (int K) const
{
    return this->operator[](K);
}

template <class FAB>
inline
FAB&
FabArray<FAB>::get (int K)
{
    return this->operator[](K);
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    for (Iterator it = m_fabs.begin(), End = m_fabs.end(); it != End; ++it)
        delete it->second;

    m_fabs.clear();
    boxarray.clear();
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
{}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc)
{
    define(bxs,nvar,ngrow,alloc);
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
    clear();
}

template <class FAB>
bool
FabArray<FAB>::ok () const
{
    long isok = true;

    for (MFIter fai(*this); fai.isValid() && isok; ++fai)
    {
        if (defined(fai.index()))
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
FabArray<FAB>::defineDoit (const BoxArray&            bxs,
                           int                        nvar,
                           int                        ngrow,
                           FabAlloc                   alloc,
                           const DistributionMapping* dm)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);

    n_grow = ngrow;
    n_comp = nvar;

    boxarray.define(bxs);

    if (dm == 0)
    {
        distributionMap.define(boxarray,ParallelDescriptor::NProcs());
    }
    else
    {
        BL_ASSERT(dm->ProcessorMap().size() == bxs.size()+1);

        distributionMap = *dm;
    }

    int cnt = 0;
    const int MyProc = ParallelDescriptor::MyProc();
    for (int i = 0, N = boxarray.size(); i < N; i++)
        if (distributionMap[i] == MyProc)
            cnt++;

    indexMap.resize(cnt);

    cnt = 0;
    for (int i = 0, N = boxarray.size(); i < N; i++)
        if (distributionMap[i] == MyProc)
            indexMap[cnt++] = i;

    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray& bxs,
                       int             nvar,
                       int             ngrow,
                       FabAlloc        alloc)
{
    defineDoit(bxs,nvar,ngrow,alloc,0);
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       int                        nvar,
                       int                        ngrow,
                       const DistributionMapping& dm,
                       FabAlloc                   alloc)
{
    defineDoit(bxs,nvar,ngrow,alloc,&dm);
}

template <class FAB>
void
FabArray<FAB>::AllocFabs ()
{
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Box tmp = BoxLib::grow(fai.validbox(), n_grow);

        m_fabs.insert(m_fabs.end(),ValueType(fai.index(),new FAB(tmp, n_comp)));
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
    BL_ASSERT(!this->defined(boxno));
    BL_ASSERT(distributionMap[boxno] == ParallelDescriptor::MyProc());

    m_fabs[boxno] = elem;
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

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp,
                     CpOp                 op)
{
    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);

    if (boxarray == src.boxarray && distributionMap == src.distributionMap)
    {
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            if (op == FabArrayBase::COPY)
            {
                get(fai).copy(src[fai],fai.validbox(),scomp,fai.validbox(),dcomp,ncomp);
            }
            else
            {
                get(fai).plus(src[fai],fai.validbox(),fai.validbox(),scomp,dcomp,ncomp);
            }
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
        for (CPC::VectorOfCopyComTags::iterator it = thecpc.m_LocTags.begin(),
                 End = thecpc.m_LocTags.end();
             it != End;
             ++it)
        {
            CopyComTag& tag = *it;

            BL_ASSERT(distributionMap[tag.fabIndex] == MyProc);
            BL_ASSERT(src.distributionMap[tag.srcIndex] == MyProc);

            if (op == FabArrayBase::COPY)
            {
                get(tag.fabIndex).copy(src[tag.srcIndex],tag.box,scomp,tag.box,dcomp,ncomp);
            }
            else
            {
                get(tag.fabIndex).plus(src[tag.srcIndex],tag.box,tag.box,scomp,dcomp,ncomp);
            }
        }
    }
    else
    {
        CopyComTag tag;

        std::vector< std::pair<int,Box> > isects;

        for (int i = 0, N = size(); i < N; i++)
        {
            src.boxarray.intersections(boxarray[i],isects);

            for (int ii = 0, M = isects.size(); ii < M; ii++)
            {
                const Box& bx  = isects[ii].second;
                const int  iii = isects[ii].first;

                if (distributionMap[i] == MyProc)
                {
                    if (src.distributionMap[iii] == MyProc)
                    {
                        if (op == FabArrayBase::COPY)
                        {
                            get(i).copy(src[iii],bx,scomp,bx,dcomp,ncomp);
                        }
                        else
                        {
                            get(i).plus(src[iii],bx,bx,scomp,dcomp,ncomp);
                        }

                        tag.box      = bx;
                        tag.fabIndex = i;
                        tag.srcIndex = iii;

                        thecpc.m_LocTags.push_back(tag);
                    }
                    else
                    {
                        tag.box      = bx;
                        tag.fabIndex = i;

                        thecpc.m_RcvTags[src.distributionMap[iii]].push_back(tag);

                        thecpc.m_rcvs_numpts += bx.numPts();
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
        //
        // Squish out any spare capacity.
        //
        if (!thecpc.m_LocTags.empty() && thecpc.m_LocTags.capacity() > thecpc.m_LocTags.size())
        {
            CPC::VectorOfCopyComTags tmp(thecpc.m_LocTags);
            tmp.swap(thecpc.m_LocTags);
        }

        for (CPC::MapOfCopyComTagVectors::iterator m_it = thecpc.m_SndTags.begin(),
                 m_End = thecpc.m_SndTags.end();
             m_it != m_End;
             ++m_it)
        {
            if (!m_it->second.empty() && m_it->second.capacity() > m_it->second.size())
            {
                CPC::VectorOfCopyComTags tmp(m_it->second);
                tmp.swap(m_it->second);
            }
        }

        for (CPC::MapOfCopyComTagVectors::iterator m_it = thecpc.m_RcvTags.begin(),
                 m_End = thecpc.m_RcvTags.end();
             m_it != m_End;
             ++m_it)
        {
            if (!m_it->second.empty() && m_it->second.capacity() > m_it->second.size())
            {
                CPC::VectorOfCopyComTags tmp(m_it->second);
                tmp.swap(m_it->second);
            }
        }
    }
    
    if (NProcs == 1) return;

#ifdef BL_USE_MPI

    const Real strt = ParallelDescriptor::second();

    if (FabArrayBase::do_alltoallv)
    {
        size_t Nsend = 0, Nrecv = 0;

        Array<int> sendcount(NProcs,0), recvcount(NProcs,0);
        Array<int> senddspls(NProcs,0), recvdspls(NProcs,0);
        
        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_SndTags.begin(),
                 m_End = thecpc.m_SndTags.end();
             m_it != m_End;
             ++m_it)
        {
            int Ns = 0;
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                Ns += it->box.numPts()*ncomp;
            }
            Nsend                 += Ns;
            sendcount[m_it->first] = Ns;
        }

        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_RcvTags.begin(),
                 m_End = thecpc.m_RcvTags.end();
             m_it != m_End;
             ++m_it)
        {
            int Nr = 0;
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                Nr += it->box.numPts()*ncomp;
            }
            Nrecv                 += Nr;
            recvcount[m_it->first] = Nr;
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

        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_SndTags.begin(),
                 m_End = thecpc.m_SndTags.end();
             m_it != m_End;
             ++m_it)
        {
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->box;
                fab.resize(bx, ncomp);
                fab.copy(src[it->fabIndex],bx,scomp,bx,0,ncomp);
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

        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_RcvTags.begin(),
                 m_End = thecpc.m_RcvTags.end();
             m_it != m_End;
             ++m_it)
        {
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->box;
                fab.resize(bx, ncomp);
                const int N = bx.numPts()*ncomp;
                memcpy(fab.dataPtr(), recvdptr, N*sizeof(value_type));

                if (op == FabArrayBase::COPY)
                {
                    get(it->fabIndex).copy(fab,bx,0,bx,dcomp,ncomp);
                }
                else
                {
                    get(it->fabIndex).plus(fab,bx,bx,0,dcomp,ncomp);
                }
                recvdptr += N;
            }
        }

        BL_ASSERT(recvdata+Nrecv == recvdptr);

        BoxLib::The_Arena()->free(recvdata);
    }
    else
    {
        const int N_rcvs = thecpc.m_RcvTags.size();
        const int N_snds = thecpc.m_SndTags.size();

        Array<MPI_Status>  status(N_rcvs);
        Array<int>         recv_from(N_rcvs), index(N_rcvs);
        Array<value_type*> recv_data(N_rcvs), send_data;
        Array<MPI_Request> recv_reqs(N_rcvs), send_reqs;

        if (FabArrayBase::do_async_sends && N_snds > 0)
        {
            send_data.resize(N_snds);
            send_reqs.resize(N_snds);
        }
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        BL_ASSERT((thecpc.m_rcvs_numpts*ncomp*sizeof(value_type)) < std::numeric_limits<int>::max());

        value_type* the_recv_data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(thecpc.m_rcvs_numpts*ncomp*sizeof(value_type)));

        const int seqno = ParallelDescriptor::SeqNum();

        int idx = 0, offset = 0;

        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_RcvTags.begin(),
                 m_End = thecpc.m_RcvTags.end();
             m_it != m_End;
             ++m_it, ++idx)
        {
            int N = 0;
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                N += it->box.numPts()*ncomp;
            }
            BL_ASSERT(N < std::numeric_limits<int>::max());

            recv_data[idx] = &the_recv_data[offset];
            recv_from[idx] = m_it->first;
            recv_reqs[idx] = ParallelDescriptor::Arecv(recv_data[idx],N,m_it->first,seqno).req();

            offset += N;
        }

        BL_ASSERT(offset == thecpc.m_rcvs_numpts*ncomp);

        FAB fab;
        //
        // Send the FAB data.
        //
        idx = 0;

        for (CPC::MapOfCopyComTagVectors::const_iterator m_it = thecpc.m_SndTags.begin(),
                 m_End = thecpc.m_SndTags.end();
             m_it != m_End;
             ++m_it, ++idx)
        {
            int N = 0;
            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                N += it->box.numPts()*ncomp;
            }
            BL_ASSERT(N < std::numeric_limits<int>::max());

            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
            value_type* dptr = data;

            for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->box;
                fab.resize(bx, ncomp);
                fab.copy(src[it->fabIndex],bx,scomp,bx,0,ncomp);
                const int count = bx.numPts()*ncomp;
                memcpy(dptr, fab.dataPtr(), count*sizeof(value_type));
                dptr += count;
            }
            BL_ASSERT(data+N == dptr);

            if (FabArrayBase::do_async_sends)
            {
                send_data[idx] = data;
                send_reqs[idx] = ParallelDescriptor::Asend(data,N,m_it->first,seqno).req();
            }
            else
            {
                ParallelDescriptor::Send(data,N,m_it->first,seqno);
                BoxLib::The_Arena()->free(data);
            }
        }
        //
        // Now receive and unpack FAB data as it becomes available.
        //
        CPC::MapOfCopyComTagVectors::const_iterator m_it;

        for (int NWaits = N_rcvs, completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(recv_reqs, completed, index, status);

            for (int k = 0; k < completed; k++)
            {
                value_type* dptr = recv_data[index[k]];

                BL_ASSERT(dptr != 0);

                m_it = thecpc.m_RcvTags.find(recv_from[index[k]]);

                BL_ASSERT(m_it != thecpc.m_RcvTags.end());

                for (CPC::VectorOfCopyComTags::const_iterator it = m_it->second.begin(),
                         End = m_it->second.end();
                     it != End;
                     ++it)
                {
                    const Box& bx = it->box;
                    fab.resize(bx, ncomp);
                    const int N = bx.numPts() * ncomp;
                    memcpy(fab.dataPtr(), dptr, N*sizeof(value_type));

                    if (op == FabArrayBase::COPY)
                    {
                        get(it->fabIndex).copy(fab,bx,0,bx,dcomp,ncomp);
                    }
                    else
                    {
                        get(it->fabIndex).plus(fab,bx,bx,0,dcomp,ncomp);
                    }
                    dptr += N;
                }
            }
        }

        BoxLib::The_Arena()->free(the_recv_data);

        if (FabArrayBase::do_async_sends && N_snds > 0)
        {
            //
            // Now grok the asynchronous send buffers.
            //
            status.resize(N_snds);

            BL_MPI_REQUIRE( MPI_Waitall(N_snds,
                                        send_reqs.dataPtr(),
                                        status.dataPtr()) );
            //
            // And free up send buffer space.
            //
            for (int i = 0; i < N_snds; i++)
                BoxLib::The_Arena()->free(send_data[i]);
        }
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      tot    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(tot,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "FabArray::copy() time = " << tot << std::endl;
    }

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src, CpOp op)
{
    copy(src,0,0,nComp(),op);
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

    if (ParallelDescriptor::NProcs() == 1)
    {
        for (int j = 0, N = size(); j < N; ++j)
        {
            if (boxarray[j].intersects(subbox))
            {
                Box destbox = boxarray[j] & subbox;

                dest.copy(get(j),destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    FAB ovlp;

    for (int i = 0, N = size(); i < N; i++)
    {
        if (subbox.intersects(boxarray[i]))
        {
            Box bx = subbox & boxarray[i];

            ovlp.resize(bx,ncomp);

            if (ParallelDescriptor::MyProc() == distributionMap[i])
            {
                ovlp.copy(get(i),bx,scomp,bx,0,ncomp);
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

    FabArrayId result(fabArrays.size());

    fabArrays.push_back(fabarray);  /* Bump size() by one */

    fabCopyDescList.push_back(FCDMap());

    return result;
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
            FabArrayBase::FabComTag fabComTag;

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

        if (returnedUnfilledBoxes != 0)
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

    if (returnedUnfilledBoxes != 0)
    {
        unfilledBoxDomain.add(destFabBox);
    }

    std::vector< std::pair<int,Box> > isects;

    fabArrays[fabarrayid.Id()]->boxArray().intersections(destFabBox,isects);

    for (int ii = 0, N = isects.size(); ii < N; ii++)
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

    if (returnedUnfilledBoxes != 0)
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

    if (returnedUnfilledBoxes != 0)
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

    if (returnedUnfilledBoxes != 0)
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
    for (unsigned int i = 0, N = fabCopyDescList.size(); i < N; ++i)
    {
        for (FCDMapIter fmi = fabCopyDescList[i].begin(), End = fabCopyDescList[i].end();
             fmi != End;
             ++fmi)
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
FabArrayCopyDescriptor<FAB>::CollectData (Array<int>*                  snd_cache,
                                          FabArrayBase::CommDataCache* cd_cache)
{
    typedef typename FAB::value_type value_type;

    dataAvailable = true;

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1) return;

#if BL_USE_MPI
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(ParallelDescriptor::CommData) == ParallelDescriptor::CommData::DIM*sizeof(int));

    const int MyProc = ParallelDescriptor::MyProc();

    int Total_Rcvs = 0;

    Array<int> Snds(NProcs,0), Rcvs(NProcs,0), N_per_cpu(NProcs,0);
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (std::deque<FabArrayBase::FabComTag>::const_iterator it = fabComTagList.begin(),
             End = fabComTagList.end();
         it != End;
         ++it)
    {
        BL_ASSERT(it->box.ok());
        BL_ASSERT(it->procThatNeedsData == MyProc);
        BL_ASSERT(it->procThatHasData   != MyProc);

        const int count = (it->box.numPts())*(it->nComp);

        Total_Rcvs += count;

        Rcvs[it->procThatHasData]++;

        N_per_cpu[it->procThatHasData] += count;
    }
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_ASSERT((Total_Rcvs*sizeof(value_type)) < std::numeric_limits<int>::max());

    if (snd_cache == 0 || snd_cache->size() == 0 || FabArrayBase::do_not_use_cache)
    {
        BL_MPI_REQUIRE( MPI_Alltoall(Rcvs.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<int>::type(),
                                     Snds.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<int>::type(),
                                     ParallelDescriptor::Communicator()) );
        BL_ASSERT(Snds[MyProc] == 0);

        if (snd_cache && !FabArrayBase::do_not_use_cache)
        {
            *snd_cache = Snds;
        }
    }
    else
    {
        Snds = *snd_cache;
    }

    Array<ParallelDescriptor::CommData> cd_others_need, cd_that_i_need;

    if (cd_cache == 0 || !cd_cache->isValid() || FabArrayBase::do_not_use_cache)
    {
        int send_sz = 0, recv_sz = 0;

        Array<int> sdispls(NProcs,0), rdispls(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recv_sz += Snds[i];
            send_sz += Rcvs[i];

            if (i > 0)
            {
                rdispls[i] = rdispls[i-1] + Snds[i-1];
                sdispls[i] = sdispls[i-1] + Rcvs[i-1];
            }
        }

        Array<int> index(sdispls);

        cd_others_need.resize(recv_sz+1); // +1 so there's always at least one element.
        cd_that_i_need.resize(send_sz+1); // +1 so there's always at least one element.

        for (std::deque<FabArrayBase::FabComTag>::const_iterator it = fabComTagList.begin(),
                 End = fabComTagList.end();
             it != End;
             ++it)
        {
            ParallelDescriptor::CommData data(0,
                                              it->fabIndex,
                                              MyProc,
                                              0,
                                              it->nComp,
                                              it->srcComp,
                                              it->fabArrayId,
                                              it->box);

            cd_that_i_need[index[it->procThatHasData]++] = data;
        }
        //
        // Now increment displacements to indicate integers not CommData.
        //
        Array<int> scnts(Rcvs), rcnts(Snds);

        for (int i = 0; i < NProcs; i++)
        {
            scnts[i] *= ParallelDescriptor::CommData::DIM;
            rcnts[i] *= ParallelDescriptor::CommData::DIM;

            if (i > 0)
            {
                sdispls[i] *= ParallelDescriptor::CommData::DIM;
                rdispls[i] *= ParallelDescriptor::CommData::DIM;
            }
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(cd_that_i_need.dataPtr(),
                                      scnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      cd_others_need.dataPtr(),
                                      rcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );

        if (cd_cache && !FabArrayBase::do_not_use_cache)
        {
            *cd_cache = cd_others_need;
        }
    }
    else
    {
        cd_others_need = cd_cache->theCommData();
    }

    if (FabArrayBase::do_alltoallv)
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

            int send_N = 0;
            for (int j = 0; j < Snds[i]; j++)
            {
                const ParallelDescriptor::CommData& cd = cd_others_need[strt[i]+j];

                send_N += cd.box().numPts() * cd.nComp();
            }
            send_sz    += send_N;
            sendcnts[i] = send_N;

            if (i < NProcs-1)
            {
                rdispls[i+1] = rdispls[i] + recvcnts[i];
                sdispls[i+1] = sdispls[i] + sendcnts[i];
            }
        }

        FAB fab;

        BL_ASSERT((send_sz*sizeof(value_type)) < std::numeric_limits<int>::max());
        BL_ASSERT((recv_sz*sizeof(value_type)) < std::numeric_limits<int>::max());

        value_type* sendbuf = static_cast<value_type*>(BoxLib::The_Arena()->alloc(send_sz*sizeof(value_type)));
        value_type* recvbuf = static_cast<value_type*>(BoxLib::The_Arena()->alloc(recv_sz*sizeof(value_type)));

        for (int i = 0; i < NProcs; i++) 
        {
            value_type* dptr = &sendbuf[sdispls[i]];

            for (int j = 0; j < Snds[i]; j++)
            {
                const ParallelDescriptor::CommData& cd = cd_others_need[strt[i]+j];
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
            if (Rcvs[i] <= 0) continue;

            value_type* dptr = &recvbuf[rdispls[i]];

            for (std::deque<FabArrayBase::FabComTag>::const_iterator it = fabComTagList.begin(),
                     End = fabComTagList.end();
                 it != End;
                 ++it)
            {
                if (it->procThatHasData == i)
                {
                    const FabArrayBase::FabComTag& tag = *it;

                    match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                    for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
                    {
                        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                        BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                        if (fcdp->subBox == tag.box)
                        {
                            value_type* dataPtr = fcdp->localFabSource->dataPtr();
                            BL_ASSERT(dataPtr != 0);
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

        BoxLib::The_Arena()->free(recvbuf);
    }
    else
    {
        int N_rcvs = 0, N_snds = 0;

        for (int i = 0; i < NProcs; i++)
            if (Rcvs[i] > 0)
                N_rcvs++;
        for (int i = 0; i < NProcs; i++)
            if (Snds[i] > 0)
                N_snds++;

        Array<int> index(N_rcvs);
        Array<int> roffset(N_rcvs);

        Array<int> who_R(N_rcvs), n_R(N_rcvs), count_R(N_rcvs);
        Array<int> who_S(N_snds), n_S(N_snds), count_S(N_snds);

        Array<MPI_Status>  status(N_rcvs);
        Array<value_type*> send_data;
        Array<MPI_Request> recv_reqs(N_rcvs), send_reqs;

        if (FabArrayBase::do_async_sends && N_snds > 0)
        {
            send_data.resize(N_snds);
            send_reqs.resize(N_snds);
        }

        for (int i = 0, j = 0; i < NProcs; i++)
        {
            if (Rcvs[i] <= 0) continue;

            n_R[j]     = Rcvs[i];
            who_R[j]   = i;
            count_R[j] = N_per_cpu[i];

            j++;
        }

        for (int i = 0, j = 0, idx = 0; i < NProcs; i++)
        {
            if (Snds[i] <= 0) continue;

            n_S[j]   = Snds[i];
            who_S[j] = i;

            int N = 0;
            for (int k = 0; k < Snds[i]; k++)
            {
                const ParallelDescriptor::CommData& cd = cd_others_need[idx+k];

                N += cd.box().numPts()*cd.nComp();
            }
            BL_ASSERT(N < std::numeric_limits<int>::max());

            count_S[j] = N;
            idx       += Snds[i];

            j++;
        }

        Rcvs.clear(); Snds.clear(); N_per_cpu.clear();
        //
        // Allocate data for rcvs as one big chunk.
        //
        value_type* recv_data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(Total_Rcvs*sizeof(value_type)));
        //
        // Post receives.
        //
        const int seqno = ParallelDescriptor::SeqNum();

        for (int i = 0, idx = 0; i < N_rcvs; i++)
        {
            BL_ASSERT(count_R[i] > 0);
            BL_ASSERT(count_R[i] < std::numeric_limits<int>::max());

            recv_reqs[i] = ParallelDescriptor::Arecv(&recv_data[idx],count_R[i],who_R[i],seqno).req();
            roffset[i]   = idx;
            idx         += count_R[i];
        }
        //
        // Send the FAB data.
        //
        FAB fab;

        for (int i = 0, idx = 0; i < N_snds; i++)
        {
            const int   N    = count_S[i];
            const int   cpu  = who_S[i];
            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
            value_type* dptr = data;

            for (int j = 0; j < n_S[i]; j++)
            {
                const ParallelDescriptor::CommData& cd = cd_others_need[idx+j];
                BL_ASSERT(cd.id() == 0);
                BL_ASSERT(cd.fromproc() == cpu);
                fab.resize(cd.box(), cd.nComp());
                const int count = fab.box().numPts()*fab.nComp();
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

            idx += n_S[i];

            if (FabArrayBase::do_async_sends)
            {
                send_data[i] = data;
                send_reqs[i] = ParallelDescriptor::Asend(data,N,cpu,seqno).req();
            }
            else
            {
                ParallelDescriptor::Send(data,N,cpu,seqno);
                BoxLib::The_Arena()->free(data);
            }
        }
        //
        // Now receive and unpack FAB data as it becomes available.
        //
        std::pair<FCDMapIter,FCDMapIter> match;

        for (int NWaits = N_rcvs, completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(recv_reqs, completed, index, status);

            for (int k = 0; k < completed; k++)
            {
                const int    who       = index[k];
                int          Processed = 0;
                value_type*  dptr      = &recv_data[roffset[who]];

                BL_ASSERT(dptr != 0);

                for (std::deque<FabArrayBase::FabComTag>::const_iterator it = fabComTagList.begin(),
                         End = fabComTagList.end();
                     it != End;
                     ++it)
                {
                    if (it->procThatHasData != who_R[who]) continue;

                    const FabArrayBase::FabComTag& tag = *it;

                    match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                    for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
                    {
                        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                        BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                        if (fcdp->subBox == tag.box)
                        {
                            value_type* dataPtr = fcdp->localFabSource->dataPtr();
                            BL_ASSERT(dataPtr != 0);
                            BL_ASSERT(fcdp->localFabSource->box() == tag.box);
                            const int N = tag.box.numPts()*tag.nComp;
                            memcpy(dataPtr, dptr, N*sizeof(value_type));
                            dptr += N;
                            Processed++;
                            break;
                        }
                    }
                }

                BL_ASSERT(Processed == n_R[who]);
            }
        }

        BoxLib::The_Arena()->free(recv_data);

        if (FabArrayBase::do_async_sends && N_snds > 0)
        {
            //
            // Now grok the asynchronous send buffers.
            //
            status.resize(N_snds);

            BL_MPI_REQUIRE( MPI_Waitall(N_snds,
                                        send_reqs.dataPtr(),
                                        status.dataPtr()) );
            //
            // And free up send buffer space.
            //
            for (int i = 0; i < N_snds; i++)
                BoxLib::The_Arena()->free(send_data[i]);
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

    BL_ASSERT(fmi != fabCopyDescList[faid.Id()].end());

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
    if ( n_grow <= 0 ) return;

    FabArrayCopyDescriptor<FAB> mfcd;
    FabArrayBase::SI&           si   = FabArrayBase::TheFBsirec(scomp,ncomp,*this);
    const FabArrayId            mfid = mfcd.RegisterFabArray(this);
    //
    // Add boxes we need to collect.
    //
    for (SI::SIRecVector::iterator it = si.m_sirec.begin(),
             End = si.m_sirec.end();
         it != End;
         ++it)
    {
        it->m_fbid = mfcd.AddBox(mfid,
                                 it->m_bx,
                                 0,
                                 it->m_j,
                                 scomp,
                                 scomp,
                                 ncomp);
    }

    mfcd.CollectData(&si.m_cache,&si.m_commdata);

    for (SI::SIRecVector::const_iterator it = si.m_sirec.begin(),
             End = si.m_sirec.end();
         it != End;
         ++it)
    {
        BL_ASSERT(DistributionMap()[it->m_i] == ParallelDescriptor::MyProc());
        //
        // Directly fill the FAB.
        //
        mfcd.FillFab(mfid,it->m_fbid,(*this)[it->m_i]);
    }
}

#endif /*BL_FABARRAY_H*/
