
#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

//
// $Id: FabArray.H,v 1.117 2001-04-11 18:11:04 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <climits>
#include <cstring>
#include <map>
#include <utility>
#include <vector>
using std::multimap;
using std::pair;
using std::vector;
#else
#include <string.h>
#include <limits.h>
#include <vector.h>
#ifdef BL_OLD_STL
#include <pair.h>
#include <multimap.h>
#else
#include <utility.h>
#include <map.h>
#endif
#endif

#include <BLassert.H>
#include <aString.H>
#include <PArray.H>
#include <Array.H>

#include <Box.H>
#include <BoxLib.H>
#include <BoxArray.H>
#include <BoxDomain.H> 
#include <FArrayBox.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <RunStats.H>
#include <Tracer.H>
#include <ccse-mpi.H>

#ifdef BL3_PROFILING
#include <BoxLib3/Profiler.H>
#endif

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

class FillBoxId
{
  public:

    FillBoxId ()
        :
        m_fillBoxId(-1),
        m_fabIndex(-1)
    {}

    FillBoxId (int        newid,
               const Box& fillbox)
        :
        m_fillBox(fillbox),
        m_fillBoxId(newid),
        m_fabIndex(-1)
    {}

    FillBoxId (const FillBoxId& rhs)
        :
        m_fillBox(rhs.m_fillBox),
        m_fillBoxId(rhs.m_fillBoxId),
        m_fabIndex(rhs.m_fabIndex)
    {}

    FillBoxId& operator= (const FillBoxId& rhs)
    {
        if (!(this == &rhs))
        {
            m_fillBox   = rhs.m_fillBox;
            m_fillBoxId = rhs.m_fillBoxId;
            m_fabIndex  = rhs.m_fabIndex;
        }
        return *this;
    }

    ~FillBoxId() {}

    int        Id () const             { return m_fillBoxId;    }
    int        FabIndex () const       { return m_fabIndex;     }
    void       FabIndex (int fabindex) { m_fabIndex = fabindex; }
    const Box& box () const            { return m_fillBox;      }

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

    CommDataCache ()
        :
        m_valid(false)
    {}

    CommDataCache (const CommDataCache& rhs)
        :
        m_commdata(rhs.m_commdata),
        m_valid(rhs.m_valid)
    {}

    void operator= (const CommDataCache& rhs)
    {
        if (this != &rhs)
        {
            m_commdata = rhs.m_commdata;
            m_valid    = rhs.m_valid;
        }
    }

    void operator= (const Array<CommData>& rhs)
    {
        m_commdata = rhs;
        m_valid    = true;
    }

    bool isValid () const { return m_valid; }

    const Array<CommData>& theCommData () const { return m_commdata; }

private:

    Array<CommData> m_commdata;
    bool            m_valid;
};

inline
ostream&
operator<< (ostream&               os,
            const Array<CommData>& cd)
{
    for (int i = 0; i < cd.length(); i++)
        os << cd[i] << '\n';
    return os;
}

//
// A forward declaration.
//
template <class T, class FAB> class FabArrayCopyDescriptor;

/*@ManDoc: An enumumeration that controls whether or not the memory for a FAB
           will actually be allocated on construction of a FabArray.
           Possible values are: Fab\_noallocate and Fab\_allocate.
*/

enum FabAlloc
{
    Fab_noallocate = 0,
    Fab_allocate
};

//
//@Man:
//@Memo: A Collection of Fortran Array-like Objects
/*@Doc:

  The FabArray<T,FAB> class implements a collection (stored as an array) of
  Fortran array-like objects.  The parameterized type FAB is intended to be
  any class derived from BaseFab<T>.  For example, FAB may be a BaseFab of
  integers, so we could write:

    FabArray< int,BaseFab<int> > int\_fabs;

  Then int\_fabs is a FabArray that can hold a collection of BaseFab<int>
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


template <class T, class FAB>
class FabArray
{

public:
    //
    //@ManDoc: Constructs an empty FabArray<T,FAB>.
    //
    FabArray ();

    /*@ManDoc: Construct a FabArray<T,FAB> with a valid region defined by bxs
               and a region of definition defined by the grow factor ngrow.
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
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);

    /*@ManDoc: This constructor is similar to constructor number two except
               that it also specifies the number of components for each FAB.
               The FAB constructor used here will take both a Box and nvar
               (the number of components) as arguments.
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
                 int             ngrow,
                 FabAlloc        mem_mode);

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

    //
    //@ManDoc: Returns the grow factor that defines the region of definition.
    //
    int nGrow () const;

    /*@ManDoc: Returns a constant reference to the BoxArray that defines the
               valid region associated with this FabArray.
    */
    const BoxArray& boxArray () const;

    /*@ManDoc: Returns a constant reference to the Kth Box in the BoxArray.
               That is, the valid region of the Kth grid.
    */
    virtual const Box& box (int K) const;

    /*@ManDoc: Returns the Kth FABs Box in the FabArray.
               That is, the region the Kth fab is actually defined on.
    */
    virtual Box fabbox (int K) const;
    //
    //@ManDoc: Returns the length of this FabArray; i.e. the number of FABs.
    //
    int length () const;

    /*@ManDoc: Returns a constant reference to the FAB associated with the
               Kth element.
    */
    const FAB& operator[] (int K) const;
    //
    //@ManDoc: Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (int K);
    //
    //@ManDoc: Returns number of variables associated with each point (nvar).
    //
    int nComp () const;
    //
    //@ManDoc: Explicitly set the Kth FAB in the FabArray to point to elem.
    //
    void setFab (int  K,
                 FAB* elem);
    //
    //@ManDoc: Releases FAB memory in the FabArray.
    //
    void clear ();
    //
    //@ManDoc: Set all components in the valid region of each FAB to val.
    //
    void setVal (T val);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp to val.
               Also set the value of nghost boundary cells.
    */
    void setVal (T   val,
                 int comp,
                 int num_comp,
                 int nghost = 0);

    /*@ManDoc: Set the value of num\_comp components in the valid region of
               each FAB in the FabArray, starting at component comp, as well
               as nghost boundary cells, to val, provided they also intersect
               with the Box region.
    */
    void setVal (T          val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells.
    */
    void setVal (T   val,
                 int nghost);

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val, including nghost boundary cells, that also
               intersect the Box region.
    */
    void setVal (T          val,
                 const Box& region,
                 int        nghost);
    //
    //@ManDoc: Set all values in the boundary region to val.
    //
    void setBndry (T val);

    /*@ManDoc: Set ncomp values in the boundary region, starting at
               start\_comp to val.
    */
    void setBndry (T   val,
                   int strt_comp,
                   int ncomp);

    /*@ManDoc: This function copies data from fa to this FabArray.  Each FAB
               in fa is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the valid region of each FAB.
    */
    void copy (const FabArray<T,FAB>& fa);

    /*@ManDoc: This function copies data from src to this FabArray.  Each FAB
               in src is intersected with all FABs in this FabArray and a copy
               is performed on the region of intersection.  The intersection
               is restricted to the num\_comp components starting at src\_comp
               in the FabArray src, with the destination components in this
               FabArray starting at dest\_comp. This assumes that
               the source and destination FabArray have identical valid
               regions.
    */
    void copy (const FabArray<T,FAB>& src,
               int                    src_comp,
               int                    dest_comp,
               int                    num_comp);

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
    //
    //@ManDoc: Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const
    {
        return distributionMap;
    }
protected:
    //
    // The data ...
    //
    mutable BoxArray    boxarray;   /* So FabSets can modify'm */
    PArray<FAB>         fabparray;
    DistributionMapping distributionMap;
    int                 n_grow;
    int                 n_comp;

private:
    //
    // These are disallowed.
    //
    FabArray (const FabArray<T,FAB>&);
    FabArray<T,FAB>& operator= (const FabArray<T,FAB>&);
    //
    // This is used locally in all define functions.
    //
    void AllocFabs ();
};

template <class T, class FAB> class ConstFabArrayIterator;

//
//@Man:
//@Memo: An iterator over the FABs in a FabArray.
/*@Doc:

  This class encapsulates an iterator over the FABs in a FabArray such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the FabArray that are actually owned by that CPU.
*/

template <class T, class FAB>
class FabArrayIterator
{
  public:
    //
    //@ManDoc: Construct a FabArrayIterator.
    //
    FabArrayIterator(FabArray<T,FAB>& fabarray);
    //
    //@ManDoc: Returns a reference to the FAB at which we currently point.
    //
    FAB& operator() ();
    //
    //@ManDoc: Returns a constant reference to the FAB at which we currently point.
    //
    const FAB& operator() () const;
    //
    //@ManDoc: Returns a reference to the FAB at which we currently point.
    //
    FAB& operator* ();
    //
    //@ManDoc: Returns a constant reference to the FAB at which we currently point.
    //
    const FAB& operator* () const;
    //
    //@ManDoc: Returns a pointer to the FAB at which we currently point.
    //
    FAB* operator-> ();
    //
    //@ManDoc: Returns a constant pointer to the FAB at which we currently point.
    //
    const FAB* operator-> () const;
    //
    //@ManDoc: Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    //@ManDoc: Returns the Box of the FAB at which we currently point.
    //
    const Box& fabbox () const;
    //
    //@ManDoc: Increments iterator to the next FAB in the FabArray that we own.
    //
    FabArrayIterator<T,FAB>& operator++ ();
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
    const FabArray<T,FAB>& theFabArray () const { return fabArray; }

  protected:

    FabArray<T,FAB>& fabArray;
    int              currentIndex;

  private:
    //
    // These are disallowed.
    //
    FabArrayIterator ();
    FabArrayIterator (const FabArrayIterator<T,FAB>&);
    FabArrayIterator<T,FAB>& operator= (const FabArrayIterator<T,FAB>&);
};

//
//@Man:
//@Memo: An iterator over a dependent FabArray.
/*@Doc:

  This class is used to iterator over a FabArray (the dependent FabArray) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a FabArrayIterator.  This allows us to iterate
  over a FabArray and access the associated FAB(s) in one or more dependent
  FabArrays.
*/

template <class T, class FAB>
class DependentFabArrayIterator
{
  public:
    //
    //@ManDoc: Constructor over a mutable dependent FabArray.
    //
    DependentFabArrayIterator(FabArrayIterator<T,FAB>& controllerfaiter,
                              FabArray<T,FAB>&         dependentfabarray);
    //
    //@ManDoc: Constructor over an immutable dependent FabArray.
    //
    DependentFabArrayIterator(FabArrayIterator<T,FAB>& controllerfaiter,
                              const FabArray<T,FAB>&   dependentfabarray);
    //
    //@ManDoc: Returns a reference to the FAB at which we currently point.
    //
    FAB& operator() ();
    //
    //@ManDoc: Returns constant reference to the FAB at which we currently point.
    //
    const FAB& operator() () const;
    //
    //@ManDoc: Returns a reference to the FAB at which we currently point.
    //
    FAB& operator* ();
    //
    //@ManDoc: Returns constant reference to FAB at which we currently point.
    //
    const FAB& operator* () const;
    //
    //@ManDoc: Returns a pointer to the FAB at which we currently point.
    //
    FAB* operator-> ();
    //
    //@ManDoc: Returns a constant pointer to the FAB at which we currently point.
    //
    const FAB* operator-> () const;
    //
    //@ManDoc: Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    //@ManDoc: Returns the Box of the FAB at which we currently point.
    //
    const Box& fabbox () const;
    //
    //@ManDoc: The index into the underlying BoxArray of the current FAB.
    //
    int index () const;

  protected:

    FabArray<T,FAB>&         fabArray;
    FabArrayIterator<T,FAB>* controller;

  private:
    //
    // These are disallowed.
    //
    DependentFabArrayIterator ();
    DependentFabArrayIterator (const DependentFabArrayIterator<T,FAB> &);
    DependentFabArrayIterator<T,FAB>& operator=(const DependentFabArrayIterator<T,FAB> &);
};

//
//@Man:
//@Memo: A constant iterator over the FABs in a FabArray.
/*@Doc:

  This class encapsulates an iterator over the FABs in a FabArray such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the FabArray that are actually owned by that CPU.  It is
  a constant iterator; you can access the FABs but you cannot change them.
*/

template <class T, class FAB>
class ConstFabArrayIterator
{
  public:
    //
    //@ManDoc: Construct a ConstFabArrayIterator.
    //
    ConstFabArrayIterator (const FabArray<T,FAB>& fabarray);
    //
    //@ManDoc: Returns a constant reference to the FAB at which we currently point.
    //
    const FAB& operator() () const;
    //
    //@ManDoc: Returns a constant reference to the FAB at which we currently point.
    //
    const FAB* operator-> () const;
    //
    //@ManDoc: Returns a constant pointer to the FAB at which we currently point.
    //
    const FAB& operator* () const;
    //
    //@ManDoc: Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    //@ManDoc: Returns the Box of the FAB at which we currently point.
    //
    const Box& fabbox () const;
    //
    //@ManDoc: Increments iterator to the next FAB in the FabArray that we own.
    //
    ConstFabArrayIterator<T,FAB>& operator++ ();
    //
    //@ManDoc: Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid () const;
    //
    //@ManDoc: The index into the underlying BoxArray of the current FAB.
    //
    int index () const;
    //
    //@ManDoc: Constant reference to FabArray over which we're iterating.
    //
    const FabArray<T,FAB>& theFabArray () const { return fabArray; }

  protected:

    const FabArray<T,FAB>& fabArray;
    int                    currentIndex;

  private:
    //
    // These are disallowed.
    //
    ConstFabArrayIterator ();
    ConstFabArrayIterator (const ConstFabArrayIterator<T,FAB> &);
    ConstFabArrayIterator<T,FAB>& operator= (const ConstFabArrayIterator<T,FAB> &);
};

//
//@Man:
//@Memo: A constant iterator over a dependent FabArray.
/*@Doc:

  This class is used to iterator over a FabArray (the dependent FabArray) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a ConstFabArrayIterator.  This allows us to
  iterate over a constant FabArray and access the associated FAB(s) in one or
  more dependent FabArrays.
*/

template <class T, class FAB>
class ConstDependentFabArrayIterator
{
  public:
    //
    //@ManDoc: The constructor.
    //
    ConstDependentFabArrayIterator(ConstFabArrayIterator<T,FAB>& cfai,
                                   const FabArray<T,FAB>&        fa);
    //
    //@ManDoc: Returns constant reference to the FAB at which we currently point.
    //
    const FAB& operator() () const;
    //
    //@ManDoc: Returns constant reference to the FAB at which we currently point.
    //
    const FAB& operator* () const;
    //
    //@ManDoc: Returns a constant pointer to the FAB at which we currently point.
    //
    const FAB* operator-> () const;
    //
    //@ManDoc: Returns the Box at the current index in the underlying BoxArray.
    //
    const Box& validbox () const;
    //
    //@ManDoc: Returns the Box of the FAB at which we currently point.
    //
    const Box& fabbox () const;
    //
    //@ManDoc: The index into the underlying BoxArray of the current FAB.
    //
    int index () const;

  protected:

    const FabArray<T,FAB>&        fabArray;
    ConstFabArrayIterator<T,FAB>* controller;

  private:
    //
    // These are disallowed.
    //
    ConstDependentFabArrayIterator ();
    ConstDependentFabArrayIterator (const ConstDependentFabArrayIterator<T,FAB> &);
    ConstDependentFabArrayIterator<T,FAB>& operator= (const ConstDependentFabArrayIterator<T,FAB> &);
};

class FabArrayId
{
  public:

    FabArrayId (int newid = -1) : fabArrayId(newid) {}

    FabArrayId (const FabArrayId& rhs) : fabArrayId(rhs.fabArrayId) {}

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

template <class T, class FAB>
class FabArrayCopyDescriptor
{
  typedef multimap<int,FabCopyDescriptor<FAB>*> FCDMap;
  typedef typename FCDMap::value_type           FCDMapValueType;
  typedef typename FCDMap::iterator             FCDMapIter;
  typedef typename FCDMap::const_iterator       FCDMapConstIter;

  public:

    FabArrayCopyDescriptor ();

    ~FabArrayCopyDescriptor ();

    FabArrayId RegisterFabArray(FabArray<T,FAB> *fabarray);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      bool       bUseValidBox = true);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp,
                      bool       bUseValidBox = true);
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
                     BoxDomain& unfilledBoxDomain,
                     BoxList&   filledBoxes);

    vector<FabArray<T,FAB>*> fabArrays;
    vector<FCDMap>           fabCopyDescList;
    vector<FabComTag>        fabComTagList;
    int                      nextFillBoxId;
    bool                     dataAvailable;

  private:
    //
    // These are disallowed.
    //
    FabArrayCopyDescriptor (const FabArrayCopyDescriptor<T,FAB>&);
    FabArrayCopyDescriptor<T,FAB>& operator= (const FabArrayCopyDescriptor<T,FAB> &);
};

//
// Inlines.
//

template <class T, class FAB>
inline
int
FabArray<T,FAB>::nGrow () const
{
    return n_grow;
}

template <class T, class FAB>
inline
const BoxArray&
FabArray<T,FAB>::boxArray () const
{
    return boxarray;
}

template <class T, class FAB>
inline
const Box&
FabArray<T,FAB>::box (int K) const
{
    return boxarray[K];
}

template <class T, class FAB>
inline
Box
FabArray<T,FAB>::fabbox (int K) const
{
    //
    // Do not use fabparray[K] because it may not be valid in parallel.
    //
#ifndef BL_NAMESPACE
    return ::grow(boxarray[K], n_grow);
#else
    return BL_NAMESPACE::grow(boxarray[K], n_grow);
#endif
}

template <class T, class FAB>
inline
int
FabArray<T,FAB>::length () const
{
    return boxarray.length();
}

template <class T, class FAB>
inline
const FAB&
FabArray<T,FAB>::operator[] (int K) const
{
    return fabparray[K];
}

template <class T, class FAB>
inline
FAB&
FabArray<T,FAB>::operator[] (int K)
{
    return fabparray[K];
}

template <class T, class FAB>
inline
int
FabArray<T,FAB>::nComp () const
{
    return n_comp;
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clear ()
{
    fabparray.clear();
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setVal (T   val,
                         int nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setVal (T          val,
                         const Box& region,
                         int        nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

template <class T, class FAB>
FabArray<T,FAB>::FabArray ()
    :
    fabparray(0,PArrayManage),
    n_grow(0),
    n_comp(0)
{}

template <class T, class FAB>
FabArray<T,FAB>::FabArray (const BoxArray& bxs,
                           int             ngrow,
                           FabAlloc        alloc)
    :
    fabparray(0,PArrayManage)
{
    define(bxs,ngrow,alloc);
}

template <class T, class FAB>
FabArray<T,FAB>::FabArray (const BoxArray& bxs,
                           int             nvar,
                           int             ngrow,
                           FabAlloc        alloc)
    :
    fabparray(0, PArrayManage)
{
    define(bxs,nvar,ngrow,alloc);
}

template <class T, class FAB>
FabArray<T,FAB>::~FabArray ()
{}

template <class T, class FAB>
bool
FabArray<T,FAB>::ok () const
{
    BL_ASSERT(boxarray.ready());

    long isok = 1;

    ConstFabArrayIterator<T,FAB> fai(*this);

    for ( ; fai.isValid() && isok; ++fai)
    {
        if (fabparray.defined(fai.index()))
        {
#ifndef BL_NAMESPACE
            if (fai().box() != ::grow(box(fai.index()),n_grow))
                isok = 0;
#else
            if (fai().box() != BL_NAMESPACE::grow(box(fai.index()),n_grow))
                isok = 0;
#endif
        }
        else
        {
            isok = 0;
        }
    }

    ParallelDescriptor::ReduceLongAnd(isok);

    return isok != 0;
}

template <class T, class FAB>
void
FabArray<T,FAB>::define (const BoxArray& bxs,
                         int             ngrow,
                         FabAlloc        alloc)
{
    BL_ASSERT(!boxarray.ready());
    n_grow = ngrow;
    n_comp = 0;
    boxarray.define(bxs);
    distributionMap.define(boxarray,ParallelDescriptor::NProcsCFD());
    fabparray.resize(boxarray.length());
    if (alloc == Fab_allocate)
    {
        n_comp = 1;
        AllocFabs();
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::define (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc)
{
    BL_ASSERT(!boxarray.ready());
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap.define(boxarray,ParallelDescriptor::NProcsCFD());
    fabparray.resize(bxs.length());
    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class T, class FAB>
void
FabArray<T,FAB>::define (const BoxArray&            bxs,
                         int                        nvar,
                         int                        ngrow,
			 const DistributionMapping& dm,
                         FabAlloc                   alloc)
{
    BL_ASSERT(!boxarray.ready());
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap = dm;
    fabparray.resize(bxs.length());
    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class T, class FAB>
void
FabArray<T,FAB>::AllocFabs ()
{
    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
#ifndef BL_NAMESPACE
        Box tmp = ::grow(fai.validbox(), n_grow);
#else
        Box tmp = BL_NAMESPACE::grow(fai.validbox(), n_grow);
#endif
        fabparray.set(fai.index(), new FAB(tmp, n_comp));
        RunStats::addNumPts(tmp.numPts());
    }
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setFab (int  boxno,
                         FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();

    BL_ASSERT(n_comp == elem->nComp());
    BL_ASSERT(boxarray.ready());
#ifndef BL_NAMESPACE
    BL_ASSERT(elem->box() == ::grow(boxarray[boxno],n_grow));
#else
    BL_ASSERT(elem->box() == BL_NAMESPACE::grow(boxarray[boxno],n_grow));
#endif
    BL_ASSERT(!fabparray.defined(boxno));
    BL_ASSERT(distributionMap[boxno] == ParallelDescriptor::MyProc());

    fabparray.set(boxno,elem);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::setBndry (T val)
{
    setBndry(val, 0, n_comp);
}

template <class T, class FAB>
void
FabArray<T,FAB>::setBndry (T   val,
                           int strt_comp,
                           int ncomp)
{
    if (n_grow > 0)
    {
        for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
        {
            fai().setComplement(val, fai.validbox(), strt_comp, ncomp);
        }
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::copy (const FabArray<T,FAB>& src,
                       int                    scomp,
                       int                    dcomp,
                       int                    ncomp)
{
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::copy()");
#endif

#ifdef BL_USE_MPI

    if (boxarray == src.boxarray)
    {
        BL_ASSERT(distributionMap == src.distributionMap);

        for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
        {
            DependentFabArrayIterator<T,FAB> dfai(fai,src);
            BL_ASSERT(fai.validbox() == dfai.validbox());
            const Box& bx = fai.validbox();
            fai().copy(dfai(),bx,scomp,bx,dcomp,ncomp);
        }

        return;
    }

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    static Array<T*>          fab_data(NProcs);
    static Array<int>         indx(NProcs);
    static Array<MPI_Status>  status(NProcs);
    static Array<MPI_Request> reqs(NProcs);

    Array< vector<FabComTag> > SndTags(NProcs);
    Array< vector<FabComTag> > RcvTags(NProcs);

    FabComTag tag;

    for (int i = 0; i < length(); i++)
    {
        if (distributionMap[i] == MyProc)
        {
            for (int ii = 0; ii < src.boxarray.length(); ii++)
            {
                if (src.boxarray[ii].intersects(boxarray[i]))
                {
                    Box bx = src.boxarray[ii] & boxarray[i];

                    if (src.distributionMap[ii] == MyProc)
                    {
                        fabparray[i].copy(src[ii],bx,scomp,bx,dcomp,ncomp);
                    }
                    else
                    {
                        tag.box      = bx;
                        tag.fabIndex = i;

                        RcvTags[src.distributionMap[ii]].push_back(tag);
                    }
                }
            }
        }
        else
        {
            for (int ii = 0; ii < src.boxarray.length(); ii++)
            {
                if (src.distributionMap[ii] == MyProc)
                {
                    if (src.boxarray[ii].intersects(boxarray[i]))
                    {
                        tag.box      = src.boxarray[ii] & boxarray[i];
                        tag.fabIndex = ii;

                        SndTags[distributionMap[i]].push_back(tag);
                    }
                }
            }
        }
    }

    if (NProcs == 1) return;

    static RunStats stats("fabarray_copy");

    stats.start();

    const int seqno = ParallelDescriptor::SeqNum();

    int rc, NWaits = 0;
    //
    // Post one receive for each chunk being sent by other CPUs.
    //
    for (int i = 0; i < NProcs; i++)
    {
        reqs[i] = MPI_REQUEST_NULL;

        if (!RcvTags[i].empty())
        {
            NWaits++;

            size_t N = 0;

            for (int j = 0; j < RcvTags[i].size(); j++)
                N += RcvTags[i][j].box.numPts() * ncomp;

            BL_ASSERT(N < INT_MAX);

            fab_data[i] = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));

            if ((rc = MPI_Irecv(fab_data[i],
                                int(N),
                                mpi_data_type(fab_data[i]),
                                i,
                                seqno,
                                ParallelDescriptor::Communicator(),
                                &reqs[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
    }

    FAB fab;
    //
    // Send the agglomerated FAB data.
    //
    for (int i = 0; i < NProcs; i++)
    {
        if (!SndTags[i].empty())
        {
            size_t N = 0;

            for (int j = 0; j < SndTags[i].size(); j++)
                    N += SndTags[i][j].box.numPts() * ncomp;

            BL_ASSERT(N < INT_MAX);

            T* data = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));
            T* dptr = data;

            for (int j = 0; j < SndTags[i].size(); j++)
            {
                const Box& bx = SndTags[i][j].box;
                fab.resize(bx, ncomp);
                fab.copy(src[SndTags[i][j].fabIndex],bx,scomp,bx,0,ncomp);
                int count = bx.numPts() * ncomp;
                memcpy(dptr, fab.dataPtr(), count*sizeof(T));
                dptr += count;
            }

            BL_ASSERT(data+N == dptr);

            if ((rc = MPI_Send(data,
                               int(N),
                               mpi_data_type(data),
                               i,
                               seqno,
                               ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            The_FAB_Arena->free(data);
        }
    }
    //
    // Now receive and unpack FAB data.
    //
    for (int completed; NWaits > 0; NWaits -= completed)
    {
        if ((rc = MPI_Waitsome(NProcs,
                               reqs.dataPtr(),
                               &completed,
                               indx.dataPtr(),
                               status.dataPtr())) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);

        for (int k = 0; k < completed; k++)
        {
            T* dptr = fab_data[indx[k]];

            BL_ASSERT(!(dptr == 0));

            for (int j = 0; j < RcvTags[indx[k]].size(); j++)
            {
                const Box& bx = RcvTags[indx[k]][j].box;
                fab.resize(bx, ncomp);
                int N = bx.numPts() * ncomp;
                memcpy(fab.dataPtr(), dptr, N*sizeof(T));
                fabparray[RcvTags[indx[k]][j].fabIndex].copy(fab,bx,0,bx,dcomp,ncomp);
                dptr += N;
            }

            The_FAB_Arena->free(fab_data[indx[k]]);
        }
    }

    stats.end();

#else

    if (boxarray == src.boxarray)
    {
        BL_ASSERT(distributionMap == src.distributionMap);

        for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
        {
            DependentFabArrayIterator<T,FAB> dfai(fai,src);
            BL_ASSERT(fai.validbox() == dfai.validbox());
            const Box& bx = fai.validbox();
            fai().copy(dfai(),bx,scomp,bx,dcomp,ncomp);
        }
    }
    else
    {
        for (int i = 0; i < length(); ++i)
        {
            for (int ii = 0; ii < src.boxarray.length(); ++ii)
            {
                if (src.boxarray[ii].intersects(boxarray[i]))
                {
                    Box bx = src.boxarray[ii] & boxarray[i];

                    fabparray[i].copy(src[ii],bx,scomp,bx,dcomp,ncomp);
                }
            }
        }
    }
#endif /*BL_USE_MPI*/
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::copy (const FabArray<T,FAB>& src)
{
    copy(src,0,0,nComp());
}

//
// Copies to FABs, note that destination is first arg.
//

template <class T, class FAB>
inline
void
FabArray<T,FAB>::copy (FAB& dest) const
{
    copy(dest, dest.box(), 0, 0, dest.nComp());
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::copy (FAB&       dest,
                       const Box& subbox) const
{
    copy(dest, subbox, 0, 0, dest.nComp());
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::copy (FAB& dest,
                       int  scomp,
                       int  dcomp,
                       int  ncomp) const
{
    copy(dest, dest.box(), scomp, dcomp, ncomp);
}

template <class T, class FAB>
void
FabArray<T,FAB>::copy (FAB&       dest,
                       const Box& subbox,
                       int        scomp,
                       int        dcomp,
                       int        ncomp) const
{
#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::copy(Fab&, ...)");
#endif
    BL_ASSERT(dcomp + ncomp <= dest.nComp());

    if (ParallelDescriptor::NProcs() == 1)
    {
        for (int j = 0; j < length(); ++j)
        {
            if (boxarray[j].intersects(subbox))
            {
                Box destbox = boxarray[j] & subbox;

                dest.copy(fabparray[j],destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

#ifdef BL_USE_MPI
    FArrayBox overlapData;

    for (int i = 0; i < length(); i++)
    {
        if (subbox.intersects(boxarray[i]))
        {
            Box bx = subbox & boxarray[i];

            overlapData.resize(bx,ncomp);

            if (ParallelDescriptor::MyProc() == distributionMap[i])
                overlapData.copy(fabparray[i],bx,scomp,bx,0,ncomp);

            int rc = MPI_Bcast(overlapData.dataPtr(),
                               bx.numPts()*ncomp,
                               mpi_data_type(overlapData.dataPtr()),
                               distributionMap[i],
                               ParallelDescriptor::Communicator());

            if (rc != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            dest.copy(overlapData,bx,0,bx,dcomp,ncomp);
        }
    }
#endif /*BL_USE_MPI*/
}

template <class T, class FAB>
void
FabArray<T,FAB>::setVal (T val)
{
    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        fai().setVal(val);
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::setVal (T   val,
                         int comp,
                         int ncomp,
                         int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
#ifndef BL_NAMESPACE
        fai().setVal(val,::grow(fai.validbox(),nghost), comp, ncomp);
#else
        fai().setVal(val, BL_NAMESPACE::grow(fai.validbox(),nghost), comp, ncomp);
#endif
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::setVal (T          val,
                         const Box& region,
                         int        comp,
                         int        ncomp,
                         int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(fai.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(fai.validbox(),nghost) & region;
#endif

        if (b.ok())
            fai().setVal(val, b, comp, ncomp);
    }
}

template <class T, class FAB>
inline
FabArrayIterator<T,FAB>::FabArrayIterator (FabArray<T,FAB>& fabarray)
    :
    fabArray(fabarray),
    currentIndex(0)
{
#ifdef BL_USE_MPI
    //
    // Increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    while (fabArray.DistributionMap()[currentIndex] != MyProc)
    {
        ++currentIndex;
    }
#endif
}

template <class T, class FAB>
inline
FAB&
FabArrayIterator<T,FAB>::operator() ()
{
  BL_ASSERT(fabArray.DistributionMap()[currentIndex] == ParallelDescriptor::MyProc());

  return fabArray[currentIndex];
}

template <class T, class FAB>
inline
FAB&
FabArrayIterator<T,FAB>::operator* ()
{
    return FabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
const FAB&
FabArrayIterator<T,FAB>::operator* () const
{
    return FabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
FAB*
FabArrayIterator<T,FAB>::operator-> ()
{
    return &FabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
const FAB*
FabArrayIterator<T,FAB>::operator-> () const
{
    return &FabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
const FAB&
FabArrayIterator<T,FAB>::operator() () const
{
  BL_ASSERT(fabArray.DistributionMap()[currentIndex] == ParallelDescriptor::MyProc());

  return fabArray[currentIndex];
}

template <class T, class FAB>
inline
FabArrayIterator<T,FAB>&
FabArrayIterator<T,FAB>::operator++ ()
{
#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (fabArray.DistributionMap()[currentIndex] != MyProc);
#else
    ++currentIndex;
#endif

    return *this;
}

template <class T, class FAB>
inline
int
FabArrayIterator<T,FAB>::index () const
{
    return currentIndex;
}

template <class T, class FAB>
inline
const Box&
FabArrayIterator<T,FAB>::validbox () const
{
    return fabArray.box(currentIndex);
}

template <class T, class FAB>
inline
const Box&
FabArrayIterator<T,FAB>::fabbox () const
{
    return fabArray[currentIndex].box();
}

template <class T, class FAB>
inline
bool
FabArrayIterator<T,FAB>::isValid ()
{
    BL_ASSERT(currentIndex >= 0);

#if defined(BL_USE_MPI) && defined(HG_DEBUG)
    extern bool HG_is_debugging;
    if (HG_is_debugging)
    {
	if (currentIndex < fabArray.length()) return true;
	int rc = MPI_Barrier(ParallelDescriptor::Communicator());
	if (rc != MPI_SUCCESS) ParallelDescriptor::Abort(rc);
	return false;
    }
#endif
    return currentIndex < fabArray.length();
}

template <class T, class FAB>
inline
DependentFabArrayIterator<T,FAB>::DependentFabArrayIterator(
    FabArrayIterator<T,FAB>& controllerfaiter,
    FabArray<T,FAB>&         dependentfabarray)
    :
    fabArray(dependentfabarray),
    controller(&controllerfaiter)
{
    BL_ASSERT(controllerfaiter.theFabArray().DistributionMap() ==
           dependentfabarray.DistributionMap());
}

template <class T, class FAB>
inline
DependentFabArrayIterator<T,FAB>::DependentFabArrayIterator(
    FabArrayIterator<T,FAB>& controllerfaiter,
    const FabArray<T,FAB>&   dependentfabarray)
    :
    fabArray(const_cast<FabArray<T,FAB>&>(dependentfabarray)),
    controller(&controllerfaiter)
{
    BL_ASSERT(controllerfaiter.theFabArray().DistributionMap() ==
           dependentfabarray.DistributionMap());
}

template <class T, class FAB>
inline
int
DependentFabArrayIterator<T,FAB>::index () const
{
    return controller->index();
}

template <class T, class FAB>
inline
const Box&
DependentFabArrayIterator<T,FAB>::validbox() const
{
    return fabArray.box(controller->index());
}

template <class T, class FAB>
inline
const Box&
DependentFabArrayIterator<T,FAB>::fabbox() const
{
    return fabArray[controller->index()].box();
}

template <class T, class FAB>
inline
FAB&
DependentFabArrayIterator<T,FAB>::operator() ()
{
    return fabArray[controller->index()];
}

template <class T, class FAB>
inline
const FAB&
DependentFabArrayIterator<T,FAB>::operator() () const
{
    return fabArray[controller->index()];
}

template <class T, class FAB>
inline
FAB&
DependentFabArrayIterator<T,FAB>::operator* ()
{
    return DependentFabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
const FAB&
DependentFabArrayIterator<T,FAB>::operator* () const
{
    return DependentFabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
FAB*
DependentFabArrayIterator<T,FAB>::operator-> ()
{
    return &DependentFabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
const FAB*
DependentFabArrayIterator<T,FAB>::operator-> () const
{
    return &DependentFabArrayIterator<T,FAB>::operator()();
}

template <class T, class FAB>
inline
ConstFabArrayIterator<T,FAB>::ConstFabArrayIterator(const FabArray<T,FAB>& fabarray)
    :
    fabArray(fabarray),
    currentIndex(0)
{
#ifdef BL_USE_MPI
    //
    // Increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc
    //
    const int MyProc = ParallelDescriptor::MyProc();

    while (fabArray.DistributionMap()[currentIndex] != MyProc)
    {
        ++currentIndex;
    }
#endif
}

template <class T, class FAB>
inline
const FAB&
ConstFabArrayIterator<T,FAB>::operator() () const
{
    return fabArray[currentIndex];
}

template <class T, class FAB>
inline
const FAB&
ConstFabArrayIterator<T,FAB>::operator* () const
{
    return ConstFabArrayIterator<T,FAB>::operator() ();
}

template <class T, class FAB>
inline
const FAB*
ConstFabArrayIterator<T,FAB>::operator-> () const
{
    return &ConstFabArrayIterator<T,FAB>::operator() ();
}

template <class T, class FAB>
inline
ConstFabArrayIterator<T,FAB>&
ConstFabArrayIterator<T,FAB>::operator++ ()
{
#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (fabArray.DistributionMap()[currentIndex] != MyProc);
#else
    ++currentIndex;
#endif

    return *this;
}

template <class T, class FAB>
inline
int
ConstFabArrayIterator<T,FAB>::index () const
{
    return currentIndex;
}

template <class T, class FAB>
inline
const Box&
ConstFabArrayIterator<T,FAB>::validbox () const
{
    return fabArray.box(currentIndex);
}

template <class T, class FAB>
inline
const Box&
ConstFabArrayIterator<T,FAB>::fabbox() const
{
    return fabArray[currentIndex].box();
}

template <class T, class FAB>
inline
bool
ConstFabArrayIterator<T,FAB>::isValid () const
{
    BL_ASSERT(currentIndex >= 0);

#if defined(BL_USE_MPI) && defined(HG_DEBUG)
    extern bool HG_is_debugging;
    if (HG_is_debugging)
    {
	if (currentIndex < fabArray.length()) return true;
	int rc = MPI_Barrier(ParallelDescriptor::Communicator());
	if (rc != MPI_SUCCESS) ParallelDescriptor::Abort(rc);
	return false;
    }
#endif
    return currentIndex < fabArray.length();
}

template <class T, class FAB>
inline
ConstDependentFabArrayIterator<T,FAB>::ConstDependentFabArrayIterator(
    ConstFabArrayIterator<T,FAB>& controllerfaiter,
    const FabArray<T,FAB>&        dependentfabarray)
    :
    fabArray(dependentfabarray),
    controller(&controllerfaiter)
{
    BL_ASSERT(controllerfaiter.theFabArray().DistributionMap() ==
           dependentfabarray.DistributionMap());
}

template <class T, class FAB>
inline
int
ConstDependentFabArrayIterator<T,FAB>::index () const
{
    return controller->index();
}

template <class T, class FAB>
inline
const Box&
ConstDependentFabArrayIterator<T,FAB>::validbox () const
{
    return fabArray.box(controller->index());
}

template <class T, class FAB>
inline
const Box&
ConstDependentFabArrayIterator<T,FAB>::fabbox() const
{
    return fabArray[controller->index()].box();
}

template <class T, class FAB>
inline
const FAB&
ConstDependentFabArrayIterator<T,FAB>::operator()() const
{
    return fabArray[controller->index()];
}

template <class T, class FAB>
inline
const FAB&
ConstDependentFabArrayIterator<T,FAB>::operator* () const
{
    return ConstDependentFabArrayIterator<T,FAB>::operator() ();
}

template <class T, class FAB>
inline
const FAB*
ConstDependentFabArrayIterator<T,FAB>::operator-> () const
{
    return &ConstDependentFabArrayIterator<T,FAB>::operator() ();
}

template <class T, class FAB>
FabArrayCopyDescriptor<T,FAB>::FabArrayCopyDescriptor ()
    :
    nextFillBoxId(0),
    dataAvailable(false)
{}

template <class T, class FAB>
FabArrayId
FabArrayCopyDescriptor<T,FAB>::RegisterFabArray(FabArray<T,FAB>* fabarray)
{
    BL_ASSERT(fabArrays.size() == fabCopyDescList.size());

    fabArrays.push_back(fabarray);  /* Bump size() by one */

    fabCopyDescList.resize(fabArrays.size(), FCDMap());

    return FabArrayId(fabArrays.size() - 1);
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::AddBoxDoIt (FabArrayId fabarrayid,
                                           const Box& destFabBox,
                                           BoxList*   returnedUnfilledBoxes,
                                           int        faindex,
                                           int        srccomp,
                                           int        destcomp,
                                           int        numcomp,
                                           bool       bUseValidBox,
                                           BoxDomain& unfilledBoxDomain,
                                           BoxList&   filledBoxes)
{
    const int MyProc = ParallelDescriptor::MyProc();

    FabArray<T,FAB>* fabArray = fabArrays[fabarrayid.Id()];

    BL_ASSERT(faindex >= 0 && faindex < fabArray->length());

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
        filledBoxes.add(intersect);

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

template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox (FabArrayId fabarrayid,
                                       const Box& destFabBox,
                                       BoxList*   returnedUnfilledBoxes,
                                       int        srccomp,
                                       int        destcomp,
                                       int        numcomp,
                                       bool       bUseValidBox)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());
    BoxList   filledBoxes(destFabBox.ixType());

    if (!(returnedUnfilledBoxes == 0))
    {
        unfilledBoxDomain.add(destFabBox);
    }

    for (int i = 0, N = fabArrays[fabarrayid.Id()]->length(); i < N; i++)
    {
        AddBoxDoIt(fabarrayid,
                   destFabBox,
                   returnedUnfilledBoxes,
                   i,
                   srccomp,
                   destcomp,
                   numcomp,
                   bUseValidBox,
                   unfilledBoxDomain,
                   filledBoxes);
    }

    if (!(returnedUnfilledBoxes == 0))
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox (FabArrayId fabarrayid,
                                       const Box& destFabBox,
                                       BoxList*   returnedUnfilledBoxes,
                                       int        fabarrayindex,
                                       int        srccomp,
                                       int        destcomp,
                                       int        numcomp,
                                       bool       bUseValidBox)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());
    BoxList   filledBoxes(destFabBox.ixType());

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
               unfilledBoxDomain,
               filledBoxes);

    if (!(returnedUnfilledBoxes == 0))
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class T, class FAB>
inline
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox (FabArrayId fabarrayid,
                                       const Box& destFabBox,
                                       BoxList*   returnedUnfilledBoxes,
                                       bool       bUseValidBox)
{
    return AddBox(fabarrayid,
                  destFabBox,
                  returnedUnfilledBoxes,
                  0,
                  0,
                  fabArrays[fabarrayid.Id()]->nComp(),
                  bUseValidBox);
}

template <class T, class FAB>
FabArrayCopyDescriptor<T,FAB>::~FabArrayCopyDescriptor()
{
   clear();
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::clear ()
{
   for (int i = 0; i < fabCopyDescList.size(); ++i)
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

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::CollectData (Array<int>*    snd_cache,
                                            CommDataCache* cd_cache)
{
    dataAvailable = true;

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1) return;

#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::CollectData()");
#endif

#ifdef BL_USE_MPI
    static RunStats stats("collect_data");

    static Array<int>         Snds(NProcs);
    static Array<int>         Rcvs(NProcs);
    static Array<int>         indx(NProcs);
    static Array<T*>          fab_data(NProcs);
    static Array<MPI_Request> req_data(NProcs);
    static Array<MPI_Request> req_cd(NProcs);
    static Array<MPI_Status>  status(NProcs);
    static Array<CommData>    senddata;
    static Array<CommData>    recv_cd;

    stats.start();

    for (int i = 0; i < NProcs; i++)
        Snds[i] = Rcvs[i] = 0;

    for (int i = 0; i < NProcs; i++)
        req_data[i] = req_cd[i] = MPI_REQUEST_NULL;

    const int seqno_1 = ParallelDescriptor::SeqNum();
    const int seqno_2 = ParallelDescriptor::SeqNum();

    int rc, idx = 0, NumReqs = 0, NWaits;
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (int i = 0; i < fabComTagList.size(); i++)
    {
        BL_ASSERT(fabComTagList[i].box.ok());
        BL_ASSERT(fabComTagList[i].procThatNeedsData == MyProc);
        BL_ASSERT(fabComTagList[i].procThatHasData   != MyProc);

        Rcvs[fabComTagList[i].procThatHasData]++;
    }
    BL_ASSERT(Rcvs[MyProc] == 0);
    //
    // Set Snds[i] to # of blocks we must send to CPU i ...
    //
#ifdef NDEBUG
    if (snd_cache == 0 || snd_cache->length() == 0)
#endif
    {
        for (int i = 0; i < NProcs; i++)
            if ((rc = MPI_Gather(&Rcvs[i],
                                 1,
                                 MPI_INT,
                                 Snds.dataPtr(),
                                 1,
                                 MPI_INT,
                                 i,
                                 ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

        BL_ASSERT(Snds[MyProc] == 0);
    }

    if (snd_cache)
    {
        if (snd_cache->length() > 0)
        {
            BL_ASSERT(Snds == *snd_cache);

            Snds = *snd_cache;
        }
        else
            *snd_cache = Snds;
    }

    for (int i = 0; i < NProcs; i++)
        NumReqs += Snds[i];

    recv_cd.resize(NumReqs);
    //
    // Post one receive for each chunk being requested by other CPUs.
    // This is the CommData describing what FAB data needs to be sent.
    //
#ifdef NDEBUG
    if (cd_cache == 0 || !cd_cache->isValid())
#endif
    {
        //
        // Make sure we can treat CommData as a stream of integers.
        //
        BL_ASSERT(sizeof(CommData) == CommData::DIM*sizeof(int));

        for (int i = 0; i < NProcs; i++)
        {
            if (Snds[i] > 0)
            {
                if ((rc = MPI_Irecv(&recv_cd[idx],
                                    Snds[i] * CommData::DIM,
                                    MPI_INT,
                                    i,
                                    seqno_1,
                                    ParallelDescriptor::Communicator(),
                                    &req_cd[i])) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                idx += Snds[i];
            }
        }

        BL_ASSERT(idx == NumReqs);
        //
        // Make one send to each CPU from which we want data.
        // Start Send()ing to upper neighbor.
        //
        for (int k = 0, i = MyProc+1; k < NProcs; k++, i++)
        {
            i %= NProcs;

            if (Rcvs[i] > 0)
            {
                senddata.resize(Rcvs[i]);

                int Processed = 0;

                for (int j = 0; j < fabComTagList.size(); j++)
                {
                    if (fabComTagList[j].procThatHasData == i)
                    {
                        CommData data(0,
                                      fabComTagList[j].fabIndex,
                                      MyProc,
                                      0,
                                      fabComTagList[j].nComp,
                                      fabComTagList[j].srcComp,
                                      fabComTagList[j].fabArrayId,
                                      fabComTagList[j].box);

                        senddata[Processed++] = data;
                    }
                }

                BL_ASSERT(Processed == Rcvs[i]);

                if ((rc = MPI_Send(senddata.dataPtr(),
                                   senddata.length() * CommData::DIM,
                                   MPI_INT,
                                   i,
                                   seqno_1,
                                   ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);
            }
        }

        NWaits = 0;
        for (int i = 0; i < NProcs; i++)
            if (req_cd[i] != MPI_REQUEST_NULL)
                NWaits++;

        for (int completed; NWaits > 0; NWaits -= completed)
        {
            if ((rc = MPI_Waitsome(NProcs,
                                   req_cd.dataPtr(),
                                   &completed,
                                   indx.dataPtr(),
                                   status.dataPtr())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
    }

    if (cd_cache)
    {
        if (cd_cache->isValid())
        {
            BL_ASSERT(recv_cd == cd_cache->theCommData());

            recv_cd = cd_cache->theCommData();
        }
        else
            *cd_cache = recv_cd;
    }
    //
    // Post one receive for data being requested by other CPUs.
    //
    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            //
            // Got to figure out # of T's to expect from each CPU.
            //
            size_t N = 0;

            for (int j = 0; j < fabComTagList.size(); j++)
                if (fabComTagList[j].procThatHasData == i)
                    N += fabComTagList[j].box.numPts()*fabComTagList[j].nComp;

            BL_ASSERT(N < INT_MAX);

            fab_data[i] = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));

            if ((rc = MPI_Irecv(fab_data[i],
                                int(N),
                                mpi_data_type(fab_data[i]),
                                i,
                                seqno_2,
                                ParallelDescriptor::Communicator(),
                                &req_data[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
    }
    //
    // Send the agglomerated FAB data.  Start Send()ing to upper neighbor.
    //
    idx = 0;

    FAB fab;

    for (int k = 0, i = MyProc+1; k < NProcs; k++, i++)
    {
        i %= NProcs;
        
        int strt = 0;
        for (int j = 0; j < i; j++)
            strt += Snds[j];

        if (Snds[i] > 0)
        {
            size_t N = 0;

            for (int j = 0; j < Snds[i]; j++)
                N += recv_cd[strt+j].box().numPts() * recv_cd[strt+j].nComp();

            BL_ASSERT(N < INT_MAX);

            T* data = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));
            T* dptr = data;

            for (int j = 0; j < Snds[i]; j++)
            {
                const CommData& cd = recv_cd[strt+j];

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

                memcpy(dptr, fab.dataPtr(), count*sizeof(T));

                dptr += count;
            }

            BL_ASSERT(data+N == dptr);

            if ((rc = MPI_Send(data,
                               int(N),
                               mpi_data_type(data),
                               i,
                               seqno_2,
                               ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            The_FAB_Arena->free(data);

            idx += Snds[i];
        }
    }

    BL_ASSERT(idx == NumReqs);
    //
    // Now receive and unpack FAB data.
    //
    pair<FCDMapIter,FCDMapIter> match;

    NWaits = 0;
    for (int i = 0; i < NProcs; i++)
        if (req_data[i] != MPI_REQUEST_NULL)
            NWaits++;

    for (int completed; NWaits > 0; NWaits -= completed)
    {
        if ((rc = MPI_Waitsome(NProcs,
                               req_data.dataPtr(),
                               &completed,
                               indx.dataPtr(),
                               status.dataPtr())) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);

        for (int k = 0; k < completed; k++)
        {
            int Processed = 0;
            T*  dptr      = fab_data[indx[k]];

            BL_ASSERT(!(dptr == 0));

            for (int j = 0; j < fabComTagList.size(); j++)
            {
                if (fabComTagList[j].procThatHasData == indx[k])
                {
                    const FabComTag& tag = fabComTagList[j];

                    match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                    FCDMapIter fmi = match.first;

                    for ( ; fmi != match.second; ++fmi)
                    {
                        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                        BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                        if (fcdp->subBox == tag.box)
                        {
                            T* dataPtr = fcdp->localFabSource->dataPtr();
                            BL_ASSERT(!(dataPtr == 0));
                            BL_ASSERT(fcdp->localFabSource->box() == tag.box);
                            int N = tag.box.numPts()*tag.nComp;
                            memcpy(dataPtr, dptr, N*sizeof(T));
                            dptr += N;
                            Processed++;
                            break;
                        }
                    }

                    BL_ASSERT(!(fmi == match.second));
                }
            }

            BL_ASSERT(Processed == Rcvs[indx[k]]);

            The_FAB_Arena->free(fab_data[indx[k]]);

            fab_data[indx[k]] = 0;
        }
    }

    stats.end();
#endif /*BL_USE_MPI*/
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::FillFab (FabArrayId       faid,
                                        const FillBoxId& fillboxid,
                                        FAB&             destFab)
{
    BL_ASSERT(dataAvailable);

    pair<FCDMapIter,FCDMapIter> match = fabCopyDescList[faid.Id()].equal_range(fillboxid.Id());

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

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::FillFab (FabArrayId       faid,
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

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::PrintStats () const
{
    const int MyProc = ParallelDescriptor::MyProc();

    cout << "----- "
         << MyProc
         << ":  Parallel stats for FabArrayCopyDescriptor:" << '\n';

    for (int fa = 0; fa < fabArrays.size(); ++fa)
    {
        cout << "fabArrays["
             << fa
             << "]->boxArray() = "
             << fabArrays[fa]->boxArray()
             << '\n';
    }
}

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_FABARRAY_H*/
