//BL_COPYRIGHT_NOTICE

#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

//
// $Id: FabArray.H,v 1.7 1997-09-17 22:51:16 lijewski Exp $
//

#include <limits.h>

#include <Assert.H>
#include <aString.H>
#include <PArray.H>
#include <Array.H>

#include <Boolean.H>
#include <Box.H>
#include <BoxLib.H>
#include <BoxArray.H>
#include <BoxDomain.H> 
#include <FArrayBox.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <Utility.H>

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

    /*@ManDoc: Returns a constant reference to the Kth Fabs Box in the FabArray.
               That is, the region the Kth fab is actually defined on.
    */
    virtual const Box fabbox (int K) const;
    //
    //@ManDoc: Returns the length of this FabArray; i.e. the number of FABs.
    //
    int length () const;

    /*@ManDoc: Returns a constant reference to the FAB associated with the
               Kth element.
    */
//protected:
    const FAB& operator[] (int K) const;
    //
    //@ManDoc: Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (int K);

public:
    //@ManDoc: Returns number of variables associated with each point (nvar).
    //
    int nComp () const;
    //
    //@ManDoc: Explicitly set the Kth FAB in the FabArray to point to elem.
    //
//private:
    void setFab (int  K,
		 FAB* elem);
public:

    /*@ManDoc: Removes the Kth FAB from the FabArray and returns it to the
               caller.  It leaves the Kth element undefined.  It is the
               callers responsibility to delete the pointer when no
               longer needed.
    */
    FAB* remove (int K);

    /*@ManDoc: Delete the Kth element in the FabArray only if the array is
	       managed.  It leaves the Kth element in the FabArray undefined
	       for both managed and unmanaged FabArrays.
    */
    void clear (int K);
    void clear ();

    /*@ManDoc: Delete all elements in the FabArray.  Use this only
	       on unmanaged FabArrays.  It leaves all elements
	       in the FabArray undefined.
    */
    void clearUnmanaged ();

    /*@ManDoc: Set all components in the valid region of each FAB in the
               FabArray to val.
    */
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

    /*@ManDoc: Copy on intersection within FabArray.  Data is copied from
               valid regions to intersecting regions of definition.  The
               purpose is to fill in the boundary regions of each FAB in
               the FabArray.
    */
    void FillBoundary ();

    /*@ManDoc: Same as FillBoundary(), but only copies
               num\_comp components starting at start\_comp.
    */
    void FillBoundary (int start_comp,
                       int num_comp);

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
               FabArray starting at dest\_comp.  nghost boundary zones are
               included in the intersection as well.  This assumes that
               the source and destination FabArray have identical valid
               regions.

    */
    void copy (const FabArray<T,FAB>& src,
               int                    src_comp,
               int                    dest_comp,
               int                    num_comp,
               int                    nghost = 0);

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

    const DistributionMapping& DistributionMap () const
    {
        return distributionMap;
    }
    aString FAName ()
    {
        return faName;
    }
    void SetFAName (const aString& faname)
    {
        faName = faname;
    }
protected:

    int         n_grow;
    int         n_comp;
    BoxArray    boxarray;
    PArray<FAB> fabparray;
    DistributionMapping distributionMap;
    aString faName;

private:
    //
    // These are disallowed.
    //
    FabArray (const FabArray<T,FAB>&);
    FabArray<T,FAB>& operator= (const FabArray<T,FAB>&);
    //
    // This is used locally in all define functions.
    //
    void AllocFabs (int nvar);
};

template <class T, class FAB> class ConstFabArrayIterator;

template <class T, class FAB>
class FabArrayIterator
{
  public:
    //
    // 
    //
    FabArrayIterator(FabArray<T,FAB>& fabarray,
                     const int        nghost = 0);
    //
    //
    //
    const FAB& operator() () const;
    //
    //
    //
    FAB& operator() ();
    //
    //
    //
    const Box& validbox () const;
    //
    //
    //
    const Box& fabbox () const;
    //
    //
    //
    FabArrayIterator<T,FAB>& operator++ ();
    //
    //
    //
    bool isValid ();
    //
    //
    //
    int index () const;
    //
    //
    //
    int NGhost() const;

  protected:

    FabArray<T,FAB>& fabArray;
    int              nGhost;
    int              currentIndex;

  private:
    //
    // These are disallowed.
    //
    FabArrayIterator ();
    FabArrayIterator (const FabArrayIterator<T,FAB>&);
    FabArrayIterator<T,FAB>& operator= (const FabArrayIterator<T,FAB>&);
};

template <class T, class FAB>
class DependentFabArrayIterator
{
  public:
    //
    //
    //
    DependentFabArrayIterator(FabArrayIterator<T,FAB>& controllerfaiter,
                              FabArray<T,FAB>&         dependentfabarray);
    //
    //
    //
    DependentFabArrayIterator(FabArrayIterator<T,FAB>& controllerfaiter,
                              const FabArray<T,FAB>&   dependentfabarray);
    //
    //
    //
    const FAB& operator() () const;
    //
    //
    //
    FAB& operator() ();
    //
    //
    //
    const Box& validbox () const;
    //
    //
    //
    const Box& fabbox () const;
    //
    //
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

template <class T, class FAB>
class ConstFabArrayIterator
{
  public:
    //
    //
    //
    ConstFabArrayIterator(const FabArray<T,FAB>& fabarray,
                          const int              nghost = 0);
    //
    //
    //
    const FAB& operator() () const;
    //
    //
    //
    const Box& validbox () const;
    //
    //
    //
    const Box& fabbox () const;
    //
    //
    //
    ConstFabArrayIterator<T,FAB>& operator++ ();
    //
    //
    //
    int index () const;
    //
    //
    //
    bool isValid () const;
    //
    //
    //
    int NGhost () const;

  protected:

    const FabArray<T,FAB>& fabArray;
    int                    nGhost;
    int                    currentIndex;

  private:
    //
    // These are disallowed.
    //
    ConstFabArrayIterator ();
    ConstFabArrayIterator (const ConstFabArrayIterator<T,FAB> &);
    ConstFabArrayIterator<T,FAB>& operator= (const ConstFabArrayIterator<T,FAB> &);
};

template <class T, class FAB>
class ConstDependentFabArrayIterator
{
  public:
    //
    //
    //
    ConstDependentFabArrayIterator(ConstFabArrayIterator<T,FAB>& cfai,
                                   const FabArray<T,FAB>&        fa);
    //
    //
    //
    const FAB& operator() () const;
    //
    //
    //
    const Box& validbox () const;
    //
    //
    //
    const Box& fabbox () const;
    //
    //
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

  protected:

    void SetId (int id) { fabArrayId = id; }

  private:

    int fabArrayId;
};

class FillBoxId
{
  public:

    FillBoxId () : fillBoxId(-1) {}

    FillBoxId (int            newid,
               const Box&     fillbox,
               const BoxList& filledboxes)
        :
        fillBoxId(newid),
        fillBox(fillbox),
        filledBoxes(filledboxes) {}

    FillBoxId (const FillBoxId& rhs)
        :
        fillBoxId(rhs.fillBoxId),
        fillBox(rhs.fillBox),
        filledBoxes(rhs.filledBoxes) {}

    ~FillBoxId() {}

    int Id () const                       { return fillBoxId;   }
    const Box& box () const               { return fillBox;     }
    const BoxArray& FilledBoxes () const  { return filledBoxes; }

  protected:

    void SetId (int id) { fillBoxId = id; }

  private:

    int      fillBoxId;
    Box      fillBox;
    BoxArray filledBoxes;
};

//
// This enum and the FabCopyDescriptor class should really be nested
// in FabArrayCopyDescriptor  (not done for portability reasons).
//

enum FillType { FillLocally, FillRemotely, Unfillable };

template <class FAB>
struct FabCopyDescriptor
{
    FabCopyDescriptor()
        :
        myProc(-1),
        copyFromProc(-1),
        copyFromIndex(-1),
        fabArrayId(-1),
        fillBoxId(-1),
        srcComp(-1),
        destComp(-1),
        nComp(-1),
        fillType(Unfillable),
        cacheDataAllocated(false),
        localFabSource(NULL)
    {}

    ~FabCopyDescriptor()
    {
        if (cacheDataAllocated)
            delete localFabSource;
    }

    int      myProc;
    int      copyFromProc;
    int      copyFromIndex;
    int      fabArrayId;
    int      fillBoxId;
    int      srcComp;
    int      destComp;
    int      nComp;
    Box      subBox;
    FillType fillType;
    bool     cacheDataAllocated;
    FAB*     localFabSource;  // Either a local fab from this fabarray or a
                              // temp local fab with data from a remote copy.
};

//
// This class orchestrates filling a destination fab of size destFabBox
// from fabarray on the local processor (myProc).
//

template <class T, class FAB>
class FabArrayCopyDescriptor
{
  public:

    FabArrayCopyDescriptor (bool cacheremotedata);

    ~FabArrayCopyDescriptor ();

    FabArrayId RegisterFabArray(FabArray<T,FAB> *fabarray);

    FillBoxId AddBox(const FabArrayId& fabarrayid,
                     const Box&        destFabBox,
		     BoxList&          unfilledBoxes);

    FillBoxId AddBox(const FabArrayId& fabarrayid,
                     const Box&        destFabBox,
		     BoxList&          unfilledBoxes,
		     int               srccomp,
                     int               destcomp,
                     int               numcomp);
    //
    // Add a box but only from FabArray[fabarrayindex].
    //
    FillBoxId AddBox (const FabArrayId& fabarrayid,
                      const Box&        destFabBox,
                      BoxList&          unfilledBoxes,
                      int               fabarrayindex,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp);

    //Array<FillBoxId> AddFabArrayBoxes(const FabArrayId &fabarrayid,
		     //int srccomp, int destcomp, int numcomp);

    void CollectData ();

    void FillFab (const FabArrayId& fabarrayid,
                  const FillBoxId&  fillboxid,
                  FAB&              destFab);

    void PrintStats ();

    //const BoxList &UnfilledBoxes()  { return unfilledBoxes;    }
    //const Box &box() const          { return destBox;          }

    bool DataAvailable () const      { return dataAvailable;    }

    //bool CompletelyFilled()         { return completelyFilled; }

  protected:

    Array< FabArray<T,FAB>*>      fabArrays;
    Array<BoxList>                unfilledBoxes;     // for each fabArray
    Array<bool>                   completelyFilled;  // for each fabArray
    List<FabCopyDescriptor<FAB>*> fabCopyDescList;
    List<FabComTag>               fabComTagList;
    int                           nextFabArrayId;
    int                           nextFillBoxId;
    bool                          dataAvailable;
    bool                          cacheRemoteData;
    //
    // stats
    //
    Array<int>                    nRemoteBoxes;  // per fabArray[index]
    Array<long>                   nRemoteBytes;  // per fabArray[index]
    int                           totalRemoteBoxes;
    //
    // Total bytes to communicate to fill destFab.
    //
    long                          totalRemoteBytes;
    //int nSubBoxes;             // how many pieces is destFabBox broken into

  private:
    //
    // These are disallowed.
    //
    FabArrayCopyDescriptor();
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
const Box&
FabArray<T,FAB>::box (int K) const
{
    return boxarray[K];
}

template <class T, class FAB>
const Box
FabArray<T,FAB>::fabbox (int K) const
{
    //
    // Do not use fabparray[K] because it may not be valid in parallel.
    //
    return (::grow(boxarray[K], n_grow));
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
FAB*
FabArray<T,FAB>::remove (int K)
{
    return fabparray.remove(K);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clear (int K)
{
    fabparray.clear(K);
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::clear ()
{
    fabparray.clear();
}

template <class T, class FAB>
void
FabArray<T,FAB>::clearUnmanaged ()
{
  for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
  {
      FAB *fab = fabparray.remove(fai.index());
      delete fab;
  }
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

//
// Non-inlined stuff.
//

template <class T, class FAB>
FabArray<T,FAB>::FabArray ()
    :
    n_grow(0),
    n_comp(0),
    fabparray(0,PArrayManage)
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
    assert(boxarray.ready());
    bool isok = true;

    for (ConstFabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        if (fabparray.defined(fai.index()))
        {
            const FAB& f = fai();
            if (f.box() != ::grow(box(fai.index()), n_grow))
            {
                isok = false;
            }
        }
        else
        {
            isok = false;
        }
    }

    ParallelDescriptor::ReduceBoolAnd(isok);
    return isok;
}

template <class T, class FAB>
void
FabArray<T,FAB>::define (const BoxArray& bxs,
                         int             ngrow,
                         FabAlloc        alloc)
{
    assert( ! boxarray.ready());
    n_grow = ngrow;
    n_comp = 0;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
    int nbox = boxarray.length();
    fabparray.resize(nbox);
    if (alloc == Fab_allocate)
    {
	n_comp = 1;
	AllocFabs(n_comp);
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::define (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc)
{
    assert(!boxarray.ready());
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
    int nbox = bxs.length();
    fabparray.resize(nbox);
    if (alloc == Fab_allocate)
    {
	AllocFabs(nvar);
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::AllocFabs (int nvar)
{
  for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
  {
      Box tmp(::grow(fai.validbox(), n_grow));
      FAB *f = new FAB(tmp, nvar);
      if (f == 0)
          BoxLib::OutOfMemory(__FILE__, __LINE__);
      fabparray.set(fai.index(), f);
  }
}

template <class T, class FAB>
void
FabArray<T,FAB>::setFab (int  boxno,
                         FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();
    assert(n_comp == elem->nComp());
    assert(boxarray.ready());
    assert(elem->box() == ::grow(boxarray[boxno],n_grow));
    assert(!fabparray.defined(boxno));
    if (distributionMap.ProcessorMap()[boxno] == ParallelDescriptor::MyProc())
    {
        fabparray.set(boxno,elem);
    }
    else
    {
        ParallelDescriptor::Abort("Error in FabArray<T,FAB>::setFab:  trying to setFab remotely");
    }
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
FabArray<T,FAB>::setBndry(T   val,
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

//
// Self copy for given component only.
//

template <class T, class FAB>
void
FabArray<T,FAB>::copy (const FabArray<T,FAB>& farray)
{
    if (farray.boxarray == boxarray)
    {
        for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
        {
            fai().copy(farray[fai.index()]);
        }
    }
    else
    {
        //
        // the index into these arrays is the same index of the FabArray
        // these are intersections relative to this FabArray.
        //
        Array<BoxArray> boxIntersections(length());
        Array<BoxArray> boxIntersectionsLocal(length());
        Array<BoxArray> boxIntersectionsRemote(length());
        Array<BoxArray> boxNonIntersections(length());
        int myproc = ParallelDescriptor::MyProc();
        FabComTag fabComTag;
        int tagSize = sizeof(FabComTag);
        ParallelDescriptor::SetMessageHeaderSize(tagSize);

        int i, ii, iblocal;
        for (i = 0; i < length(); ++i)
        {
            //
            // loop over fabparray.
            //
            boxIntersections[i]    = intersect(farray.boxarray, boxarray[i]);
            boxNonIntersections[i] = complementIn(boxarray[i], farray.boxarray);

            //boxIntersectionsLocal[i].reserve(boxIntersections[i].length());
            if (distributionMap[i] == myproc)
            {
                //
                // fab[i] is local
                //
                iblocal = 0;
                for (ii = 0; ii < farray.boxarray.length(); ++ii)
                {
                    if (farray.distributionMap[ii] == myproc)
                    {
                        //
                        // both fabs are local
                        //
                        Box intersectBox(farray.boxarray[ii] & boxarray[i]);

                        if (intersectBox.ok())
                        {
                            boxIntersectionsLocal[i].resize(iblocal + 1);
                            boxIntersectionsLocal[i].set(iblocal, intersectBox);
                            ++iblocal;
                            fabparray[i].copy(farray[ii]);  // do the local copy
                        }
                    }
                }
            }
            else
            {
                //
                // fab[i] is remote
                //
                iblocal = 0;

                for (ii = 0; ii < farray.boxarray.length(); ++ii)
                {
                    if (farray.distributionMap[ii] == myproc)
                    {
                        //
                        // farray[ii] is local
                        //
                        Box intersectBox(farray.boxarray[ii] & boxarray[i]);

                        if (intersectBox.ok())
                        {
                            boxIntersectionsRemote[i].resize(iblocal + 1);
                            boxIntersectionsRemote[i].set(iblocal, intersectBox);

                            // send the data
                            // send the intersection of farray[ii] to fab[i]s processor

                            FAB tempFab(boxIntersectionsRemote[i][iblocal], farray.nComp());
                            int tempDestComp = 0;
                            int tempSrcComp  = 0;
                            Box srcBox(boxIntersectionsRemote[i][iblocal]);
                            Box destBox(boxIntersectionsRemote[i][iblocal]);
                            tempFab.copy(farray[ii], srcBox, tempSrcComp,
                                         destBox, tempDestComp, farray.nComp());

                            fabComTag.fromProc = myproc;
                            fabComTag.toProc   = distributionMap[i];
                            fabComTag.fabIndex = i;
                            fabComTag.destComp = tempDestComp;
                            fabComTag.nComp    = farray.nComp();
                            fabComTag.box      = destBox;

                            ParallelDescriptor::SendData(fabComTag.toProc, fabComTag,
                                                         tempFab.dataPtr(),
                                                         fabComTag.box.numPts() *
                                                         fabComTag.nComp * sizeof(T));

                            ++iblocal;
                        }  // end if(intersectBox.ok())

                    }  // end if(farray[ii] is local)
                }  // end for(ii...)
            }  // end if(fab[i] is local)
        }  // end for(i...)

        ParallelDescriptor::Synchronize();

        // now receive data if any was sent

        int dataWaitingSize;

        while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
        {
            //
            // Data was sent to this processor.
            //
            if (myproc != fabComTag.toProc)
            {
                cerr << "Error:  _in FabArray::copy(FabArray):  myproc!=fabComTag.toProc : "
                     << myproc << " != " << fabComTag.toProc << endl;
                ParallelDescriptor::Abort("Bad fabComTag.toProc");
            }

            long tmp_long = fabComTag.box.numPts()*fabComTag.nComp*sizeof(T);

            assert(tmp_long < INT_MAX);
            
            int shouldReceiveBytes = int(tmp_long);

            if (dataWaitingSize != shouldReceiveBytes)
            {
                cerr << "Error:  _in FabArray::copy(FabArray):  "
                     << "dataWaitingSize != shouldReceiveBytes:  = "
                     << dataWaitingSize << " != " << shouldReceiveBytes << endl;
                ParallelDescriptor::Abort("Bad receive nbytes");
            }
            if (!fabComTag.box.ok())
            {
                cerr << "Error:  _in FabArray::copy(FabArray):  bad fabComTag.box" << endl;
                cerr << "fabComTag.box = " << fabComTag.box << endl;
                ParallelDescriptor::Abort("Bad receive nbytes");
            }

            FAB tempFab(fabComTag.box, fabComTag.nComp);
            ParallelDescriptor::ReceiveData(tempFab.dataPtr(),
                                            fabComTag.box.numPts() * fabComTag.nComp * sizeof(T));
            int srcComp = 0;
            fabparray[fabComTag.fabIndex].copy(tempFab, fabComTag.box, srcComp,
                                               fabComTag.box, fabComTag.destComp, fabComTag.nComp);
        }

        ParallelDescriptor::Synchronize();  // do we need this here?

/*
    if(myproc == 0) {
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersections[" << i << "] = "
             << boxIntersections[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersectionsLocal[" << i << "] = "
             << boxIntersectionsLocal[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxIntersectionsRemote[" << i << "] = "
             << boxIntersectionsRemote[i] << endl;
      }
      for(i = 0; i < length(); i++) {  // loop over fabparray
        cout << "boxNonIntersections[" << i << "] = "
             << boxNonIntersections[i] << endl;
      }
    }
*/

  }  // end if(same boxarray...)
}

template <class T, class FAB>
void
FabArray<T,FAB>::copy (const FabArray<T,FAB>& src,
                       int                    src_comp,
                       int                    dest_comp,
                       int                    num_comp,
                       int                    nghost)
{
    assert(src.boxarray == boxarray);
    assert(nghost <= n_grow);
    assert(nghost <= src.n_grow);
    int myproc = ParallelDescriptor::MyProc();
    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        Box b(grow(fai.validbox(), nghost));
        fai().copy(src[fai.index()], b, src_comp, b, dest_comp, num_comp);
    }
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
                       int  src_comp,
                       int  dest_comp,
                       int  num_comp) const
{
    copy(dest, dest.box(), src_comp, dest_comp, num_comp);
}

template <class T, class FAB>
void
FabArray<T,FAB>::copy (FAB&       dest,
                       const Box& subbox,
                       int        src_comp,
                       int        dest_comp,
                       int        num_comp) const
{
    assert(dest_comp + num_comp <= dest.nComp());

    int i, overlapLength = 0, overlapIndex = 0;
    int myproc = ParallelDescriptor::MyProc();
    for (i = 0; i < length(); i++)
    {
        if (subbox.intersects(boxarray[i]))
        {
            overlapLength++;
        }
    }

    BoxArray overlapBoxArray;  // make a boxarray of overlaps
    Array<int> overlapProcessors;
    Array<int> overlapFAIndex;
    PArray<FAB> overlapData(PArrayManage);

    if (overlapLength > 0)
    {
        overlapBoxArray.resize(overlapLength);  // make a boxarray of overlaps
        overlapProcessors.resize(overlapLength);
        overlapFAIndex.resize(overlapLength);
        for (i = 0; i < length(); i++)
        {
            Box overlapBox = subbox & boxarray[i];
            if (overlapBox.ok())
            {
                overlapBoxArray.set(overlapIndex, overlapBox);
                overlapProcessors[overlapIndex] = distributionMap[i];
                overlapFAIndex[overlapIndex] = i;
                overlapIndex++;
            }
        }
        //
        // Cannot use a FabArray because they are distributed.
        //
        overlapData.resize(overlapBoxArray.length());
        for (int tfi = 0; tfi < overlapData.length(); tfi++)
        {
            int tfigrow = 0;
            Box tfibox(::grow(overlapBoxArray[tfi], tfigrow));
            FAB *f = new FAB(tfibox, num_comp);
            if (f == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            overlapData.set(tfi, f);
        }

        for (i = 0; i < overlapBoxArray.length(); i++)
        {
            //
            // register for sharing
            //
            ParallelDescriptor::ShareVar(overlapData[i].dataPtr(),
                                         overlapBoxArray[i].numPts() * overlapData[i].nComp() * sizeof(T));
        }

        for (i = 0; i < overlapBoxArray.length(); i++)
        {
            if (myproc == overlapProcessors[i])
            {
                //
                // local data
                //
                int overlapDestComp = 0;
                overlapData[i].copy(fabparray[overlapFAIndex[i]],
                                    src_comp, overlapDestComp, num_comp);
            }
        }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();  // for ShareVar

    if (overlapLength > 0)
    {
        for (i = 0; i < overlapBoxArray.length(); i++)
        {
            bool overlapDataIsLocal = (overlapProcessors[i] == myproc);
            for (int np = 0; np < ParallelDescriptor::NProcs(); np++)
            {
                if (overlapDataIsLocal && myproc != np)
                {
                    //
                    // broadcast data
                    //
                    int doffset = 0;
                    ParallelDescriptor::WriteData(np, overlapData[i].dataPtr(),
                                                  overlapData[i].dataPtr(), doffset,
                                                  overlapData[i].box().numPts() * overlapData[i].nComp() * sizeof(T));
                }
            }
        }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();  // make sure data has arrived

    if (overlapLength > 0)
    {
        for (i = 0; i < overlapBoxArray.length(); i++)
        {
            int overlapSrcComp = 0;
            dest.copy(overlapData[i], overlapBoxArray[i], overlapSrcComp,
                      overlapBoxArray[i], dest_comp, num_comp);
        }

        for (i = overlapBoxArray.length() - 1; i >= 0; i--)
        {
            //
            // unshare in
            //
            ParallelDescriptor::UnshareVar(overlapData[i].dataPtr());  // reverse order
        }
    }  // end if(overlapLength > 0)

    ParallelDescriptor::Synchronize();
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
                         int num_comp,
                         int nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);
    int k;
    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        Box b(::grow(fai.validbox(), nghost));
        for (k = 0; k < num_comp; k++)
        {
            fai().setVal(val, b, comp + k);
	}
    }
}

template <class T, class FAB>
void
FabArray<T,FAB>::setVal (T          val,
                         const Box& region,
                         int        comp,
                         int        num_comp,
                         int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);
    int k;
    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        Box b(::grow(fai.validbox(), nghost));
        b &= region;
        if (b.ok())
        {
            for (k = 0; k < num_comp; k++)
            {
                fai().setVal(val, b, comp + k);
	    }
        }
    }
}

template <class T, class FAB>
inline
void
FabArray<T,FAB>::FillBoundary ()
{
    FillBoundary(0, n_comp);
}

template <class T, class FAB>
void
FabArray<T,FAB>::FillBoundary (int start_comp,
                               int num_comp)
{
/*
  original code:
    
    for (int i = 0; i < length(); ++i) {
        FAB& fab = fs[i];
        for (int j=0; j < ba.nborMax(i); j++) {
            int nbor = ba.nborIndex(i,j);
            // For self copy, don't copy into self.
            if (nbor == i)
                continue;
            Box destbox(ba[nbor]);
            destbox &= fab.box();
            fab.copy(fs[nbor],destbox,start_comp,destbox,start_comp,num_comp);
        }
    }
*/

    FabArrayCopyDescriptor<T,FAB> facd(true);
    FabArrayId faid = facd.RegisterFabArray(this);
    List<FillBoxId> fillBoxIdList;

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        for (int j = 0; j < length(); j++)
        {
            if (j == fai.index())
            {
                //
                // Dont copy into self.
                //
                continue;
            }
            Box destbox(boxarray[j]);
            destbox &= fai().box();
            if (destbox.ok())
            {
                //fai().copy(fs[j],destbox,start_comp,destbox,start_comp,num_comp);
                IndexType boxType(destbox.ixType());
                BoxList unfilledBoxes(boxType);  // unused here
                FillBoxId fbid;
                fbid = facd.AddBox(faid, destbox, unfilledBoxes,
                                   start_comp, start_comp, num_comp);
                fillBoxIdList.append(fbid);
            }
        }
    }

    facd.CollectData();

    ListIterator<FillBoxId> fbidli(fillBoxIdList);

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        for (int j = 0; j < length(); j++)
        {
            if (j == fai.index())
            {
                //
                // Dont copy into self.
                //
                continue;
            }
            Box destbox(boxarray[j]);
            destbox &= fai().box();
            if (destbox.ok())
            {
                assert(fbidli);
                FillBoxId fbid = fbidli();
                ++fbidli;

                assert(destbox == fbid.box());
                FArrayBox overlapFab(fbid.box(), fai().nComp());
                facd.FillFab(faid, fbid, overlapFab);

                fai().copy(overlapFab, overlapFab.box(), start_comp,
                           destbox, start_comp, num_comp);
            }
        }
    }
}

template <class T, class FAB>
FabArrayIterator<T,FAB>::FabArrayIterator(FabArray<T,FAB>& fabarray,
                                          const int        nghost)
    :
    fabArray(fabarray),
    nGhost(nghost),
    currentIndex(0)
{
    //
    // increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc
    //
    int myproc = ParallelDescriptor::MyProc();
    while (currentIndex < fabArray.length() &&
          fabArray.DistributionMap().ProcessorMap()[currentIndex] != myproc)
    {
        ++currentIndex;
    }
}

template <class T, class FAB>
inline
FAB&
FabArrayIterator<T,FAB>::operator() ()
{
  int myproc = ParallelDescriptor::MyProc();
  assert(fabArray.DistributionMap().ProcessorMap()[currentIndex] == myproc);
  return fabArray[currentIndex];
}

template <class T, class FAB>
inline
const FAB&
FabArrayIterator<T,FAB>::operator() () const
{
  int myproc = ParallelDescriptor::MyProc();
  assert(fabArray.DistributionMap().ProcessorMap()[currentIndex] == myproc);
  return fabArray[currentIndex];
}

template <class T, class FAB>
FabArrayIterator<T,FAB>&
FabArrayIterator<T,FAB>::operator++ ()
{
    int myproc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (currentIndex < fabArray.length() &&
           fabArray.DistributionMap().ProcessorMap()[currentIndex] != myproc);

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
  bool valid = (currentIndex >= 0 && currentIndex < fabArray.length());
  if (!valid)
  {
      ParallelDescriptor::Synchronize();
  }
  return valid;
}

template <class T, class FAB>
inline
DependentFabArrayIterator<T,FAB>::DependentFabArrayIterator(
    FabArrayIterator<T,FAB>& controllerfaiter,
    FabArray<T,FAB>&         dependentfabarray)
    :
    fabArray(dependentfabarray),
    controller(&controllerfaiter)
{}

template <class T, class FAB>
inline
DependentFabArrayIterator<T,FAB>::DependentFabArrayIterator(
    FabArrayIterator<T,FAB>& controllerfaiter,
    const FabArray<T,FAB>&   dependentfabarray)
    :
    fabArray(const_cast<FabArray<T,FAB>&>(dependentfabarray)),
    controller(&controllerfaiter)
{}

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
ConstFabArrayIterator<T,FAB>::ConstFabArrayIterator(
    const FabArray<T,FAB> &fabarray,
    const int nghost)
    :
    fabArray(fabarray),
    nGhost(nghost),
    currentIndex(0)
{
    //
    // increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc
    //
    int myproc = ParallelDescriptor::MyProc();
    while (currentIndex < fabArray.length() &&
          fabArray.DistributionMap()[currentIndex] != myproc)
    {
        ++currentIndex;
    }
}

template <class T, class FAB>
inline
const FAB&
ConstFabArrayIterator<T,FAB>::operator() () const
{
    return fabArray[currentIndex];
}

template <class T, class FAB>
ConstFabArrayIterator<T,FAB>&
ConstFabArrayIterator<T,FAB>::operator++ ()
{
    int myproc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (currentIndex < fabArray.length() &&
           fabArray.DistributionMap()[currentIndex] != myproc);

    return *this;
}

template <class T, class FAB>
inline
int
ConstFabArrayIterator<T,FAB>::index () const
{
    return (currentIndex);
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
    bool valid = (currentIndex >= 0 && currentIndex < fabArray.length());
    if (!valid)
    {
        ParallelDescriptor::Synchronize();
    }
    return valid;
}

template <class T, class FAB>
inline
ConstDependentFabArrayIterator<T,FAB>::ConstDependentFabArrayIterator(
    ConstFabArrayIterator<T,FAB>& controllerfaiter,
    const FabArray<T,FAB>&        dependentfabarray)
    :
    fabArray(dependentfabarray),
    controller(&controllerfaiter)
{}

template <class T, class FAB>
inline
int
ConstDependentFabArrayIterator<T,FAB>::index () const
{
    return (controller->index());
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
FabArrayCopyDescriptor<T,FAB>::FabArrayCopyDescriptor (bool cacheremotedata)
    :
    nextFabArrayId(0),
    nextFillBoxId(0),
    dataAvailable(false),
    cacheRemoteData(cacheremotedata),
    completelyFilled(false),
    totalRemoteBoxes(0),
    //nSubBoxes(0),
    totalRemoteBytes(0)
{
    if (cacheRemoteData == false)
    {
        cerr << "Error in FabArrayCopyDescriptor(bool):  "
             << "uncached remote data not implemented." << endl;
        ParallelDescriptor::Abort("Exiting.");
    }

    //nRemoteBoxes.resize(fabArray.length(), 0);
    //nRemoteBytes.resize(fabArray.length(), 0);

    FabComTag fabComTag;
    int tagSize = sizeof(FabComTag);
    ParallelDescriptor::SetMessageHeaderSize(tagSize);

    fabArrays.reserve(16);
    fabComTagList.clear();
}

template <class T, class FAB>
FabArrayId
FabArrayCopyDescriptor<T,FAB>::RegisterFabArray(FabArray<T,FAB>* fabarray)
{
    assert(nextFabArrayId == fabArrays.length());
    fabArrays.resize(nextFabArrayId + 1);
    fabArrays[nextFabArrayId] = fabarray;
    return FabArrayId(nextFabArrayId++);
}

template <class T, class FAB>
inline
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox(const FabArrayId& fabarrayid,
                                      const Box&        destFabBox,
                                      BoxList&          returnedUnfilledBoxes)
{
    return AddBox(fabarrayid, destFabBox, returnedUnfilledBoxes,
                  0, 0, fabArrays[fabarrayid.Id()]->nComp());
}

template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox (const FabArrayId& fabarrayid,
                                       const Box&        destFabBox,
                                       BoxList&          returnedUnfilledBoxes,
                                       int               srccomp,
                                       int               destcomp,
                                       int               numcomp)
{
    int myproc = ParallelDescriptor::MyProc();
    FabComTag fabComTag;

    //nRemoteBoxes.resize(fabArray.length(), 0);
    //nRemoteBytes.resize(fabArray.length(), 0);

    BoxDomain unfilledBoxDomain(destFabBox.ixType());
    BoxList filledBoxes(destFabBox.ixType());
    unfilledBoxDomain.add(destFabBox);
    FabArray<T,FAB> *fabArray = fabArrays[fabarrayid.Id()];
    for (int i = 0; i < fabArray->length(); ++i)
    {
        const Box intersectBox = destFabBox & fabArray->box(i);
        if (intersectBox.isValid())
        {
            filledBoxes.add(intersectBox);
            int remoteProc = fabArray->DistributionMap().ProcessorMap()[i];
            FabCopyDescriptor<FAB> *fcd = new FabCopyDescriptor<FAB>;
            fcd->fabArrayId = fabarrayid.Id();
            fcd->fillBoxId  = nextFillBoxId;
            fcd->subBox = intersectBox;
            fcd->myProc = myproc;
            fcd->copyFromProc = remoteProc;
            fcd->copyFromIndex = i;
            fcd->srcComp  = srccomp;
            fcd->destComp = destcomp;
            fcd->nComp    = numcomp;
            if (myproc == remoteProc)
            {
                //
                // Data is local.
                //
                fcd->fillType = FillLocally;
                fcd->localFabSource = &((*fabArray)[i]);
            }
            else
            {
                //
                // Data is remote.
                //
                dataAvailable = false;
                fcd->fillType = FillRemotely;
                if (cacheRemoteData)
                {
                    fcd->localFabSource = new FAB(intersectBox, numcomp);
                    //fcd->localFabSource->setVal(1.e28);
                    fcd->cacheDataAllocated = true;
                    // send request to send data
                    fabComTag.fabArrayId = fabarrayid.Id();
                    fabComTag.fillBoxId  = nextFillBoxId;
                    fabComTag.fabIndex   = i;
                    fabComTag.procThatNeedsData = myproc;
                    fabComTag.procThatHasData   = remoteProc;
                    fabComTag.subBox     = intersectBox;
                    fabComTag.srcComp    = srccomp;
                    fabComTag.destComp   = destcomp;
                    fabComTag.nComp      = numcomp;
                    //
                    // Do not SendData yet.
                    //
                    fabComTagList.append(fabComTag);
                }
                //nRemoteBoxes[i]++;
                //nRemoteBytes[i]  += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
                totalRemoteBoxes++;
                totalRemoteBytes += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
            }

            fabCopyDescList.append(fcd);
            unfilledBoxDomain.rmBox(intersectBox);
            //nSubBoxes++;
        }
    }
    returnedUnfilledBoxes.clear();
    returnedUnfilledBoxes = unfilledBoxDomain.boxList();
    //if(unfilledBoxes.length() == 0) {
    //completelyFilled = true;
    //}

    return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
}

template <class T, class FAB>
FillBoxId
FabArrayCopyDescriptor<T,FAB>::AddBox(const FabArrayId& fabarrayid,
                                      const Box&        destFabBox,
                                      BoxList&          returnedUnfilledBoxes,
                                      int               fabarrayindex,
                                      int               srccomp,
                                      int               destcomp,
                                      int               numcomp)
{
    int myproc = ParallelDescriptor::MyProc();
    FabComTag fabComTag;

    BoxDomain unfilledBoxDomain(destFabBox.ixType());
    BoxList filledBoxes(destFabBox.ixType());
    unfilledBoxDomain.add(destFabBox);
    FabArray<T,FAB> *fabArray = fabArrays[fabarrayid.Id()];
    assert(fabarrayindex >= 0 && fabarrayindex < fabArray->length());
    //for(int i = 0; i < fabArray->length(); ++i) {
    int i = fabarrayindex;
    const Box &intersectBox = destFabBox & fabArray->box(i);
    if (intersectBox.isValid())
    {
        filledBoxes.add(intersectBox);
        int remoteProc = fabArray->DistributionMap().ProcessorMap()[i];
        FabCopyDescriptor<FAB> *fcd = new FabCopyDescriptor<FAB>;
        fcd->fabArrayId = fabarrayid.Id();
        fcd->fillBoxId  = nextFillBoxId;
        fcd->subBox = intersectBox;
        fcd->myProc = myproc;
        fcd->copyFromProc = remoteProc;
        fcd->copyFromIndex = i;
        fcd->srcComp  = srccomp;
        fcd->destComp = destcomp;
        fcd->nComp    = numcomp;
        if (myproc == remoteProc)
        {
            //
            // Data is local.
            //
            fcd->fillType = FillLocally;
            fcd->localFabSource = &((*fabArray)[i]);
        }
        else
        {
            //
            // Data is remote.
            //
            dataAvailable = false;
            fcd->fillType = FillRemotely;
            if (cacheRemoteData)
            {
                fcd->localFabSource = new FAB(intersectBox, numcomp);
                //fcd->localFabSource->setVal(1.e28);
                fcd->cacheDataAllocated = true;
                // send request to send data
                fabComTag.fabArrayId = fabarrayid.Id();
                fabComTag.fillBoxId  = nextFillBoxId;
                fabComTag.fabIndex   = i;
                fabComTag.procThatNeedsData = myproc;
                fabComTag.procThatHasData   = remoteProc;
                fabComTag.subBox     = intersectBox;
                fabComTag.srcComp    = srccomp;
                fabComTag.destComp   = destcomp;
                fabComTag.nComp      = numcomp;
                //
                // Do not SendData yet.
                //
                fabComTagList.append(fabComTag);
            }
            //nRemoteBoxes[i]++;
            //nRemoteBytes[i]  += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
            totalRemoteBoxes++;
            totalRemoteBytes += intersectBox.numPts() * fabArray->nComp() * sizeof(T);
        }

        fabCopyDescList.append(fcd);
        unfilledBoxDomain.rmBox(intersectBox);
        //nSubBoxes++;
    }
    returnedUnfilledBoxes.clear();
    returnedUnfilledBoxes = unfilledBoxDomain.boxList();
    //if(unfilledBoxes.length() == 0) {
    //completelyFilled = true;
    //}

    return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
}

template <class T, class FAB>
FabArrayCopyDescriptor<T,FAB>::~FabArrayCopyDescriptor()
{
   for (ListIterator<FabCopyDescriptor<FAB>*> fli(fabCopyDescList); fli; ++fli)
   { 
       delete fli();
   }
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::CollectData ()
{
    int dataWaitingSize;
    int myproc = ParallelDescriptor::MyProc();
    FabComTag fabComTag;

    int tagSize = sizeof(FabComTag);
    ParallelDescriptor::SetMessageHeaderSize(tagSize);
    int *nullptr = NULL;
    //
    // Go through the fabComTagList and send all fab requests.
    //
    if (ParallelDescriptor::NProcs() == 1)
    {
        assert(fabComTagList.length() == 0);
    }
    for (ListIterator<FabComTag> fctli(fabComTagList); fctli; ++fctli)
    {
        assert(fctli().procThatNeedsData == myproc);
        assert(fctli().procThatHasData   != myproc);
        assert(fctli().subBox.ok());
        ParallelDescriptor::SendData(fctli().procThatHasData, fctli(), NULL, 0);
    }
    fabComTagList.clear();

    ParallelDescriptor::Synchronize();
    //
    // Check for consistency and correct # of received and expected messages.
    //
    if (!dataAvailable)
    {
        while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
        {
            //
            // Data was sent to this processor.
            //
            assert(dataWaitingSize == 0);
            assert(fabComTag.procThatHasData == myproc);
            assert(fabComTag.subBox.ok());
            assert(fabArrays[fabComTag.fabArrayId]->
                   box(fabComTag.fabIndex).contains(fabComTag.subBox));

            ParallelDescriptor::ReceiveData(nullptr, 0);  // Advance msg hdr.
            //
            // Now send the requested data from the fabs on this processor.
            //
            FAB tempFab(fabComTag.subBox, fabComTag.nComp);
            //tempFab.setVal(1.e29);
            tempFab.copy((*fabArrays[fabComTag.fabArrayId])[fabComTag.fabIndex],
                         fabComTag.subBox, fabComTag.srcComp,
                         fabComTag.subBox, 0, fabComTag.nComp);
            //
            // ^ copy to zero component of temp fab
            //
            ParallelDescriptor::SendData(fabComTag.procThatNeedsData, fabComTag,
                                         tempFab.dataPtr(),
                                         tempFab.box().numPts() *
                                         tempFab.nComp() * sizeof(T));

        }
    }

    ParallelDescriptor::Synchronize();
    //
    // Now collect all remote data into local fab data caches.
    //
    while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
    {
        //
        // Data was sent to this processor.
        //
        assert(fabComTag.procThatHasData != myproc);
        assert(fabComTag.subBox.ok());
        //
        // Check consistency and correct # of received and expected messages.
        //
        if (myproc != fabComTag.procThatNeedsData)
        {
            cerr << "myproc = " << myproc << endl;
            cerr << "fabComTag.procThatNeedsData = " << fabComTag.procThatNeedsData << endl;
            ParallelDescriptor::Abort("FillFab:  bad procThatNeedsData");
        }
        //
        // Find the box in the list and move the data.
        //
        bool matchFound = false;
        for (ListIterator<FabCopyDescriptor<FAB> *> fli(fabCopyDescList);
            fli && ! matchFound;
            ++fli)
        {
            if (fli()->subBox == fabComTag.subBox &&
                fabComTag.fabArrayId == fli()->fabArrayId &&
                fabComTag.fillBoxId == fli()->fillBoxId)
            {
                //
                // copy the data
                //
                ParallelDescriptor::ReceiveData(fli()->localFabSource->dataPtr(),
                                                fabComTag.subBox.numPts() * fabComTag.nComp * sizeof(T));
                matchFound = true;
            }
        }
        if (matchFound == false)
        {
            ParallelDescriptor::Abort("FillFab:  match not found");
        }
    }

    dataAvailable = true;

    ParallelDescriptor::Synchronize();

}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::FillFab (const FabArrayId& fabarrayid,
                                        const FillBoxId&  fillboxid,
                                        FAB&              destFab)
{
    //
    // At this point, all remote data should be in local caches.
    //
    if (!dataAvailable)
    {
        cerr << "Error in FabArrayCopyDescriptor<T,FAB>::FillFab:  ";
        cerr << "data not available" << endl;
        ParallelDescriptor::Abort("FillFab:  data not available");
    }
    //
    // TODO -- optimize this.
    //
    for (ListIterator<FabCopyDescriptor<FAB>*> fli(fabCopyDescList); fli; ++fli)
    {
        if (fli()->fabArrayId == fabarrayid.Id() && fli()->fillBoxId == fillboxid.Id())
        {
            int localSrcComp;
            if (fli()->fillType == FillLocally)
            {
                localSrcComp = fli()->srcComp;
            }
            else
            {
                //
                // FillRemotely
                //
                localSrcComp = 0; // Copy from zero component of local temp fab
            }
            destFab.copy(*(fli()->localFabSource), fli()->subBox, localSrcComp,
                         fli()->subBox, fli()->destComp, fli()->nComp);
        }
    }
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::PrintStats ()
{
  int myproc = ParallelDescriptor::MyProc();

  cout << "----- " << myproc << ":  "
       << "Parallel stats for FabArrayCopyDescriptor:" << '\n';
  cout << "----- " << myproc << ":  "
       << "    totalRemoteBoxes = " << totalRemoteBoxes << '\n';
  cout << "----- " << myproc << ":  "
       << "    totalRemoteBytes = " << totalRemoteBytes << '\n';
  for (int fa = 0; fa < fabArrays.length(); ++fa)
  {
      cout << "fabArrays[" << fa << "]->boxArray() = " << fabArrays[fa]->boxArray()
           << '\n';
  }

  cout << flush();

  /*
  int i;
  cout << "           nSubBoxes = " << nSubBoxes << endl;
  for(i = 0; i < nRemoteBoxes.length(); ++i) {
    cout << "    nRemoteBoxes[" << i << "] = " << nRemoteBoxes[i] << endl;
  }
  for(i = 0; i < nRemoteBytes.length(); ++i) {
    cout << "    nRemoteBytes[" << i << "] = " << nRemoteBytes[i] << endl;
  }
  */
}

#endif /*BL_FABARRAY_H*/
