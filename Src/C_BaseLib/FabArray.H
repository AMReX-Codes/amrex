//BL_COPYRIGHT_NOTICE

#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

//
// $Id: FabArray.H,v 1.92 1999-03-18 22:12:01 lijewski Exp $
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

#include <Assert.H>
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

class FillBoxId
{
  public:

    FillBoxId () : fillBoxId(-1), fabIndex(-1) {}

    FillBoxId (int            newid,
               const Box&     fillbox,
               const BoxList& filledboxes,
               int            fabindex = -1)
        :
        filledBoxes(filledboxes),
        fillBox(fillbox),
        fillBoxId(newid),
        fabIndex(fabindex) {}

    FillBoxId (const FillBoxId& rhs)
        :
        filledBoxes(rhs.filledBoxes),
        fillBox(rhs.fillBox),
        fillBoxId(rhs.fillBoxId),
        fabIndex(rhs.fabIndex) {}

    ~FillBoxId() {}

    int Id () const                       { return fillBoxId;    }
    int FabIndex () const                 { return fabIndex;     }
    void FabIndex (int fabindex)          { fabIndex = fabindex; }
    const Box& box () const               { return fillBox;      }
    const BoxArray& FilledBoxes () const  { return filledBoxes;  }

  protected:

    void SetId (int id) { fillBoxId = id; }

  private:

    BoxArray filledBoxes;
    Box      fillBox;
    int      fillBoxId;
    int      fabIndex;
};

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

    /*@ManDoc: Removes the Kth FAB from the FabArray and returns it to the
               caller.  It leaves the Kth element undefined.  It is the
               callers responsibility to delete the pointer when no
               longer needed.  Use this at your own peril; it breaks
               the parallel paradigm.
    */
    FAB* remove (int K);

    /*@ManDoc: Delete the Kth element in the FabArray only if the array is
               managed.  It leaves the Kth element in the FabArray undefined
               for both managed and unmanaged FabArrays.   Use this at your
               own peril; it breaks the parallel paradigm.
    */
    void clear (int K);

    /*@ManDoc: Releases FAB memory in the FabArray.
    */
    void clear ();

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

  protected:

    void SetId (int id) { fabArrayId = id; }

  private:

    int fabArrayId;
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

    ~FabCopyDescriptor()
    {
        if (cacheDataAllocated)
            delete localFabSource;
    }

    FAB*     localFabSource;  // Either a local fab from this fabarray or a
                              // temp local fab with data from a remote copy.
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
};

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

    void CollectData ();

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
    return ::grow(boxarray[K], n_grow);
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
    assert(boxarray.ready());

    long isok = 1;

    ConstFabArrayIterator<T,FAB> fai(*this);

    for ( ; fai.isValid() && isok; ++fai)
    {
        if (fabparray.defined(fai.index()))
        {
            if (fai().box() != ::grow(box(fai.index()),n_grow))
                isok = 0;
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
    assert(!boxarray.ready());
    n_grow = ngrow;
    n_comp = 0;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
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
    assert(!boxarray.ready());
    n_grow = ngrow;
    n_comp = nvar;
    boxarray.define(bxs);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
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
        Box tmp = ::grow(fai.validbox(), n_grow);
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

    assert(n_comp == elem->nComp());
    assert(boxarray.ready());
    assert(elem->box() == ::grow(boxarray[boxno],n_grow));
    assert(!fabparray.defined(boxno));
    assert(distributionMap[boxno] == ParallelDescriptor::MyProc());

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
                       int                    src_comp,
                       int                    dest_comp,
                       int                    num_comp)
{
    if (src.boxarray == boxarray)
    {
        for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
        {
            DependentFabArrayIterator<T,FAB> dfai(fai, src);

            if (fai.validbox().intersects(dfai.validbox()))
            {
                Box isect = fai.validbox() & dfai.validbox();

                fai().copy(dfai(),isect,src_comp,isect,dest_comp,num_comp);
            }
        }
        return;
    }

    const int MyProc = ParallelDescriptor::MyProc();

#ifdef BL_USE_MPI
    FabComTag         tag;
    vector<FabComTag> SndTags;
    Array<int>        Snds(ParallelDescriptor::NProcs(),0);
    Array<int>        Rcvs(ParallelDescriptor::NProcs(),0);
#endif

    for (int i = 0; i < length(); ++i)
    {
        if (distributionMap[i] == MyProc)
        {
            for (int ii = 0; ii < src.boxarray.length(); ++ii)
            {
                if (src.distributionMap[ii] == MyProc)
                {
                    //
                    // Both fabs are local
                    //
                    if (src.boxarray[ii].intersects(boxarray[i]))
                    {
                        Box intersect = src.boxarray[ii] & boxarray[i];

                        fabparray[i].copy(src[ii],
                                          intersect,
                                          src_comp,
                                          intersect,
                                          dest_comp,
                                          num_comp);
                    }
                }
            }
        }
#ifdef BL_USE_MPI
        else
        {
            for (int ii = 0; ii < src.boxarray.length(); ++ii)
            {
                if (src.distributionMap[ii] == MyProc)
                {
                    if (src.boxarray[ii].intersects(boxarray[i]))
                    {
                        tag.fromProc = MyProc;
                        tag.toProc   = distributionMap[i];
                        tag.fabIndex = i;
                        tag.box      = src.boxarray[ii] & boxarray[i];
                        //
                        // Use tag.fineIndex to store index into `src'.
                        //
                        tag.fineIndex = ii;
                        SndTags.push_back(tag);
                        Snds[distributionMap[i]]++;
                    }
                }
            }
        }
#endif /*BL_USE_MPI*/
    }

#ifdef BL_USE_MPI
    if (ParallelDescriptor::NProcs() == 1)
        return;

    assert(Snds[MyProc] == 0);

    int rc;

    static RunStats stats("fabarray_copy");
    static RunStats mpi_recv("mpi_recv");
    static RunStats mpi_send("mpi_send");
    static RunStats mpi_wait("mpi_wait");
    static RunStats mpi_gath("mpi_gather");
    static RunStats mpi_redu("mpi_reduce");

    stats.start();
    //
    // Set Rcvs[i] to # of blocks we expect to get from CPU i ...
    //
    mpi_gath.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if ((rc = MPI_Gather(&Snds[i],
                             1,
                             MPI_INT,
                             Rcvs.dataPtr(),
                             1,
                             MPI_INT,
                             i,
                             MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }
    mpi_gath.end();

    assert(Rcvs[MyProc] == 0);

    int NumRcvs = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        NumRcvs += Rcvs[i];

    Array<MPI_Request> req_cd(ParallelDescriptor::NProcs());
    Array<CommData>    rcv_cd(NumRcvs);
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    assert(sizeof(CommData) == CommData::DIM*sizeof(int));
    //
    // Post one receive for each chunk being sent by other CPUs.
    // This is the CommData describing the FAB data that will be sent.
    //
    int idx = 0;

    mpi_recv.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            if ((rc = MPI_Irecv(&rcv_cd[idx],
                                Rcvs[i] * CommData::DIM,
                                MPI_INT,
                                i,
                                139,
                                MPI_COMM_WORLD,
                                &req_cd[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            idx += Rcvs[i];
        }
    }
    mpi_recv.end();

    assert(idx == NumRcvs);
    //
    // Now send the CommData.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Snds[i] > 0)
        {
            Array<CommData> senddata(Snds[i]);

            int Processed = 0;

            for (int j = 0; j < SndTags.size(); j++)
            {
                if (SndTags[j].toProc == i)
                {
                    CommData data(0,                    // Not used.
                                  SndTags[j].fabIndex,
                                  MyProc,
                                  0,                    // Not used.
                                  0,                    // Not used.
                                  0,                    // Not used.
                                  0,                    // Not used.
                                  SndTags[j].box);

                    senddata[Processed++] = data;
                }
            }

            assert(Processed == Snds[i]);
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(senddata.dataPtr(),
                                senddata.length() * CommData::DIM,
                                MPI_INT,
                                i,
                                139,
                                MPI_COMM_WORLD)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();
        }
    }
    //
    // Post one receive for data being sent by CPU i ...
    //
    MPI_Status status;

    Array<T*>          fab_data(ParallelDescriptor::NProcs());
    Array<MPI_Request> req_data(ParallelDescriptor::NProcs());

    idx = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&req_cd[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();
            //
            // Got to figure out # of Ts to expect from this CPU.
            //
            size_t N = 0;

            for (int j = 0; j < Rcvs[i]; j++)
                N += rcv_cd[idx+j].box().numPts() * num_comp;

            assert(N < INT_MAX);
            assert(!(The_FAB_Arena == 0));

            fab_data[i] = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));

            mpi_recv.start();
            if ((rc = MPI_Irecv(fab_data[i],
                                int(N),
                                mpi_data_type(fab_data[i]),
                                i,
                                537,
                                MPI_COMM_WORLD,
                                &req_data[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_recv.end();

            idx += Rcvs[i];
        }
    }

    FAB fab;

    assert(idx == NumRcvs);
    //
    // Send the agglomerated FAB data.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Snds[i] > 0)
        {
            size_t N = 0;

            for (int j = 0; j < SndTags.size(); j++)
                if (SndTags[j].toProc == i)
                    N += SndTags[j].box.numPts() * num_comp;

            assert(N < INT_MAX);
            assert(!(The_FAB_Arena == 0));

            T* data = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));
            T* dptr = data;

            for (int j = 0; j < SndTags.size(); j++)
            {
                if (SndTags[j].toProc == i)
                {
                    fab.resize(SndTags[j].box, num_comp);

                    assert(src.distributionMap[SndTags[j].fineIndex] == MyProc);

                    fab.copy(src[SndTags[j].fineIndex],
                             SndTags[j].box,
                             src_comp,
                             SndTags[j].box,
                             0,
                             num_comp);

                    int count = SndTags[j].box.numPts() * num_comp;

                    memcpy(dptr, fab.dataPtr(), count*sizeof(T));

                    dptr += count;
                }
            }

            assert(data + N == dptr);
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(data,
                                int(N),
                                mpi_data_type(data),
                                i,
                                537,
                                MPI_COMM_WORLD)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();

            assert(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(data);
        }
    }
    //
    // Now receive and unpack FAB data.
    //
    idx = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&req_data[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();

            int Processed = 0;

            T* dptr = fab_data[i];

            assert(!(dptr == 0));

            for (int j = 0; j < Rcvs[i]; j++)
            {
                const CommData& cd = rcv_cd[idx+j];

                fab.resize(cd.box(), num_comp);

                int N = fab.box().numPts() * fab.nComp();

                memcpy(fab.dataPtr(), dptr, N*sizeof(T));

                assert(distributionMap[cd.fabindex()] == MyProc);

                fabparray[cd.fabindex()].copy(fab,
                                              fab.box(),
                                              0,
                                              fab.box(),
                                              dest_comp,
                                              num_comp);
                dptr += N;

                Processed++;
            }

            assert(Processed == Rcvs[i]);
            assert(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(fab_data[i]);

            fab_data[i] = 0;

            idx += Rcvs[i];
        }
    }

    assert(idx == NumRcvs);

    stats.end();
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

#ifdef BL_USE_MPI
    static RunStats mpi_bcast("mpi_broadcast");

    const int MyProc = ParallelDescriptor::MyProc();

    int overlapLength = 0, overlapIndex = 0;

    for (int i = 0; i < length(); i++)
        if (subbox.intersects(boxarray[i]))
            overlapLength++;

    BoxArray    overlapBA(overlapLength);
    Array<int>  overlapProc(overlapLength);
    Array<int>  overlapFAIdx(overlapLength);
    PArray<FAB> overlapData(overlapLength,PArrayManage);

    if (overlapLength > 0)
    {
        for (int i = 0; i < length(); i++)
        {
            if (subbox.intersects(boxarray[i]))
            {
                overlapBA.set(overlapIndex, subbox & boxarray[i]);
                overlapProc[overlapIndex] = distributionMap[i];
                overlapFAIdx[overlapIndex++]  = i;
            }
        }
    }
    for (int i = 0; i < overlapLength; i++)
    {
        overlapData.set(i, new FAB(overlapBA[i], num_comp));

        if (MyProc == overlapProc[i])
        {
            overlapData[i].copy(fabparray[overlapFAIdx[i]],
                                src_comp,
                                0,
                                num_comp);
        }
    }

    mpi_bcast.start();
    for (int i = 0, rc; i < overlapLength; i++)
    {
        if ((rc = MPI_Bcast(overlapData[i].dataPtr(),
                            overlapData[i].box().numPts()*num_comp,
                            mpi_data_type(overlapData[i].dataPtr()),
                            overlapProc[i],
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }
    mpi_bcast.end();

    for (int i = 0; i < overlapLength; i++)
    {
        dest.copy(overlapData[i],
                  overlapBA[i],
                  0,
                  overlapBA[i],
                  dest_comp,
                  num_comp);
    }
#else
    for (int j = 0; j < length(); ++j)
    {
        //
        // Test all distributed objects.
        //
        if (boxarray[j].intersects(subbox))
        {
            //
            // Restrict copy to domain of validity of source.
            //
            Box destbox = boxarray[j] & subbox;
            dest.copy(fabparray[j],
                      destbox,
                      src_comp,
                      destbox,
                      dest_comp,
                      num_comp);
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
                         int num_comp,
                         int nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        fai().setVal(val, ::grow(fai.validbox(), nghost), comp, num_comp);
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

    for (FabArrayIterator<T,FAB> fai(*this); fai.isValid(); ++fai)
    {
        Box b = ::grow(fai.validbox(),nghost) & region;

        if (b.ok())
            fai().setVal(val, b, comp, num_comp);
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
  assert(fabArray.DistributionMap()[currentIndex] == ParallelDescriptor::MyProc());

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
  assert(fabArray.DistributionMap()[currentIndex] == ParallelDescriptor::MyProc());

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
    assert(currentIndex >= 0);

#if defined(BL_USE_MPI) && defined(HG_DEBUG)
    extern bool HG_is_debugging;
    if (HG_is_debugging)
    {
	if (currentIndex < fabArray.length()) return true;
	int rc = MPI_Barrier(MPI_COMM_WORLD);
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
    assert(controllerfaiter.theFabArray().DistributionMap() ==
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
    assert(controllerfaiter.theFabArray().DistributionMap() ==
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
    assert(currentIndex >= 0);

#if defined(BL_USE_MPI) && defined(HG_DEBUG)
    extern bool HG_is_debugging;
    if (HG_is_debugging)
    {
	if (currentIndex < fabArray.length()) return true;
	int rc = MPI_Barrier(MPI_COMM_WORLD);
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
    assert(controllerfaiter.theFabArray().DistributionMap() ==
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
    assert(fabArrays.size() == fabCopyDescList.size());

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

    assert(faindex >= 0 && faindex < fabArray->length());

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
            fcd->localFabSource = &((*fabArray)[faindex]);
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

    return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
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

    return FillBoxId(nextFillBoxId++, destFabBox, filledBoxes);
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
FabArrayCopyDescriptor<T,FAB>::CollectData ()
{
    dataAvailable = true;

    if (ParallelDescriptor::NProcs() == 1)
        return;

#ifdef BL_USE_MPI
    static RunStats mpi_stat("collect_data");
    static RunStats mpi_recv("mpi_recv");
    static RunStats mpi_send("mpi_send");
    static RunStats mpi_wait("mpi_wait");
    static RunStats mpi_gath("mpi_gather");
    static RunStats mpi_redu("mpi_reduce");

    mpi_stat.start();

    const int MyProc = ParallelDescriptor::MyProc();

    assert(ParallelDescriptor::NProcs() > 1 || fabComTagList.size() == 0);

    Array<int> Rcvs(ParallelDescriptor::NProcs(),0);
    Array<int> Snds(ParallelDescriptor::NProcs(),0);
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (int i = 0; i < fabComTagList.size(); i++)
    {
        assert(fabComTagList[i].box.ok());
        assert(fabComTagList[i].procThatNeedsData == MyProc);
        assert(fabComTagList[i].procThatHasData   != MyProc);

        Rcvs[fabComTagList[i].procThatHasData]++;
    }

    assert(Rcvs[MyProc] == 0);

    int rc;
    //
    // Set Snds[i] to # of blocks we must send to CPU i ...
    //
    mpi_gath.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if ((rc = MPI_Gather(&Rcvs[i],
                             1,
                             MPI_INT,
                             Snds.dataPtr(),
                             1,
                             MPI_INT,
                             i,
                             MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }
    mpi_gath.end();

    assert(Snds[MyProc] == 0);

    int NumReqs = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        NumReqs += Snds[i];

    Array<MPI_Request> reqs_cd(ParallelDescriptor::NProcs());
    Array<CommData>    recv_cd(NumReqs);
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    assert(sizeof(CommData) == CommData::DIM*sizeof(int));
    //
    // Post one receive for each chunk being requested by other CPUs.
    // This is the CommData describing what FAB data needs to be sent.
    //
    int idx = 0;

    mpi_recv.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Snds[i] > 0)
        {
            if ((rc = MPI_Irecv(&recv_cd[idx],
                                Snds[i] * CommData::DIM,
                                MPI_INT,
                                i,
                                471,
                                MPI_COMM_WORLD,
                                &reqs_cd[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            idx += Snds[i];
        }
    }
    mpi_recv.end();

    assert(idx == NumReqs);
    //
    // Make one send to each CPU from which we want data.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            Array<CommData> senddata(Rcvs[i]);

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

            assert(Processed == Rcvs[i]);
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(senddata.dataPtr(),
                                senddata.length() * CommData::DIM,
                                MPI_INT,
                                i,
                                471,
                                MPI_COMM_WORLD)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();
        }
    }
    //
    // Post one receive for data being requested by other CPUs.
    //
    Array<T*>          fab_data(ParallelDescriptor::NProcs());
    Array<MPI_Request> reqs_data(ParallelDescriptor::NProcs());

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
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

            assert(N < INT_MAX);
            assert(!(The_FAB_Arena == 0));

            fab_data[i] = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));

            mpi_recv.start();
            if ((rc = MPI_Irecv(fab_data[i],
                                int(N),
                                mpi_data_type(fab_data[i]),
                                i,
                                917,
                                MPI_COMM_WORLD,
                                &reqs_data[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_recv.end();
        }
    }
    //
    // Send the agglomerated FAB data.
    //
    FAB        fab;
    MPI_Status status;

    idx = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Snds[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&reqs_cd[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();

            size_t N = 0;

            for (int j = 0; j < Snds[i]; j++)
                N += recv_cd[idx+j].box().numPts() * recv_cd[idx+j].nComp();

            assert(N < INT_MAX);
            assert(!(The_FAB_Arena == 0));

            T* data = static_cast<T*>(The_FAB_Arena->alloc(N*sizeof(T)));
            T* dptr = data;

            for (int j = 0; j < Snds[i]; j++)
            {
                int k = idx+j;

                assert(recv_cd[k].id() == 0);
                assert(recv_cd[k].fromproc() == i);

                fab.resize(recv_cd[k].box(), recv_cd[k].nComp());

                int count = fab.box().numPts() * fab.nComp();

                fab.copy((*fabArrays[recv_cd[k].fabarrayid()])[recv_cd[k].fabindex()],
                         fab.box(),
                         recv_cd[k].srcComp(),
                         fab.box(),
                         0,
                         recv_cd[k].nComp());

                memcpy(dptr, fab.dataPtr(), count*sizeof(T));

                dptr += count;
            }

            assert(data + N == dptr);
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(data,
                                int(N),
                                mpi_data_type(data),
                                i,
                                917,
                                MPI_COMM_WORLD)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();

            assert(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(data);

            idx += Snds[i];
        }
    }

    assert(idx == NumReqs);
    //
    // Now receive and unpack FAB data.
    //
    pair<FCDMapIter,FCDMapIter> match;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&reqs_data[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();

            int Processed = 0;

            T* dptr = fab_data[i];

            assert(!(dptr == 0));

            for (int j = 0; j < fabComTagList.size(); j++)
            {
                if (fabComTagList[j].procThatHasData == i)
                {
                    const FabComTag& tag = fabComTagList[j];

                    match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                    FCDMapIter fmi = match.first;

                    for ( ; fmi != match.second; ++fmi)
                    {
                        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                        assert(fcdp->fillBoxId == tag.fillBoxId);

                        if (fcdp->subBox == tag.box)
                        {
                            T* dataPtr = fcdp->localFabSource->dataPtr();
                            assert(!(dataPtr == 0));
                            assert(fcdp->localFabSource->box() == tag.box);
                            int N = tag.box.numPts()*tag.nComp;
                            memcpy(dataPtr, dptr, N*sizeof(T));
                            dptr += N;
                            Processed++;
                            break;
                        }
                    }

                    assert(!(fmi == match.second));
                }
            }

            assert(Processed == Rcvs[i]);
            assert(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(fab_data[i]);

            fab_data[i] = 0;
        }
    }

    mpi_stat.end();
#endif /*BL_USE_MPI*/
}

template <class T, class FAB>
void
FabArrayCopyDescriptor<T,FAB>::FillFab (FabArrayId       faid,
                                        const FillBoxId& fillboxid,
                                        FAB&             destFab)
{
    assert(dataAvailable);

    FCDMapIter fmi = fabCopyDescList[faid.Id()].lower_bound(fillboxid.Id());

    FCDMapIter end = fabCopyDescList[faid.Id()].upper_bound(fillboxid.Id());

    for ( ; fmi != end; ++fmi)
    {
        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

        assert(fcdp->fillBoxId == fillboxid.Id());

        destFab.copy(*(fcdp->localFabSource),
                     fcdp->subBox,
                     (fcdp->fillType == FillLocally) ? fcdp->srcComp : 0,
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
    assert(dataAvailable);

    FCDMapIter fmi = fabCopyDescList[faid.Id()].lower_bound(fillboxid.Id());

    assert(!(fmi == fabCopyDescList[faid.Id()].end()));

    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

    assert(fcdp->fillBoxId == fillboxid.Id());

    assert(fcdp->subBox.sameSize(destBox));

    destFab.copy(*(fcdp->localFabSource),
                 fcdp->subBox,
                 (fcdp->fillType == FillLocally) ? fcdp->srcComp : 0,
                 destBox,
                 fcdp->destComp,
                 fcdp->nComp);

    assert(++fmi == fabCopyDescList[faid.Id()].upper_bound(fillboxid.Id()));
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

#endif /*BL_FABARRAY_H*/
