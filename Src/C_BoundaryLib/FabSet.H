
#ifndef _FABSET_H_
#define _FABSET_H_

//
// $Id: FabSet.H,v 1.18 2000-10-02 20:49:03 lijewski Exp $
//

#include <FArrayBox.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>

//@Man:
/*@Memo:
        A FabSet is a group of FArrayBox's.  The grouping is designed
        specifically to represent regions along the boundary of Box's,
        and are used to implement boundary conditions to discretized
        partial differential equations.
*/        
/*@Doc:
        A FabSet is an array of pointers to FABs.  The standard FAB operators,
        however, have been modified to be more useful for maintaining
        boundary conditions for partial differential equations discretized
        on boxes.
        Under normal circumstances, a FAB will be created for each face of a
        box.  For a group of boxes, a FabSet will be the group of FABs at a
        particular orientation (ie. the lo-i side of each grid in a list).

        Since a FabSet FAB will likely be used to bound a grid box,
        FArrayBox::resize() operations are disallowed.  Also, to preserve
        flexibility in applicable boundary scenarios, intersecting
        FABs in the FabSet are not guaranteed to contain identical data--thus
        copy operations from a FabSet to any FAB-like structure may be
        order-dependent.

        FabSets are used primarily as a data storage mechanism, and are
        manipulated by more sophisticated control classes.
*/

class FabSet
    :
    private MultiFab
{

public:
    //
    //@ManDoc: The default constructor -- you must later call define().
    //
    FabSet ();
    //
    //@ManDoc: Construct a FabSet of specified number of components on 
    //         the grids.
    //
    FabSet (const BoxArray& grids, int ncomp);
    //
    //@ManDoc: The destructor.
    //
    virtual ~FabSet ();
    //
    //@ManDoc: Define a FabSet constructed via default constructor.
    //
    void define (const BoxArray& grids, int ncomp);
    //
    //@ManDoc: Returns number of variables associated with each point.
    //
    MultiFab::nComp;
    //
    //@ManDoc: Returns a constant reference to the BoxArray.
    //
    MultiFab::boxArray;
    //
    //@ManDoc: Returns a constant reference to the Kth Box in the BoxArray.
    //
    MultiFab::box;
    MultiFab::fabbox;
    //
    //@ManDoc: Promote all setVal() member functions back to public.
    //
    MultiFab::setVal;
    //
    //@ManDoc: Copy from this FabSet to the destination FAB.
    //
    const FabSet& copyTo (FArrayBox& dest) const;
    //
    //@ManDoc: Copy specified components of this FabSet to the components of the destination FAB
    //
    const FabSet& copyTo (FArrayBox& dest,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    //@ManDoc: Copy specified components of this FabSet to the components of the destination FAB over the subregion
    //
    const FabSet& copyTo (FArrayBox& dest,
                          const Box& subbox,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    //@ManDoc: Copy from the FAB to this FabSet.
    //
    FabSet& copyFrom (const FArrayBox& src);
    //
    //@ManDoc: Copy from the specified components of a source FAB to the components to this destination FabSet
    //
    FabSet& copyFrom (const FArrayBox& src,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    //@ManDoc: Copy from the specified components of a source FAB to the destination components of this FabSet over the subregion
    //
    FabSet& copyFrom (const FArrayBox& src,
                      const Box&       subbox,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    //@ManDoc: Copy from the source MultiFab to this destination FabSet.
    //
    FabSet& copyFrom (const MultiFab& src,
                      int             nghost,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp);
    //
    //@ManDoc: Copy from the FabSet to this FabSet.
    //
    FabSet& copyFrom (const FabSet& src);
    //
    //@ManDoc: Copy from the specified components of a source FabSet to the components to this destination FabSet
    //
    FabSet& copyFrom (const FabSet& src,
                      int           src_comp,
                      int           dest_comp,
                      int           num_comp);
    //
    //@ManDoc: Copy to the MultiFab from this FabSet.
    //
    void copyTo (MultiFab& dest) const;
    //
    //@ManDoc: Componentwise multiply each datum in this FabSet for specified components.
    //
    FabSet& mult (Real v,
                  int  comp,
                  int  num_comp);
    //
    //@ManDoc: Componentwise multiply each datum in this FabSet for specified components in the subregion.
    //
    FabSet& mult (Real       v,
                  const Box& subreg,
                  int        comp,
                  int        num_comp);
    //
    //@ManDoc: Componentwise add scalar to each datum in this FabSet for specified components.
    //
    FabSet& plus (Real v,
                  int  comp,
                  int  num_comp);
    //
    //@ManDoc: Componentwise add scalar to each datum in this FabSet for specified components in the subregion.
    //
    FabSet& plus (Real       v,
                  const Box& subreg,
                  int        comp,
                  int        num_comp);
    //
    //@ManDoc: Add each datum in MultiFab to each in this FabSet for specified components in the subregion, including ghost cells specified.
    //
    FabSet& plusFrom (const MultiFab& src,
                      int             nghost,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp);
    //
    //@ManDoc: Linear combination: this := a*this + b*src (FabSets must be commensurate).
    //
    FabSet& linComb (Real          a,
                     Real          b,
                     const FabSet& src,
                     int           src_comp,
                     int           dest_comp,
                     int           num_comp);
    //
    //@ManDoc: Linear combination: this := a*mfa + b*mfb on intersection with valid region.
    //
    FabSet& linComb (Real            a,
                     const MultiFab& mfa,
                     int             a_comp,
                     Real            b,
                     const MultiFab& mfb,
                     int             b_comp,
                     int             dest_comp,
                     int             num_comp,
                     int             n_ghost=0);

    friend class FabSetIterator;
    friend class DependentFabSetIterator;
    friend class ConstFabSetIterator;
    friend class ConstDependentFabSetIterator;
    friend class FabSetCopyDescriptor;
    //
    //@ManDoc: Promote operator[] back to public.
    //
    MultiFab::operator[];
    //
    //@ManDoc: Promote clear() back to public.
    //
    MultiFab::clear;
    //
    //@ManDoc: Promote DistributionMap() back to public.
    //
    MultiFab::DistributionMap;
    //
    //@ManDoc: True if there is at least one FAB in the underlying MultiFab.
    //
    bool ready () const { return fabparray.ready(); }
    //
    //@ManDoc: True if there is a FAB at specified index.
    //
    bool defined (int i) const { return fabparray.defined(i); }
    //
    //@ManDoc: Flush cached stuff between regriddings.
    //
    static void FlushCache ();

protected:
    //
    // copyFrom(MultiFab) and plusFrom(MultiFab) use common function.
    //
    enum How { COPYFROM, PLUSFROM };

    void DoIt (const MultiFab& src,
               int             ngrow,
               int             scomp,
               int             dcomp,
               int             ncomp,
               How             how);
};

//
//@Man:
//@Memo: An iterator over the FABs in a FabSet.
/*@Doc:

  This class encapsulates an iterator over the FABs in a FabSet such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the FabSet that are actually owned by that CPU.
*/

class FabSetIterator
    :
    public MultiFabIterator
{
public:
    //
    //@ManDoc: Construct a FabSetIterator.
    //
    FabSetIterator (FabSet& fabset)
        :
        MultiFabIterator(fabset) {}
    //
    //@ManDoc: The destructor.
    //
    ~FabSetIterator () {}
};

//
//@Man:
//@Memo: An iterator over a dependent FabSet.
/*@Doc:

  This class is used to iterator over a FabSet (the dependent FabSet) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a FabSetIterator.  This allows us to iterate
  over a FabSet and access the associated FAB(s) in one or more dependent
  FabSets.
*/

class DependentFabSetIterator
    :
    public DependentMultiFabIterator
{
public:
    //
    //@ManDoc: Constructor over a mutable dependent FabSet.
    //
    DependentFabSetIterator (FabSetIterator& controllerfsiter,
                             FabSet&         fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: Constructor over an immutable dependent FabSet.
    //
    DependentFabSetIterator (FabSetIterator& controllerfsiter,
                             const FabSet&   fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: Constructor over a mutable dependent FabSet initialized with a MultiFabIterator.
    //
    DependentFabSetIterator (MultiFabIterator& controllerfsiter,
                             FabSet&           fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: Constructor over an immutable dependent FabSet initialized with a MultiFabIterator.
    //
    DependentFabSetIterator (MultiFabIterator& controllerfsiter,
                             const FabSet&     fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: The destructor.
    //
    ~DependentFabSetIterator () {}
};

//
//@Man:
//@Memo: A constant iterator over the FABs in a FabSet.
/*@Doc:

  This class encapsulates an iterator over the FABs in a FabSet such that
  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
  only those FABs in the FabSet that are actually owned by that CPU.  It is
  a constant iterator; you can access the FABs but you cannot change them.
*/

class ConstFabSetIterator
    :
    public ConstMultiFabIterator
{
public:
    //
    //@ManDoc: Construct a ConstFabSetIterator.
    //
    ConstFabSetIterator (const FabSet &fabset)
        :
        ConstMultiFabIterator(fabset) {}

    //
    //@ManDoc: The destructor.
    //
    ~ConstFabSetIterator () {}
};

//
//@Man:
//@Memo: A constant iterator over a dependent FabSet.
/*@Doc:

  This class is used to iterator over a FabSet (the dependent FabSet) that
  is defined over a BoxArray that is the same as (or sufficiently similar to)
  one for which we already have a ConstFabSetIterator.  This allows us to
  iterate over a constant FabSet and access the associated FAB(s) in one or
  more dependent FabSets.
*/

class ConstDependentFabSetIterator
    :
    public ConstDependentMultiFabIterator
{
public:
    //
    //@ManDoc: Constructor taking a ConstFabSetIterator.
    //
    ConstDependentFabSetIterator (ConstFabSetIterator& controllerfsiter,
                                  const FabSet&        fabset)
        :
        ConstDependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: Constructor taking a ConstMultiFabIterator.
    //
    ConstDependentFabSetIterator (ConstMultiFabIterator& controllerfsiter,
                                  const FabSet&          fabset)
        :
        ConstDependentMultiFabIterator(controllerfsiter, fabset) {}
    //
    //@ManDoc: The destructor.
    //
    ~ConstDependentFabSetIterator () {}
};

typedef FabArrayId FabSetId;

class FabSetCopyDescriptor
    :
    public MultiFabCopyDescriptor
{
public:
    FabSetCopyDescriptor ()
        :
        MultiFabCopyDescriptor() {}

    ~FabSetCopyDescriptor () {}

    FabSetId RegisterFabSet (FabSet* fabset)
    {
        return RegisterMultiFab(fabset);
    }
private:
    //
    // These are disallowed.
    //
    FabSetCopyDescriptor (const FabSetCopyDescriptor& rhs);
    FabSetCopyDescriptor& operator= (const FabSetCopyDescriptor& rhs);
};

//
// Inlines.
//

inline
void
FabSet::define (const BoxArray& grids, int ncomp)
{
    MultiFab* tmp = this;

    tmp->define(grids, ncomp, 0, Fab_allocate);
}

inline
const FabSet&
FabSet::copyTo (FArrayBox& dest) const
{
    copy(dest);
    return *this;
}

inline
const FabSet&
FabSet::copyTo (FArrayBox& dest,
                int        src_comp,
                int        dest_comp,
                int        num_comp) const
{
    copy(dest,src_comp,dest_comp,num_comp);
    return *this;
}

inline
const FabSet&
FabSet::copyTo (FArrayBox& dest,
                const Box& subbox,
                int        src_comp,
                int        dest_comp,
                int        num_comp) const
{
    copy(dest,subbox,src_comp,dest_comp,num_comp);
    return *this;
}

inline
void
FabSet::copyTo (MultiFab& dest) const
{
    dest.copy(*this);
}

inline
FabSet&
FabSet::copyFrom (const FabSet& src)
{
    copy(src);
    return *this;
}

inline
FabSet&
FabSet::copyFrom (const FabSet& src,
                  int           src_comp,
                  int           dest_comp,
                  int           num_comp)
{
    copy(src,src_comp,dest_comp,num_comp);
    return *this;
}

//
// The following are different from MultiFab only in the return value
//

inline
FabSet&
FabSet::plus (Real v,
              int  comp,
              int  num_comp)
{
    MultiFab* tmp = this;
    tmp->plus(v, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::plus (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    MultiFab* tmp = this;
    tmp->plus(v, subreg, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::mult (Real v,
              int  comp,
              int  num_comp)
{
    MultiFab* tmp = this;
    tmp->mult(v, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::mult (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    MultiFab* tmp = this;
    tmp->mult(v, subreg, comp, num_comp);
    return *this;
}

#endif /*_FABSET_H_*/
