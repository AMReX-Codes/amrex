//BL_COPYRIGHT_NOTICE

#ifndef _FABSET_H_
#define _FABSET_H_

//
// $Id: FabSet.H,v 1.5 1998-03-30 20:05:35 lijewski Exp $
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
        FArrayBox::resize operations are disallowed.  Also, to preserve
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
    //@ManDoc: default constructor.
    //
    FabSet ();
    //
    //@ManDoc: construct _len sets of FABs.
    //
    FabSet (int _len);
    //
    //@ManDoc: destructor.
    //
    virtual ~FabSet ();
    //
    // promote these back to public.
    //
    MultiFab::setVal;
    //
    //@ManDoc: copy from this FabSet to the destination FAB.
    //
    const FabSet& copyTo (FArrayBox& dest) const;
    //
    //@ManDoc: copy specified components of this FabSet to the components of the destination FAB
    //
    const FabSet& copyTo (FArrayBox& dest,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    //@ManDoc: copy specified components of this FabSet to the compoentns of the destination FAB over the subregion
    //
    const FabSet& copyTo (FArrayBox& dest,
                          const Box& subbox,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    //@ManDoc: copy from the FAB to this FabSet.
    //
    FabSet& copyFrom (const FArrayBox& src);
    //
    //@ManDoc: copy from the specified components of a source FAB to the components to this destination FabSet
    //
    FabSet& copyFrom (const FArrayBox& src,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    //@ManDoc: copy from the specified components of a source FAB to the destination components of this FabSet over the subregion
    //
    FabSet& copyFrom (const FArrayBox& src,
                      const Box&       subbox,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    //@ManDoc: copy from the source MultiFab to this destination FabSet.
    //
    FabSet& copyFrom (const MultiFab& src,
                      int             nghost,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp);
    //
    //@ManDoc: copy from this source FabSet to the destination MultiFab.
    //
    const FabSet& copyTo (MultiFab& dest,
                          int       nghost,
                          int       src_comp,
                          int       dest_comp,
                          int       num_comp) const;
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

    MultiFab::operator[];
    MultiFab::clear;
    MultiFab::DistributionMap;

    const Box& box (int K) const { return fabboxarray[K]; }
    Box fabbox (int K) const { return fabboxarray[K]; }

    void setFab (int        boxno,
                 FArrayBox* fab);

    void resize (int newsize) { fabparray.resize(newsize); }
    void setBox (int index, const Box& b) { fabboxarray.set(index, b); }
    void DefineGrids (const BoxArray& newgrids)
    {
        boxarray.define(newgrids);
        fabboxarray.resize(newgrids.length());
    }
    bool ready ()             { return fabparray.ready();    }
    bool defined (int i)      { return fabparray.defined(i); }

    void DefineDistributionMap (const BoxArray& boxarray)
    {
        distributionMap.define(ParallelDescriptor::NProcs(), boxarray);
    }

  private:
    BoxArray fabboxarray;
};

class FabSetIterator
    :
    public MultiFabIterator
{
public:
    FabSetIterator (FabSet& fabset)
        :
        MultiFabIterator(fabset) {}
    ~FabSetIterator () {}
};

class DependentFabSetIterator
    :
    public DependentMultiFabIterator
{
public:
    DependentFabSetIterator (FabSetIterator& controllerfsiter,
                             FabSet&         fabset)
        : DependentMultiFabIterator(controllerfsiter, fabset) {}

    DependentFabSetIterator (FabSetIterator& controllerfsiter,
                             const FabSet&   fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}

    DependentFabSetIterator (MultiFabIterator& controllerfsiter,
                             FabSet&           fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}

    DependentFabSetIterator (MultiFabIterator& controllerfsiter,
                             const FabSet&     fabset)
        :
        DependentMultiFabIterator(controllerfsiter, fabset) {}

    ~DependentFabSetIterator () {}
};

class ConstFabSetIterator
    :
    public ConstMultiFabIterator
{
public:
    ConstFabSetIterator (const FabSet &fabset)
        :
        ConstMultiFabIterator(fabset) {}

    ~ConstFabSetIterator () {}
};

class ConstDependentFabSetIterator
    :
    public ConstDependentMultiFabIterator
{
public:
    ConstDependentFabSetIterator (ConstFabSetIterator& controllerfsiter,
                                  const FabSet&        fabset)
        :
        ConstDependentMultiFabIterator(controllerfsiter, fabset) {}

    ConstDependentFabSetIterator (ConstMultiFabIterator& controllerfsiter,
                                  const FabSet&          fabset)
        :
        ConstDependentMultiFabIterator(controllerfsiter, fabset) {}

    ~ConstDependentFabSetIterator () {}
};

typedef FabArrayId FabSetId;

class FabSetCopyDescriptor
    :
    public MultiFabCopyDescriptor
{
public:
    FabSetCopyDescriptor (bool cacheremotedata = true)
        :
        MultiFabCopyDescriptor(cacheremotedata) {}

    ~FabSetCopyDescriptor () {}

    FabSetId RegisterFabSet (FabSet* fabset)
    {
        return RegisterMultiFab(fabset);
    }
private:
    //
    // These are disallowed.
    //
    // FabSetCopyDescriptor ();
    FabSetCopyDescriptor (const FabSetCopyDescriptor& rhs);
    FabSetCopyDescriptor& operator= (const FabSetCopyDescriptor& rhs);
};

//
// Inlines.
//

inline
const FabSet&
FabSet::copyTo (FArrayBox& dest) const
{
    this->copy(dest);
    return *this;
}

inline
const FabSet&
FabSet::copyTo (FArrayBox& dest,
                int        src_comp,
                int        dest_comp,
                int        num_comp) const
{
    this->copy(dest, src_comp, dest_comp, num_comp);
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
    this->copy(dest, subbox, src_comp, dest_comp, num_comp);
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
    this->plus(v, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::plus (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    this->plus(v, subreg, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::mult (Real v,
              int  comp,
              int  num_comp)
{
    this->mult(v, comp, num_comp);
    return *this;
}

inline
FabSet&
FabSet::mult (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    this->mult(v, subreg, comp, num_comp);
    return *this;
}

#endif /*_FABSET_H_*/
