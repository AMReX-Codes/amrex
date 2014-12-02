
#ifndef _FABSET_H_
#define _FABSET_H_

#include <MultiFab.H>

/*
        A FabSet is a group of FArrayBox's.  The grouping is designed
        specifically to represent regions along the boundary of Box's,
        and are used to implement boundary conditions to discretized
        partial differential equations.

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
    // The default constructor -- you must later call define().
    //
    FabSet ();
    //
    // Construct a FabSet of specified number of components on the grids.
    //
    FabSet (const BoxArray& grids, int ncomp);
    //
    // The destructor.
    //
    ~FabSet ();
    //
    // Define a FabSet constructed via default constructor.
    //
    void define (const BoxArray& grids, int ncomp);
    //
    // Define a FabSet constructed via default constructor.
    //
    void define (const BoxArray&            grids,
                 int                        ncomp,
                 const DistributionMapping& dm);
    //
    // Returns number of variables associated with each point.
    //
    using MultiFab::nComp;
    //
    // Returns number of Fabs owned by this MPI rank
    //
    using MultiFab::local_size;
    //
    // Returns a constant reference to the BoxArray.
    //
    using MultiFab::boxArray;
    //
    // Returns a constant reference to the Kth Box in the BoxArray.
    //
    using MultiFab::box;
    using MultiFab::fabbox;
    //
    // Promote all setVal() member functions back to public.
    //
    using MultiFab::setVal;
    //
    // Copy from this FabSet to the destination FAB.
    //
    const FabSet& copyTo (FArrayBox& dest) const;
    //
    // Copy specified components of this FabSet to the components of the destination FAB
    //
    const FabSet& copyTo (FArrayBox& dest,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    // Copy specified components of this FabSet to the components of the destination FAB over the subregion
    //
    const FabSet& copyTo (FArrayBox& dest,
                          const Box& subbox,
                          int        src_comp,
                          int        dest_comp,
                          int        num_comp) const;
    //
    // Copy from the FAB to this FabSet.
    //
    FabSet& copyFrom (const FArrayBox& src);
    //
    // Copy from the specified components of a source FAB to the components to this destination FabSet
    //
    FabSet& copyFrom (const FArrayBox& src,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    // Copy from the specified components of a source FAB to the destination components of this FabSet over the subregion
    //
    FabSet& copyFrom (const FArrayBox& src,
                      const Box&       subbox,
                      int              src_comp,
                      int              dest_comp,
                      int              num_comp);
    //
    // Copy from the source MultiFab to this destination FabSet.
    //
    FabSet& copyFrom (const MultiFab& src,
                      int             nghost,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp);
    //
    // Copy from the FabSet to this FabSet.
    //
    FabSet& copyFrom (const FabSet& src);
    //
    // Copy from the specified components of a source FabSet to the components to this destination FabSet
    //
    FabSet& copyFrom (const FabSet& src,
                      int           src_comp,
                      int           dest_comp,
                      int           num_comp);
    //
    // Copy to the MultiFab from this FabSet.
    //
    void copyTo (MultiFab& dest) const;
    //
    // Componentwise multiply each datum in this FabSet for specified components.
    //
    FabSet& mult (Real v,
                  int  comp,
                  int  num_comp);
    //
    // Componentwise multiply each datum in this FabSet for specified components in the subregion.
    //
    FabSet& mult (Real       v,
                  const Box& subreg,
                  int        comp,
                  int        num_comp);
    //
    // Componentwise add scalar to each datum in this FabSet for specified components.
    //
    FabSet& plus (Real v,
                  int  comp,
                  int  num_comp);
    //
    // Componentwise add scalar to each datum in this FabSet for specified components in the subregion.
    //
    FabSet& plus (Real       v,
                  const Box& subreg,
                  int        comp,
                  int        num_comp);
    //
    // Add each datum in MultiFab to each in this FabSet for specified components in the subregion, including ghost cells specified.
    //
    FabSet& plusFrom (const MultiFab& src,
                      int             nghost,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp);
    //
    // Linear combination: this := a*this + b*src (FabSets must be commensurate).
    //
    FabSet& linComb (Real          a,
                     Real          b,
                     const FabSet& src,
                     int           src_comp,
                     int           dest_comp,
                     int           num_comp);
    //
    // Linear combination: this := a*mfa + b*mfb on intersection with valid region.
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

    friend class FabSetIter;
    friend class FabSetCopyDescriptor;
    //
    // Promote operator[] back to public.
    //
    using MultiFab::operator[];

    using MultiFab::get;
    //
    // Promote clear() back to public.
    //
    using MultiFab::clear;
    //
    // Promote DistributionMap() back to public.
    //
    using MultiFab::DistributionMap;

    using MultiFab::IndexMap;
    //
    // Promote size() back to public.
    //
    using MultiFab::size;
    //
    // True if there is a FAB at specified index.
    //
    bool defined (int i) const { return this->MultiFab::defined(i); }
    //
    // Write (used for writing to checkpoint)
    //
    void write (const std::string& name) const;
    //
    // Read (used for reading from checkpoint)
    //
    void read (const std::string& name);

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
//  An iterator over the FABs in a FabSet.
//
//  This class encapsulates an iterator over the FABs in a FabSet such that
//  each CPU in a Single Program Multiple Data (SPMD) parallel program accesses
//  only those FABs in the FabSet that are actually owned by that CPU.
//

class FabSetIter
    :
    public MFIter
{
public:

    FabSetIter (const FabSet& fabset)
        :
        MFIter(fabset) {}

    FabSetIter (const FabSet& fabset, bool do_tiling)
        :
        MFIter(fabset,do_tiling) {}

    FabSetIter (const FabSet& fabset, const IntVect& tilesize)
        :
        MFIter(fabset,tilesize) {}
};

typedef FabArrayId FabSetId;

class FabSetCopyDescriptor
    :
    public MultiFabCopyDescriptor
{
public:

    FabSetCopyDescriptor ();

    ~FabSetCopyDescriptor ();

    FabSetId RegisterFabSet (FabSet* fabset) { return RegisterMultiFab(fabset); }

private:
    //
    // These are disallowed.
    //
    FabSetCopyDescriptor (const FabSetCopyDescriptor& rhs);
    FabSetCopyDescriptor& operator= (const FabSetCopyDescriptor& rhs);
};

#endif /*_FABSET_H_*/
