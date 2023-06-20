
#ifndef AMREX_BNDRYREGISTER_H_
#define AMREX_BNDRYREGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_BoxArray.H>
#include <AMReX_FabSet.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Orientation.H>
#include <AMReX_Utility.H>
#include <utility>

namespace amrex {

class Orientation;

/**
* \brief A BndryRegister organizes FabSets bounding each grid in a BoxArray.
*        A FabSet is maintained for each boundary orientation, as well as
*        the BoxArray domain of definition.
*
*        A BndryRegister object contains a list of FabSets bounding the grids
*        in a BoxArray.  The FabSet FABs are at the same refinement level
*        as the grids they bound, and are accessed and modified via a variety
*        of member functions.
*
*        Non-default instantiation allocates a set of FABs, grown into and
*        out of the bounding surfaces of each box in the BoxArray.  The width of
*        the layer (in/out), as well as the "extent" of a bounding FABs (the
*        number of nodes beyond the box boundaries, parallel to the grid
*        surface) are determined by constructor argument.  All boxes and
*        FABs in this context are assumed to be cell-centered.
*
*        A small number of linear mathematical operations are provided for
*        BndryRegisters, as well as overloaded [] operators for access based
*        on grid boundary orientation.  The BoxArray domain of definition is
*        accessible, but not resettable,
*/
template <typename MF>
class BndryRegisterT
{
public:

    using value_type = typename MF::value_type;

    //! The default constructor.
    BndryRegisterT () noexcept = default;

    //! The constructor, given number of cells in/out, extent and number of components (assumes cell-centered boxes, and allocates cell-centered FABs)
    BndryRegisterT (const BoxArray& grids_, // NOLINT(modernize-pass-by-value)
                    const DistributionMapping& dmap,
                    int in_rad, int out_rad, int extent_rad, int ncomp);

    //! The destructor.
    ~BndryRegisterT () = default;

    BndryRegisterT (BndryRegisterT<MF>&& rhs) noexcept = default;

    BndryRegisterT (const BndryRegisterT<MF>& src) = delete;
    BndryRegisterT& operator= (const BndryRegisterT<MF>& src) = delete;
    BndryRegisterT& operator= (BndryRegisterT<MF>&& rhs) = delete;

    void define (const BoxArray& grids_, const DistributionMapping& dmap,
                 int in_rad, int out_rad, int extent_rad, int ncomp);

    //! Build FABs along given face, specifying the DistributionMapping.
    void define (Orientation face, IndexType typ,
                 int in_rad, int out_rad, int extent_rad, int ncomp,
                 const DistributionMapping& dm);

    void clear ();

    //! Get box domain (as an array of boxes).
    const BoxArray& boxes () const noexcept { return grids; }

    //! Return the number of grids in this domain.
    int size () const noexcept { return grids.size(); }

    //! Return const set of FABs bounding the domain grid boxes on a given orientation
    const FabSetT<MF>& operator[] (Orientation face) const noexcept { return bndry[face]; }

    //! Return set of FABs bounding the domain grid boxes on a given orientation
    FabSetT<MF>& operator[] (Orientation face) noexcept { return bndry[face]; }

    //! Set all boundary FABs to given value.
    void setVal (value_type v);

    //! register += rhs
    BndryRegisterT<MF>& operator+= (const BndryRegisterT<MF>& rhs);
    BndryRegisterT<MF>& plus (const BndryRegisterT<MF>& rhs);

    //! Fill the boundary FABs on intersection with given MF.
    BndryRegisterT<MF>& copyFrom (const MF& src, int nghost,
                                  int src_comp, int dest_comp, int num_comp,
                                  const Periodicity& period = Periodicity::NonPeriodic());

    //! Increment the boundary FABs on intersect with given MF.
    BndryRegisterT<MF>& plusFrom (const MF& src, int nghost,
                                  int src_comp, int dest_comp, int num_comp,
                                  const Periodicity& period = Periodicity::NonPeriodic());

    //! Linear combination: this := a*mfa + b*mfb on intersection of MFs with the boundary FABs
    BndryRegisterT<MF>& linComb (value_type a, const MF& mfa, int a_comp,
                                 value_type b, const MF& mfb, int b_comp,
                                 int dest_comp, int num_comp, int n_ghost = 0);

    //! Set box domain, if not set previously.
    void setBoxes (const BoxArray& grids);

    //! Returns constant reference to associated DistributionMapping.
    const DistributionMapping& DistributionMap () const noexcept { return bndry[0].DistributionMap(); }

    //! Write (used for writing to checkpoint)
    void write (const std::string& name, std::ostream& os) const;

    //! Read (used for reading from checkpoint)
    void read (const std::string& name, std::istream& is);

    //! Local copy function
    static void Copy (BndryRegisterT<MF>& dst, const BndryRegisterT<MF>& src);

protected:

    //! The data.
    FabSetT<MF> bndry[2*AMREX_SPACEDIM];
    BoxArray    grids;
};

template <typename MF>
BndryRegisterT<MF>::BndryRegisterT (const BoxArray& grids_, // NOLINT(modernize-pass-by-value)
                                    const DistributionMapping& dmap,
                                    int in_rad, int out_rad,
                                    int extent_rad, int ncomp)
    : grids(grids_)
{
    BL_ASSERT(ncomp > 0);
    BL_ASSERT(grids[0].cellCentered());

    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),in_rad,out_rad,extent_rad,ncomp,dmap);
    }
}

template <typename MF>
void
BndryRegisterT<MF>::define (const BoxArray& grids_,
                            const DistributionMapping& dmap,
                            int in_rad, int out_rad,
                            int extent_rad, int ncomp)
{
    grids = grids_;
    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),in_rad,out_rad,extent_rad,ncomp,dmap);
    }
}

template <typename MF>
void
BndryRegisterT<MF>::clear ()
{
    for (auto & i : bndry) {
        i.clear();
    }
    grids.clear();
}

template <typename MF>
void
BndryRegisterT<MF>::define (Orientation _face, IndexType _typ, int _in_rad,
                            int _out_rad, int _extent_rad, int _ncomp,
                            const DistributionMapping& dmap)
{
    BoxArray fsBA(grids, BATransformer(_face,_typ,_in_rad,_out_rad,_extent_rad));

    FabSetT<MF>& fabs = bndry[_face];

    BL_ASSERT(fabs.size() == 0);

    fabs.define(fsBA,dmap,_ncomp);
    //
    // Go ahead and assign values to the boundary register fabs
    // since in some places APPLYBC (specifically in the tensor
    // operator) the boundary registers are used for a few calculations
    // before the masks are tested to see if you need them.
    //
    fabs.setVal(std::numeric_limits<value_type>::quiet_NaN());
}

template <typename MF>
void
BndryRegisterT<MF>::setBoxes (const BoxArray& _grids)
{
    BL_ASSERT(grids.empty());
    BL_ASSERT(!_grids.empty());
    BL_ASSERT(_grids[0].cellCentered());

    grids = _grids;
    //
    // Check that bndry regions are not allocated.
    //
    for (auto const& k : bndry) {
        amrex::ignore_unused(k);
        BL_ASSERT(k.size() == 0);
    }
}

template <typename MF>
void BndryRegisterT<MF>::setVal (value_type v)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].setVal(v);
    }
}

template <typename MF>
BndryRegisterT<MF>&
BndryRegisterT<MF>::operator+= (const BndryRegisterT<MF>& rhs)
{
    BL_ASSERT(grids == rhs.grids);
    for (OrientationIter face; face; ++face) {
        const auto f = face();
        const int ncomp = bndry[f].nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (FabSetIter bfsi(rhs[f]); bfsi.isValid(); ++bfsi) {
            const Box& bx = bfsi.validbox();
            auto const sfab =   rhs[f].array(bfsi);
            auto       dfab = bndry[f].array(bfsi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dfab(i,j,k,n) += sfab(i,j,k,n);
            });
        }
    }
    return *this;
}

template <typename MF>
BndryRegisterT<MF>&
BndryRegisterT<MF>::plus (const BndryRegisterT<MF>& rhs)
{
    return operator+=(rhs);
}

template <typename MF>
BndryRegisterT<MF>&
BndryRegisterT<MF>::linComb (value_type a, const MF& mfa, int a_comp,
                             value_type b, const MF& mfb, int b_comp,
                             int dest_comp, int num_comp, int n_ghost)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].linComb(a, mfa, a_comp,
                              b, mfb, b_comp,
                              dest_comp, num_comp, n_ghost);
    }
    return *this;
}

template <typename MF>
BndryRegisterT<MF>&
BndryRegisterT<MF>::copyFrom (const MF& src, int nghost,
                              int src_comp, int dest_comp, int num_comp,
                              const Periodicity& period)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyFrom(src,nghost,src_comp,dest_comp,num_comp,period);
    }
    return *this;
}

template <typename MF>
BndryRegisterT<MF>&
BndryRegisterT<MF>::plusFrom (const MF& src, int nghost,
                              int src_comp, int dest_comp, int num_comp,
                              const Periodicity& period)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(src,nghost,src_comp,dest_comp,num_comp,period);
    }
    return *this;
}

template <typename MF>
void
BndryRegisterT<MF>::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        grids.writeOn(os);
        os << '\n';
    }

    for (OrientationIter face; face; ++face)
    {
        //
        // Take name here and make a "new" name unique to each face.
        // Simplest thing would probably to append "_n" to the name,
        // where n is the integer value of face().
        //
        const int i = face();
        BL_ASSERT(i >= 0 && i <= 7);

        std::string facename = amrex::Concatenate(name + '_', i, 1);

        bndry[face()].write(facename);
    }
}

template <typename MF>
void
BndryRegisterT<MF>::read (const std::string& name, std::istream& is)
{
    BoxArray grids_in;

    grids_in.readFrom(is);

    if (!amrex::match(grids,grids_in)) {
        amrex::Abort("BndryRegisterT<MF>::read: grids do not match");
    }

    for (OrientationIter face; face; ++face)
    {
        //
        // Take name here and make a "new" name unique to each face.
        // Simplest thing would probably to append "_n" to the name,
        // where n is the integer value of face().
        //
        const int i = face();
        BL_ASSERT(i >= 0 && i <= 7);

        std::string facename = amrex::Concatenate(name + '_', i, 1);

        bndry[face()].read(facename);
    }
}

// Local copy function
template <typename MF>
void
BndryRegisterT<MF>::Copy (BndryRegisterT<MF>& dst, const BndryRegisterT<MF>& src)
{
    for (OrientationIter face; face; ++face)
    {
        FabSetT<MF>::Copy(dst[face()], src[face()]);
    }
}

using BndryRegister = BndryRegisterT<MultiFab>;
using fBndryRegister = BndryRegisterT<fMultiFab>;

}

#endif /*_BNDRYREGISTER_H_*/
