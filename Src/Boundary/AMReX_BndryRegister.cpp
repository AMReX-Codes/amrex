
#include <AMReX_BndryRegister.H>
#include <AMReX_Orientation.H>
#include <AMReX_Utility.H>
#include <limits>

namespace amrex {

namespace {
    static const Real BL_SAFE_BOGUS = std::numeric_limits<Real>::quiet_NaN();
}

BndryRegister::BndryRegister () noexcept {}

BndryRegister::~BndryRegister () {}

BndryRegister::BndryRegister (const BoxArray& grids_,
                              const DistributionMapping& dmap,
                              int             in_rad,
                              int             out_rad,
                              int             extent_rad,
                              int             ncomp)
    :
    grids(grids_)
{
    BL_ASSERT(ncomp > 0);
    BL_ASSERT(grids[0].cellCentered());

    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),in_rad,out_rad,extent_rad,ncomp,dmap);
    }
}

void
BndryRegister::define (const BoxArray& grids_,
                       const DistributionMapping& dmap,
                       int             in_rad,
                       int             out_rad,
                       int             extent_rad,
                       int             ncomp)
{
    grids = grids_;
    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),in_rad,out_rad,extent_rad,ncomp,dmap);
    }
}

void
BndryRegister::clear ()
{
    for (int i = 0; i < 2*AMREX_SPACEDIM; ++i) {
        bndry[i].clear();
    }
    grids.clear();
}

void
BndryRegister::init (const BndryRegister& src)
{
    grids = src.grids;

    for (int idim = 0; idim < 2*AMREX_SPACEDIM; idim++)
    {
        const int ncomp = src.bndry[idim].nComp();
        bndry[idim].define(src.bndry[idim].boxArray(), src.DistributionMap(), ncomp);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (FabSetIter mfi(src.bndry[idim]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto const sfab = src.bndry[idim].array(mfi);
            auto       dfab =     bndry[idim].array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dfab(i,j,k,n) = sfab(i,j,k,n);
            });
        }
    }
}

BndryRegister::BndryRegister (const BndryRegister& src)
{
    init(src);
}

BndryRegister&
BndryRegister::operator= (const BndryRegister& src)
{
    if (this != &src)
    {
        if (grids.size() > 0)
        {
            grids.clear();

            for (int i = 0; i < 2*AMREX_SPACEDIM; i++)
                bndry[i].clear();
        }

        init(src);
    }
    return *this;
}

void
BndryRegister::define (Orientation _face,
                       IndexType   _typ,
                       int         _in_rad,
                       int         _out_rad,
                       int         _extent_rad,
                       int         _ncomp,
                       const DistributionMapping& dmap)
{
    BoxArray fsBA(grids, BATransformer(_face,_typ,_in_rad,_out_rad,_extent_rad));

    FabSet& fabs = bndry[_face];

    BL_ASSERT(fabs.size() == 0);

    fabs.define(fsBA,dmap,_ncomp);
    //
    // Go ahead and assign values to the boundary register fabs
    // since in some places APPLYBC (specifically in the tensor
    // operator) the boundary registers are used for a few calculations
    // before the masks are tested to see if you need them.
    //
    fabs.setVal(BL_SAFE_BOGUS);
}

void
BndryRegister::setBoxes (const BoxArray& _grids)
{
    BL_ASSERT(grids.size() == 0);
    BL_ASSERT(_grids.size() > 0);
    BL_ASSERT(_grids[0].cellCentered());

    grids = _grids;
    //
    // Check that bndry regions are not allocated.
    //
    for (int k = 0; k < 2*AMREX_SPACEDIM; k++)
        BL_ASSERT(bndry[k].size() == 0);
}

void BndryRegister::setVal (Real v)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].setVal(v);
    }
}

BndryRegister&
BndryRegister::operator+= (const BndryRegister& rhs)
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

BndryRegister&
BndryRegister::plus (const BndryRegister& rhs)
{
    return operator+=(rhs);
}

BndryRegister&
BndryRegister::linComb (Real            a,
                        const MultiFab& mfa,
                        int             a_comp,
                        Real            b,
                        const MultiFab& mfb,
                        int             b_comp,
                        int             dest_comp,
                        int             num_comp,
                        int             n_ghost)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].linComb(a,
                              mfa,
                              a_comp,
                              b,
                              mfb,
                              b_comp,
                              dest_comp,
                              num_comp,
                              n_ghost);
    }
    return *this;
}

BndryRegister&
BndryRegister::copyFrom (const MultiFab& src,
                         int             nghost,
                         int             src_comp,
                         int             dest_comp,
                         int             num_comp,
                         const Periodicity& period)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyFrom(src,nghost,src_comp,dest_comp,num_comp,period);
    }
    return *this;
}

BndryRegister&
BndryRegister::plusFrom (const MultiFab& src,
                         int             nghost,
                         int             src_comp,
                         int             dest_comp,
                         int             num_comp,
                         const Periodicity& period)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(src,nghost,src_comp,dest_comp,num_comp,period);
    }
    return *this;
}

void
BndryRegister::write (const std::string& name, std::ostream& os) const
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

void
BndryRegister::read (const std::string& name, std::istream& is)
{
    BoxArray grids_in;

    grids_in.readFrom(is);

    if (!amrex::match(grids,grids_in)) {
        amrex::Abort("BndryRegister::read: grids do not match");
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
void
BndryRegister::Copy (BndryRegister& dst, const BndryRegister& src)
{
    for (OrientationIter face; face; ++face)
    {
        FabSet::Copy(dst[face()], src[face()]);
    }
}

}

