
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_C.H>

namespace amrex { namespace EB2 {

void
GFab::buildTypes (EBCellFlagFab& celltype)
{
    Array4<Real const> const& s = m_levelset.const_array();
    Array4<EBCellFlag> const& cell = celltype.array();
    AMREX_D_TERM(Array4<Type_t> const& fx = m_facetype[0].array();,
                 Array4<Type_t> const& fy = m_facetype[1].array();,
                 Array4<Type_t> const& fz = m_facetype[2].array(););

    const Box& bxg2 = amrex::grow(m_validbox,2);
    const Box& nodal_box = amrex::surroundingNodes(bxg2);

#if (AMREX_SPACEDIM == 2)

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( nodal_box, tbx,
    {
        amrex_eb2_build_types(tbx, bxg2, s, cell, fx, fy);
    });

#elif (AMREX_SPACEDIM == 3)

    Array4<Type_t> const& ex = m_edgetype[0].array();
    Array4<Type_t> const& ey = m_edgetype[1].array();
    Array4<Type_t> const& ez = m_edgetype[2].array();

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( nodal_box, tbx,
    {
        amrex_eb2_build_types(tbx, bxg2, s, cell, fx, fy, fz, ex, ey, ez);
    });

#endif
}

MultiFab
MultiGFab::getLevelSet ()
{
    MultiFab r(amrex::convert(boxArray(),IntVect::TheNodeVector()),
               DistributionMap(), 1, GFab::ng, MFInfo().SetAlloc(false));

    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
        auto& fab = (*this)[mfi].getLevelSet();
        FArrayBox* p = new FArrayBox(fab.box(),1,fab.dataPtr());
        r.setFab(mfi,p);
    }

    return r;
}

}}
