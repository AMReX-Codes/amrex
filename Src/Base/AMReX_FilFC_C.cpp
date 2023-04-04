#include <AMReX_FilFC_C.H>

namespace amrex {

void fab_filfc (Box const& bx, Array4<Real> const& qn, int ncomp,
                Box const& domain, Real const* /*dx*/, Real const* /*xlo*/,
                BCRec const* bcn)
{
    Box gdomain = domain;
    const auto idxType = bx.ixType();
    const IntVect& len = bx.length();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (bcn->lo(idim) == BCType::int_dir) {
            AMREX_ASSERT(bcn->hi(idim) == BCType::int_dir);
            gdomain.grow(idim, len[idim]);
        }
    }

    FilfcFace fillfunc{};

    // fill on the box faces first
    {
        Array<Box,2*AMREX_SPACEDIM> dom_face_boxes
            = {{ AMREX_D_DECL(amrex::convert(amrex::adjCellLo(gdomain, 0, len[0]),idxType),
                              amrex::convert(amrex::adjCellLo(gdomain, 1, len[1]),idxType),
                              amrex::convert(amrex::adjCellLo(gdomain, 2, len[2]),idxType)),
                 AMREX_D_DECL(amrex::convert(amrex::adjCellHi(gdomain, 0, len[0]),idxType),
                              amrex::convert(amrex::adjCellHi(gdomain, 1, len[1]),idxType),
                              amrex::convert(amrex::adjCellHi(gdomain, 2, len[2]),idxType)) }};

        for (const Box& b : dom_face_boxes) {
            Box tmp = b & bx;
            amrex::LoopOnCpu(tmp, [=] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                IntVect const idx(AMREX_D_DECL(i,j,k));
                fillfunc(idx, qn, 0, ncomp, domain, bcn, 0);
            });
        }
    }

#if (AMREX_SPACEDIM >= 2)
    // fill on the box edges
    {
#if (AMREX_SPACEDIM == 2)
        Array<Box,4> dom_edge_boxes
            = {{ amrex::convert(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),idxType) }};
#else
        Array<Box,12> dom_edge_boxes
            = {{ amrex::convert(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),idxType),
                 //
                 amrex::convert(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),2,len[2]),idxType),
                 //
                 amrex::convert(amrex::adjCellLo(amrex::adjCellLo(gdomain,1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(gdomain,1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(gdomain,1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(gdomain,1,len[1]),2,len[2]),idxType) }};
#endif

        for (const Box& b : dom_edge_boxes) {
            Box tmp = b & bx;
            amrex::LoopOnCpu(tmp, [=] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                IntVect const idx(AMREX_D_DECL(i,j,k));
                fillfunc(idx, qn, 0, ncomp, domain, bcn, 0);
            });
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    // fill on box corners
    {
        Array<Box,8> dom_corner_boxes
            = {{ amrex::convert(amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType),
                 amrex::convert(amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),idxType) }};

        for (const Box& b : dom_corner_boxes) {
            Box tmp = b & bx;
            amrex::LoopOnCpu(tmp, [=] (int i, int j, int k) noexcept
            {
                IntVect const idx(AMREX_D_DECL(i,j,k));
                fillfunc(idx, qn, 0, ncomp, domain, bcn, 0);
            });
        }
    }
#endif
}

}
