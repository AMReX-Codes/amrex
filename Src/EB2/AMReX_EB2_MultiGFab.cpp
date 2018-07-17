
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_F.H>

namespace amrex { namespace EB2 {

void
GFab::buildTypes (EBCellFlagFab& celltype)
{
    static_assert(sizeof(Type_t) == sizeof(int), "sizeof c_int is not 32");
    
#if (AMREX_SPACEDIM == 2)
    amrex_eb2_gfab_build_types (BL_TO_FORTRAN_BOX(m_validbox),
                                BL_TO_FORTRAN_ANYD(m_levelset),
                                BL_TO_FORTRAN_ANYD(celltype),
                                BL_TO_FORTRAN_ANYD(m_facetype[0]),
                                BL_TO_FORTRAN_ANYD(m_facetype[1]));
#elif (AMREX_SPACEDIM == 3)
    amrex_eb2_gfab_build_types (BL_TO_FORTRAN_BOX(m_validbox),
                                BL_TO_FORTRAN_ANYD(m_levelset),
                                BL_TO_FORTRAN_ANYD(celltype),
                                BL_TO_FORTRAN_ANYD(m_facetype[0]),
                                BL_TO_FORTRAN_ANYD(m_facetype[1]),
                                BL_TO_FORTRAN_ANYD(m_facetype[2]),
                                BL_TO_FORTRAN_ANYD(m_edgetype[0]),
                                BL_TO_FORTRAN_ANYD(m_edgetype[1]),
                                BL_TO_FORTRAN_ANYD(m_edgetype[2]));
#endif
}

MultiFab
MultiGFab::getLevelSet ()
{
    Vector<Real*> p;
    p.reserve(local_size());
    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
        p.push_back((*this)[mfi].getLevelSet().dataPtr());
    }
    return MultiFab(amrex::convert(boxArray(),IntVect::TheNodeVector()),
                    DistributionMap(), 1, IntVect(GFab::ng), p);
}

}}
