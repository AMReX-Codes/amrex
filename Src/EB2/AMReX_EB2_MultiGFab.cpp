
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_F.H>

namespace amrex { namespace EB2 {

void
GFab::buildTypes ()
{
    static_assert(sizeof(Type_t) == sizeof(int), "sizeof c_int is not 32");
    
    amrex_eb2_gfab_build_types (BL_TO_FORTRAN_BOX(m_validbox),
                                BL_TO_FORTRAN_ANYD(m_levelset),
                                BL_TO_FORTRAN_ANYD(m_celltype),
                                AMREX_D_DECL(BL_TO_FORTRAN_ANYD(m_facetype[0]),
                                             BL_TO_FORTRAN_ANYD(m_facetype[1]),
                                             BL_TO_FORTRAN_ANYD(m_facetype[2])));
}

}}
