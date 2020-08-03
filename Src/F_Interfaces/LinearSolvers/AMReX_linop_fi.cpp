
#include <AMReX_MLLinOp.H>

using namespace amrex;

extern "C" {

    void amrex_fi_delete_linop (MLLinOp* linop)
    {
        delete linop;
    }

    void amrex_fi_linop_set_maxorder (MLLinOp* linop, int ord)
    {
        linop->setMaxOrder(ord);
    }

    void amrex_fi_linop_set_domain_bc (MLLinOp* linop, const int* ilobc, const int* ihibc)
    {
        std::array<LinOpBCType,AMREX_SPACEDIM> lobc
        {{AMREX_D_DECL(static_cast<LinOpBCType>(ilobc[0]),
                       static_cast<LinOpBCType>(ilobc[1]),
                       static_cast<LinOpBCType>(ilobc[2]))}};
        std::array<LinOpBCType,AMREX_SPACEDIM> hibc
        {{AMREX_D_DECL(static_cast<LinOpBCType>(ihibc[0]),
                       static_cast<LinOpBCType>(ihibc[1]),
                       static_cast<LinOpBCType>(ihibc[2]))}};
        linop->setDomainBC(lobc,hibc);
    }

    void amrex_fi_linop_set_coarse_fine_bc (MLLinOp* linop, const MultiFab* crse, int crse_ratio)
    {
        linop->setCoarseFineBC(crse, crse_ratio);
    }

    void amrex_fi_linop_set_level_bc (MLLinOp* linop, int amrlev, const MultiFab* levelbcdata)
    {
        linop->setLevelBC(amrlev, levelbcdata);
    }
}
