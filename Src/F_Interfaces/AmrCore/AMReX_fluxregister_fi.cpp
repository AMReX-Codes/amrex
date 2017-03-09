
#include <AMReX_FluxRegister.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_new_fluxregister (FluxRegister*& flux_reg, const BoxArray* ba, 
                                    const DistributionMapping* dm, int rr, int flev, int ncomp)
    {
        flux_reg = new FluxRegister(*ba, *dm, IntVect(D_DECL(rr,rr,rr)), flev, ncomp);
    }

    void amrex_fi_delete_fluxregister (FluxRegister* flux_reg)
    {
        delete flux_reg;
    }
}
