
#include <AMReX_FlashFluxRegister.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_flash_fluxregister (FlashFluxRegister*& flux_reg,
                                          BoxArray const* fba, BoxArray const* cba,
                                          DistributionMapping const* fdm,
                                          DistributionMapping const* cdm,
                                          Geometry const* fgm, Geometry const* cgm,
                                          int rr, int ncomp)
    {
        flux_reg = new FlashFluxRegister(*fba, *cba, *fdm, *cdm, *fgm, *cgm,
                                         IntVect(AMREX_D_DECL(rr,rr,rr)), ncomp);
    }

    void amrex_fi_delete_flash_fluxregister (FlashFluxRegister* flux_reg)
    {
        delete flux_reg;
    }

    void amrex_fi_flash_fluxregister_load (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                           Real const* fabdata, const int* flo, const int* fhi, int nc,
                                           Real scaling_factor)
    {
        Box bx;
	bx = Box(IntVect(flo), IntVect(fhi));
	bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,fabdata);
        flux_reg->load(cgid, dir, fab, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_store (FlashFluxRegister* flux_reg, int fgid, int dir,
                                            Real* fabdata, const int* flo, const int* fhi, int nc,
                                            Real scaling_factor)
    {
        Box bx;
	bx = Box(IntVect(flo), IntVect(fhi));
	bx.shiftHalf(dir,-1);
        const FArrayBox fab(bx,nc,const_cast<Real*>(fabdata));
        flux_reg->store(fgid, dir, fab, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_communicate (FlashFluxRegister* flux_reg)
    {
        flux_reg->communicate();
    }
}
