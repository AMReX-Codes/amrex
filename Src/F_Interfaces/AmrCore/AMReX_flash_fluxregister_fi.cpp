
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

    void amrex_fi_flash_fluxregister_load_1 (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                             Real* flux, const int* flo, const int* fhi, int nc,
                                             Real scaling_factor)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,flux);
        flux_reg->load(cgid, dir, fab, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_load_2 (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                             Real* flux, const int* flo, const int* fhi, int nc,
                                             Real const* cflux, Real sf_f, Real sf_c)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,flux);
        const FArrayBox cfab(bx,nc,const_cast<Real*>(cflux));
        flux_reg->load(cgid, dir, fab, cfab, sf_f, sf_c);
    }

    void amrex_fi_flash_fluxregister_load_1_area (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                                Real* flux, const int* flo, const int* fhi, int nc,
                                                Real const* area)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,flux);
        const FArrayBox areafab(bx,1,const_cast<Real*>(area));
        flux_reg->load(cgid, dir, fab, areafab);
    }

    void amrex_fi_flash_fluxregister_load_2_area (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                                Real* flux, const int* flo, const int* fhi, int nc,
                                                Real const* cflux, Real const* area,
                                                Real sf_f, Real sf_c)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,flux);
        const FArrayBox cfab(bx,nc,const_cast<Real*>(cflux));
        const FArrayBox areafab(bx,1,const_cast<Real*>(area));
        flux_reg->load(cgid, dir, fab, cfab, areafab, sf_f, sf_c);
    }

    void amrex_fi_flash_fluxregister_load_area_ifd (FlashFluxRegister const* flux_reg, int cgid, int dir,
                                                    Real* flux, const int* flo, const int* fhi, int nc,
                                                    Real const* cflux, Real const* area,
                                                    const int* ifd, Real sf_f, Real sf_c)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        FArrayBox fab(bx,nc,flux);
        const FArrayBox cfab(bx,nc,const_cast<Real*>(cflux));
        const FArrayBox areafab(bx,1,const_cast<Real*>(area));
        flux_reg->load(cgid, dir, fab, cfab, areafab, ifd, sf_f, sf_c);
    }

    void amrex_fi_flash_fluxregister_store (FlashFluxRegister* flux_reg, int fgid, int dir,
                                            Real const* flux, const int* flo, const int* fhi, int nc,
                                            Real scaling_factor)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        const FArrayBox fab(bx,nc,const_cast<Real*>(flux));
        flux_reg->store(fgid, dir, fab, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_store_area (FlashFluxRegister* flux_reg, int fgid, int dir,
                                                 Real const* flux, const int* flo, const int* fhi, int nc,
                                                 Real const* area, Real scaling_factor)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        const FArrayBox fab(bx,nc,const_cast<Real*>(flux));
        const FArrayBox areafab(bx,1,const_cast<Real*>(area));
        flux_reg->store(fgid, dir, fab, areafab, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_store_area_ifd (FlashFluxRegister* flux_reg, int fgid, int dir,
                                                     Real const* flux, const int* flo, const int* fhi, int nc,
                                                     Real const* area, const int* ifd, Real scaling_factor)
    {
        Box bx;
        bx = Box(IntVect(flo), IntVect(fhi));
        bx.shiftHalf(dir,-1);
        const FArrayBox fab(bx,nc,const_cast<Real*>(flux));
        const FArrayBox areafab(bx,1,const_cast<Real*>(area));
        flux_reg->store(fgid, dir, fab, areafab, ifd, scaling_factor);
    }

    void amrex_fi_flash_fluxregister_communicate (FlashFluxRegister* flux_reg)
    {
        flux_reg->communicate();
    }
}
