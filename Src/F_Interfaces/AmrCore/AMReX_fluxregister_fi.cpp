
#include <AMReX_FluxRegister.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_new_fluxregister (FluxRegister*& flux_reg, const BoxArray* ba, 
                                    const DistributionMapping* dm, int rr, int flev, int ncomp)
    {
        flux_reg = new FluxRegister(*ba, *dm, IntVect(AMREX_D_DECL(rr,rr,rr)), flev, ncomp);
    }

    void amrex_fi_delete_fluxregister (FluxRegister* flux_reg)
    {
        delete flux_reg;
    }

    void amrex_fi_fluxregister_fineadd (FluxRegister* flux_reg, MultiFab* flxs[], Real scale)
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            BL_ASSERT(flux_reg->nComp() == flxs[dir]->nComp());
            flux_reg->FineAdd(*flxs[dir], dir, 0, 0, flux_reg->nComp(), scale);
        }
    }

  void amrex_fi_fluxregister_fineadd_1fab_1dir (FluxRegister* flux_reg, const Real* fabdata,  const int* flo, const int* fhi, int dir, int boxno, int zeroFirst, int nfluxes, Real scale)
    {
        Box bx;
	bx = Box(IntVect(flo), IntVect(fhi));
	bx.shiftHalf(dir,-1);

	BL_ASSERT(flux_reg->nComp() == nfluxes);
	if (zeroFirst)
	  flux_reg->FineSetVal(dir, boxno, 0, flux_reg->nComp(), 0.0, RunOn::Cpu);
        const FArrayBox fab(bx, nfluxes, const_cast<Real*>(fabdata));
	flux_reg->FineAdd(fab, dir, boxno, 0, 0, flux_reg->nComp(), scale, RunOn::Cpu);
    }

    void amrex_fi_fluxregister_crseinit (FluxRegister* flux_reg, MultiFab* flxs[], Real scale)
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            BL_ASSERT(flux_reg->nComp() == flxs[dir]->nComp());
            flux_reg->CrseInit(*flxs[dir], dir, 0, 0, flux_reg->nComp(), scale);
        }
    }

    void amrex_fi_fluxregister_crseadd (FluxRegister* flux_reg, MultiFab* flxs[], Real scale,
                                        const Geometry* geom)
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            BL_ASSERT(flux_reg->nComp() == flxs[dir]->nComp());
            flux_reg->CrseAdd(*flxs[dir], dir, 0, 0, flux_reg->nComp(), scale, *geom);
        }
    }

    void amrex_fi_fluxregister_setval (FluxRegister* flux_reg, Real val)
    {
        flux_reg->setVal(val);
    }

    void amrex_fi_fluxregister_reflux (FluxRegister* flux_reg, MultiFab* mf, Real scale, 
                                       const Geometry* geom)
    {
        MultiFab vol;
        geom->GetVolume(vol, mf->boxArray(), mf->DistributionMap(), 0);
        BL_ASSERT(flux_reg->nComp() == mf->nComp());
        flux_reg->Reflux(*mf, vol, scale, 0, 0, flux_reg->nComp(), *geom);
    }

    void amrex_fi_fluxregister_overwrite (FluxRegister* flux_reg, MultiFab* crse_flxs[],
                                          Real scale, const Geometry* crse_geom)
    {
        const int ncomp = flux_reg->nComp();
        flux_reg->OverwriteFlux({AMREX_D_DECL(crse_flxs[0], crse_flxs[1], crse_flxs[2])},
                                scale, 0, 0, ncomp, *crse_geom);
    }
}
