
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

    void amrex_fi_fluxregister_fineadd_dg
      ( FluxRegister* flux_reg, MultiFab* SurfaceFluxes[],
        int nFields, int nDOFX_X1, int nDOFX_X2, int nDOFX_X3,
        Real WeightsX_X1[], Real WeightsX_X2[], Real WeightsX_X3[],
        Real LX_X1[], Real LX_X2[], Real LX_X3[] )
    {
        for (int iDimX = 0; iDimX < BL_SPACEDIM; ++iDimX)
        {
            BL_ASSERT( flux_reg->nComp() == SurfaceFluxes[iDimX]->nComp() );

            if( iDimX == 0 )
            {
              flux_reg->FineAdd_DG( *SurfaceFluxes[iDimX], iDimX,
                                    nFields, nDOFX_X1, WeightsX_X1, LX_X1 );
            }
            else if( iDimX == 1 )
            {
              flux_reg->FineAdd_DG( *SurfaceFluxes[iDimX], iDimX,
                                    nFields, nDOFX_X2, WeightsX_X2, LX_X2 );
            }
            else if( iDimX == 2 )
            {
              flux_reg->FineAdd_DG( *SurfaceFluxes[iDimX], iDimX,
                                    nFields, nDOFX_X3, WeightsX_X3, LX_X3 );
            }
            else
            {
                std::cout<< "Src/F_Interfaces/AmrCore/AMReX_fluxregister_fi.cpp"
                         << "\namrex_fi_fluxregister_fineadd_dg"
                         << "\nThis should never print.\n";
            }
        }
    } /* END void amrex_fi_fluxregister_fineadd_dg */

    void amrex_fi_fluxregister_crseinit (FluxRegister* flux_reg, MultiFab* flxs[], Real scale)
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            BL_ASSERT(flux_reg->nComp() == flxs[dir]->nComp());
            flux_reg->CrseInit(*flxs[dir], dir, 0, 0, flux_reg->nComp(), scale);
        }
    }

    void amrex_fi_fluxregister_crseinit_dg
           ( FluxRegister* flux_reg, MultiFab* SurfaceFluxes[],
             int nFields, int nDOFX_X1, int nDOFX_X2, int nDOFX_X3,
             Real WeightsX_X1[], Real WeightsX_X2[], Real WeightsX_X3[] )
    {

        for ( int iDimX = 0; iDimX < BL_SPACEDIM; ++iDimX )
        {
            BL_ASSERT( flux_reg->nComp() == SurfaceFluxes[iDimX]->nComp() );

            if( iDimX == 0 )
            {
                flux_reg->CrseInit_DG
                            ( *SurfaceFluxes[iDimX], iDimX, nFields, nDOFX_X1,
                              WeightsX_X1 );
            }
            else if( iDimX == 1 )
            {
                flux_reg->CrseInit_DG
                            ( *SurfaceFluxes[iDimX], iDimX, nFields, nDOFX_X2,
                              WeightsX_X2 );
            }
            else if( iDimX == 2 )
            {
                flux_reg->CrseInit_DG
                            ( *SurfaceFluxes[iDimX], iDimX, nFields, nDOFX_X3,
                              WeightsX_X3 );
            }
            else
            {
                std::cout<< "Src/F_Interfaces/AmrCore/AMReX_fluxregister_fi.cpp"
                         << "\namrex_fi_fluxregister_crseinit_dg"
                         << "\nThis should never print.\n";
            }
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

std::cout<<mf->nComp()<<std::endl;

        BL_ASSERT(flux_reg->nComp() == mf->nComp());
        flux_reg->Reflux(*mf, vol, scale, 0, 0, flux_reg->nComp(), *geom);
    }

    void amrex_fi_fluxregister_reflux_dg
      ( FluxRegister* flux_reg, MultiFab* mf,
        const Geometry* geom, int nFields, int nDOFX, int nNodesX[],
        Real WeightsX_q[], Real dX1[], Real dX2[], Real dX3[],
        Real LX_X1_Dn[], Real LX_X1_Up[],
        Real LX_X2_Dn[], Real LX_X2_Up[],
        Real LX_X3_Dn[], Real LX_X3_Up[] )
    {
        flux_reg->Reflux_DG( *mf, *geom, flux_reg->nComp(),
                             nFields, nDOFX, nNodesX,
                             WeightsX_q, dX1, dX2, dX3,
                             LX_X1_Dn, LX_X1_Up,
                             LX_X2_Dn, LX_X2_Up,
                             LX_X3_Dn, LX_X3_Up );
    }

    void amrex_fi_fluxregister_overwrite (FluxRegister* flux_reg, MultiFab* crse_flxs[],
                                          Real scale, const Geometry* crse_geom)
    {
        const int ncomp = flux_reg->nComp();
        flux_reg->OverwriteFlux({{AMREX_D_DECL(crse_flxs[0], crse_flxs[1], crse_flxs[2])}},
                                scale, 0, 0, ncomp, *crse_geom);
    }
}
