
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
        if (zeroFirst) {
          flux_reg->FineSetVal(dir, boxno, 0, flux_reg->nComp(), 0.0, RunOn::Cpu);
        }
        const FArrayBox fab(bx, nfluxes, const_cast<Real*>(fabdata));
        flux_reg->FineAdd(fab, dir, boxno, 0, 0, flux_reg->nComp(), scale, RunOn::Cpu);
    }

    void amrex_fi_fluxregister_fineadd_dg
           ( FluxRegister* flux_reg, MultiFab* SurfaceFluxes[],
             int nFields, Real FaceRatio,
             int nDOFX_X1, int nDOFX_X2, int nDOFX_X3,
             int nFineX_X1, int nFineX_X2, int nFineX_X3,
             Real * WeightsX_X1, Real * WeightsX_X2, Real * WeightsX_X3,
             void * vpLX_X1_Refined,
             void * vpLX_X2_Refined,
             void * vpLX_X3_Refined )

    {
        for ( int iDimX = 0; iDimX < BL_SPACEDIM; ++iDimX )
        {
            BL_ASSERT( flux_reg->nComp() == SurfaceFluxes[iDimX]->nComp() );

            flux_reg->FineAdd_DG
              ( *SurfaceFluxes[iDimX], iDimX, nFields, FaceRatio,
                nDOFX_X1, nDOFX_X2, nDOFX_X3,
                nFineX_X1, nFineX_X2, nFineX_X3,
                WeightsX_X1, WeightsX_X2, WeightsX_X3,
                vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined );
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
           ( FluxRegister * flux_reg, MultiFab * SurfaceFluxes[], int nFields,
             int nDOFX_X1, int nDOFX_X2, int nDOFX_X3,
             Real * WeightsX_X1, Real * WeightsX_X2, Real * WeightsX_X3 )
    {

        for ( int iDimX = 0; iDimX < BL_SPACEDIM; ++iDimX )
        {
            BL_ASSERT( flux_reg->nComp() == SurfaceFluxes[iDimX]->nComp() );

            flux_reg->CrseInit_DG
                        ( *SurfaceFluxes[iDimX], iDimX, nFields,
                          nDOFX_X1, nDOFX_X2, nDOFX_X3,
                          WeightsX_X1, WeightsX_X2, WeightsX_X3 );
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

    void amrex_fi_fluxregister_reflux_dg
      ( FluxRegister*   FluxReg,
        MultiFab*       MF_G,
        MultiFab*       MF_dU,
        const Geometry* geom,
        int             nDOFX,
        int             nDOFX_X1,
        int             nDOFX_X2,
        int             nDOFX_X3,
        int             nFields,
        int             iGF_SqrtGm,
        void          * vpNodeNumberTableX_X1,
        void          * vpNodeNumberTableX_X2,
        void          * vpNodeNumberTableX_X3,
        void          * vpWeightsX_q,
        void          * vpLX_X1_Up,
        void          * vpLX_X1_Dn,
        void          * vpLX_X2_Up,
        void          * vpLX_X2_Dn,
        void          * vpLX_X3_Up,
        void          * vpLX_X3_Dn,
        Real          dX1,
        Real          dX2,
        Real          dX3 )
    {

        auto *pNodeNumberTableX_X1
               = reinterpret_cast<int*>(vpNodeNumberTableX_X1);
        Array4<int> NodeNumberTableX_X1
                      ( pNodeNumberTableX_X1, {0,0,0}, {nDOFX,1,1}, 1 );

        auto *pNodeNumberTableX_X2
               = reinterpret_cast<int*>(vpNodeNumberTableX_X2);
        Array4<int> NodeNumberTableX_X2
                      ( pNodeNumberTableX_X2, {0,0,0}, {nDOFX,1,1}, 1 );

        auto *pNodeNumberTableX_X3
               = reinterpret_cast<int*>(vpNodeNumberTableX_X3);
        Array4<int> NodeNumberTableX_X3
                      ( pNodeNumberTableX_X3, {0,0,0}, {nDOFX,1,1}, 1 );

        auto *pWeightsX_q
               = reinterpret_cast<Real*>(vpWeightsX_q);
        Array4<Real> WeightsX_q
                       ( pWeightsX_q, {0,0,0}, {nDOFX,1,1}, 1 );

        auto *pLX_X1_Up
               = reinterpret_cast<Real*>(vpLX_X1_Up);
        Array4<Real> LX_X1_Up
                       ( pLX_X1_Up, {0,0,0}, {nDOFX_X1,nDOFX,1}, 1 );

        auto *pLX_X1_Dn
               = reinterpret_cast<Real*>(vpLX_X1_Dn);
        Array4<Real> LX_X1_Dn
                       ( pLX_X1_Dn, {0,0,0}, {nDOFX_X1,nDOFX,1}, 1 );

        auto *pLX_X2_Up
               = reinterpret_cast<Real*>(vpLX_X2_Up);
        Array4<Real> LX_X2_Up
                       ( pLX_X2_Up, {0,0,0}, {nDOFX_X2,nDOFX,1}, 1 );

        auto *pLX_X2_Dn
               = reinterpret_cast<Real*>(vpLX_X2_Dn);
        Array4<Real> LX_X2_Dn
                       ( pLX_X2_Dn, {0,0,0}, {nDOFX_X2,nDOFX,1}, 1 );

        auto *pLX_X3_Up
               = reinterpret_cast<Real*>(vpLX_X3_Up);
        Array4<Real> LX_X3_Up
                       ( pLX_X3_Up, {0,0,0}, {nDOFX_X3,nDOFX,1}, 1 );

        auto *pLX_X3_Dn
               = reinterpret_cast<Real*>(vpLX_X3_Dn);
        Array4<Real> LX_X3_Dn
                       ( pLX_X3_Dn, {0,0,0}, {nDOFX_X3,nDOFX,1}, 1 );

        FluxReg->Reflux_DG( *MF_G, *MF_dU, *geom,
                            nDOFX, nDOFX_X1, nDOFX_X2, nDOFX_X3,
                            nFields, iGF_SqrtGm,
                            NodeNumberTableX_X1,
                            NodeNumberTableX_X2,
                            NodeNumberTableX_X3,
                            WeightsX_q,
                            LX_X1_Up, LX_X1_Dn,
                            LX_X2_Up, LX_X2_Dn,
                            LX_X3_Up, LX_X3_Dn,
                            dX1, dX2, dX3 );
    }

    void amrex_fi_fluxregister_overwrite (FluxRegister* flux_reg, MultiFab* crse_flxs[],
                                          Real scale, const Geometry* crse_geom)
    {
        const int ncomp = flux_reg->nComp();
        flux_reg->OverwriteFlux({{AMREX_D_DECL(crse_flxs[0], crse_flxs[1], crse_flxs[2])}},
                                scale, 0, 0, ncomp, *crse_geom);
    }
}
