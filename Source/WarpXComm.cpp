
#include <WarpX.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void
WarpX::FillBoundaryB ()
{
    FillBoundaryB(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        amrex::InterpBfield( { MultiFab(*Bfield[lev-1][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                               MultiFab(*Bfield[lev-1][1], amrex::make_alias, 0, 1),
#endif
                               MultiFab(*Bfield[lev-1][2], amrex::make_alias, 0, 1) },
                             { MultiFab(*Bfield[lev  ][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                               MultiFab(*Bfield[lev  ][1], amrex::make_alias, 0, 1),
#endif
                               MultiFab(*Bfield[lev  ][2], amrex::make_alias, 0, 1) },
            Geom(lev));
    }
}

void
WarpX::FillBoundaryE ()
{
    FillBoundaryE(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        FillBoundaryE(lev, true);
    }
}

void
WarpX::FillBoundaryE(int lev, bool force)
{
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Efield[lev][0], *pml_E[0], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][1], *pml_E[1], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][2], *pml_E[2], geom[lev]);

            (*pml_E[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[2]).FillBoundary( geom[lev].periodicity() );
        }

        (*Efield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
}

void
WarpX::FillBoundaryB(int lev, bool force)
{
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Bfield[lev][0], *pml_B[0], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][1], *pml_B[1], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][2], *pml_B[2], geom[lev]);

            (*pml_B[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[2]).FillBoundary( geom[lev].periodicity() );
        }

        (*Bfield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
}

void
WarpX::AverageDownB ()
{
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const Geometry& crse_geom = Geom(lev-1);
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][1].get(), Bfield[lev][2].get()};
#else
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][2].get()};
#endif
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_faces(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            Bfield[lev-1][idim]->copy(*crse[idim], crse_geom.periodicity());
        }
#else
        Bfield[lev-1][0]->copy(*crse[0], crse_geom.periodicity());
        Bfield[lev-1][2]->copy(*crse[1], crse_geom.periodicity());
        amrex::average_down(*Bfield[lev][1], *Bfield[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
}

void
WarpX::AverageDownE ()
{
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const Geometry& crse_geom = Geom(lev-1);
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {Efield[lev][0].get(), Efield[lev][1].get(), Efield[lev][2].get()};
#else
        Array<const MultiFab*> fine {Efield[lev][0].get(), Efield[lev][2].get()};
#endif
        
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_edges(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            Efield[lev-1][idim]->copy(*crse[idim], crse_geom.periodicity());
        }
#else
        Efield[lev-1][0]->copy(*crse[0], crse_geom.periodicity());
        Efield[lev-1][2]->copy(*crse[1], crse_geom.periodicity());
//xxxxx amrex::average_down_nodes(*Efield[lev][1], *Efield[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
}

void
WarpX::AverageDownCurrent ()
{

}

void
WarpX::SyncCurrent ()
{
    for (int lev = 0; lev < finest_level; ++lev)
    {
        
    }
}
