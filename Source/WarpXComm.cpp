
#include <WarpX.H>
#include <AMReX_FillPatchUtil.H>
#include <WarpX_f.H>

#include <algorithm>
#include <cstdlib>

using namespace amrex;

void
WarpX::FillBoundaryB ()
{
    FillBoundaryB(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);

        const std::array<MultiFab,BL_SPACEDIM> crse {
                    MultiFab(*Bfield[lev-1][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Bfield[lev-1][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Bfield[lev-1][2], amrex::make_alias, 0, 1) };

        std::array<MultiFab,BL_SPACEDIM> fine {
                    MultiFab(*Bfield[lev][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Bfield[lev][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Bfield[lev][2], amrex::make_alias, 0, 1) };

        amrex::InterpCrseFineBndryEMfield(amrex::InterpB, crse, fine, Geom(lev-1), Geom(lev), ref_ratio[0]);

        FillBoundaryB(lev, true);
    }
}

void
WarpX::FillBoundaryE ()
{
    FillBoundaryE(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);

        const std::array<MultiFab,BL_SPACEDIM> crse {
                    MultiFab(*Efield[lev-1][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Efield[lev-1][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Efield[lev-1][2], amrex::make_alias, 0, 1) };

        std::array<MultiFab,BL_SPACEDIM> fine {
                    MultiFab(*Efield[lev][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Efield[lev][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Efield[lev][2], amrex::make_alias, 0, 1) };

        amrex::InterpCrseFineBndryEMfield(amrex::InterpE, crse, fine, Geom(lev-1), Geom(lev), ref_ratio[0]);

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
    AverageDownEorCurrent(Efield);
}

void
WarpX::AverageDownCurrent ()
{
    AverageDownEorCurrent(current);
}

void
WarpX::AverageDownEorCurrent (Array<std::array<std::unique_ptr<amrex::MultiFab>,3> >& data)
{
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const auto& crse_period = Geom(lev-1).periodicity();
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][1].get(), data[lev][2].get()};
#else
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][2].get()};
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
            data[lev-1][idim]->copy(*crse[idim], crse_period);
        }
#else
        data[lev-1][0]->copy(*crse[0], crse_period);
        data[lev-1][2]->copy(*crse[1], crse_period);
//xxxxx amrex::average_down_nodes(*data[lev][1], *data[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
}

void
WarpX::SyncCurrent ()
{
    for (int lev = finest_level; lev >= 0; --lev)
    {
        const auto& period = Geom(lev).periodicity();
        current[lev][0]->SumBoundary(period);
        current[lev][1]->SumBoundary(period);
        current[lev][2]->SumBoundary(period);
        if (lev > 0) {
            bndry4fine[lev-1]->addTo(*current[lev][0], *current[lev][1], *current[lev][2]);
        }
        if (lev < finest_level) {
            bndry4crse[lev+1]->addTo(*current[lev][0], *current[lev][1], *current[lev][2]);
        }
    }

    AverageDownCurrent();
}


