
#include <WarpX.H>

using namespace amrex;

void
WarpX::InitCrseFineBndryGrids (int ngrow)
{
    const int myproc = ParallelDescriptor::MyProc();

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const BoxArray& ba = grids[lev];
        const DistributionMapping& dm = dmap[lev];
        const Geometry& gm = Geom(lev);

        Box domain = gm.Domain();
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            if (Geometry::isPeriodic(idim)) {
                domain.grow(idim, ngrow);
            }
        }

        BoxList bl(ba.ixType());
        Array<int> iprocs;
        Array<int> grid_idx;

        for (int i = 0, N = ba.size(); i < N; ++i)
        {
            Box bx = ba[i];
            bx.grow(ngrow);
            bx &= domain;

            const BoxList& noncovered = ba.complementIn(bx);
            for (const Box& b : noncovered) {
                bl.push_back(b);
                iprocs.push_back(dm[i]);
                if (dm[i] = myproc) {
                    grid_idx.push(i);
                }
            }
        }

        if (!iprocs.empty()) {
            // define CrseFineBndryGrids
        }
    }
}

void
WarpX::FillBoundaryB ()
{
    FillBoundaryB(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        FillBoundaryB(lev, true);
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
