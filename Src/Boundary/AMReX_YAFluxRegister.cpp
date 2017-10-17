
#include <AMReX_YAFluxRegister.H>
#include <AMReX_YAFluxRegister_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

YAFluxRegister::YAFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                const IntVect& ref_ratio, int fine_lev, int nvar)
{
    define(fba, cba, fdm, cdm, fgeom, cgeom, ref_ratio, fine_lev, nvar);
}

void
YAFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                        const DistributionMapping& fdm, const DistributionMapping& cdm,
                        const Geometry& fgeom, const Geometry& cgeom,
                        const IntVect& ref_ratio, int fine_lev, int nvar)
{
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ratio = ref_ratio;
    m_fine_level = fine_lev;
    m_ncomp = nvar;

    m_crse_data.define(cba, cdm, nvar, 0);

    m_crse_flag.define(cba, cdm, 1, 1);

    const auto& cperiod = m_crse_geom.periodicity();
    const std::vector<IntVect>& pshifts = cperiod.shiftIntVect();

    BoxArray cfba = fba;
    cfba.coarsen(ref_ratio);

    Box cdomain = m_crse_geom.Domain();
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_crse_geom.isPeriodic(idim)) {
            cdomain.grow(idim,1);
        }
    }

    m_crse_fab_flag.resize(m_crse_flag.local_size(), crse_cell);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;

        for (MFIter mfi(m_crse_flag); mfi.isValid(); ++mfi)
        {
            auto& fab = m_crse_flag[mfi];
            const Box& bx = fab.box() & cdomain;

            fab.setVal(crse_cell);
            bool has_fine = false;

            for (const auto& iv : pshifts)
            {
                cfba.intersections(bx+iv, isects, false, 1);
                for (const auto& is : isects)
                {
                    const Box& ibx = is.second - iv;
                    fab.setVal(crse_fine_boundary_cell, ibx, 0, 1);
                    has_fine = true;
                }
            }

            for (const auto& iv : pshifts)
            {
                cfba.intersections(bx+iv, isects);
                for (const auto& is : isects)
                {
                    const Box& ibx = is.second - iv;
                    fab.setVal(fine_cell, ibx, 0, 1);
                }
            }

            if (has_fine) {
                m_crse_fab_flag[mfi.LocalIndex()] = fine_cell;
            }
        }
    }

    BoxList cfp_bl;
    Vector<int> cfp_procmap;
    int nlocal = 0;
    int myproc = ParallelDescriptor::MyProc();
    
    for (int i = 0; i < cfba.size(); ++i)
    {
        Box bx = amrex::grow(cfba[i], 1);
        bx &= cdomain;

        const BoxList& bl = cfba.complementIn(bx);
        cfp_bl.join(bl);

        int proc = fdm[i];
        for (int j = 0; j < bl.size(); ++j) {
            cfp_procmap.push_back(proc);
            if (proc == myproc) {
                m_cfp_localindex.push_back(nlocal);  // This Array store local index in fine ba/dm.
                                                          // Its size is local size of cfp.
            }
        }

        if (proc == myproc) {
            ++nlocal;
        }
    }

    // It's safe even if cfp_bl is empty.

    BoxArray cfp_ba(std::move(cfp_bl));
    DistributionMapping cfp_dm(cfp_procmap);
    m_cfpatch.define(cfp_ba, cfp_dm, nvar, 0);

    m_cfp_fab.resize(nlocal);
    for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
    {
        const int li = mfi.LocalIndex();
        const int flgi = m_cfp_localindex[li];
        FArrayBox& fab = m_cfpatch[mfi];
        m_cfp_fab[flgi].push_back(&fab);
    }

    bool is_periodic = m_fine_geom.isAnyPeriodic();
    if (is_periodic) {
        m_cfp_mask.define(cfp_ba, cfp_dm, 1, 0);

        const Box& domainbox = m_crse_geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(m_cfp_mask); mfi.isValid(); ++mfi)
            {
                auto& fab = m_cfp_mask[mfi];
                fab.setVal(1.0);
                const Box& bx = fab.box();
                if (!domainbox.contains(bx))  // part of the box is outside periodic boundary
                {
                    for (const auto& iv : pshifts)
                    {
                        if (iv != IntVect::TheZeroVector())
                        {
                            cfba.intersections(bx+iv, isects);
                            for (const auto& is : isects)
                            {
                                const Box& ibx = is.second - iv;
                                fab.setVal(0.0, ibx, 0, 1);
                            }
                        }
                    }
                }
            }
        }
    }
}


void
YAFluxRegister::reset ()
{
    m_crse_data.setVal(0.0);
    m_cfpatch.setVal(0.0);
}


void
YAFluxRegister::CrseAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt)
{
    BL_ASSERT(m_crse_data.nComp() == flux[0]->nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    FArrayBox& fab = m_crse_data[mfi];
    const Box& bx = mfi.tilebox();
    const int nc = fab.nComp();

    const IArrayBox& flag = m_crse_flag[mfi];

    amrex_ya_flux_reg_crseadd(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD(fab),
                              BL_TO_FORTRAN_ANYD(flag),
                              AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*flux[0]),
                                           BL_TO_FORTRAN_ANYD(*flux[1]),
                                           BL_TO_FORTRAN_ANYD(*flux[2])),
                              dx, &dt,&nc);
}


void
YAFluxRegister::FineAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt)
{
    BL_ASSERT(m_cfpatch.nComp() == flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Vector<FArrayBox*>& fabs = m_cfp_fab[li];
    if (fabs.empty()) return;

    BL_ASSERT(mfi.tilebox().cellCentered());
    AMREX_ALWAYS_ASSERT(mfi.tilebox().coarsenable(m_ratio));
    const Box& bx = amrex::coarsen(mfi.tilebox(), m_ratio);
    const int nc = m_cfpatch.nComp();

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
    {
        const Box& lobx = amrex::adjCellLo(bx, idim);
        const Box& hibx = amrex::adjCellHi(bx, idim);
        for (FArrayBox* cfp : m_cfp_fab[li])
        {
            {
                const Box& lobx_is = lobx & cfp->box();
                const int side = 0;
                if (lobx_is.ok())
                {
                    amrex_ya_flux_reg_fineadd(BL_TO_FORTRAN_BOX(lobx_is),
                                              BL_TO_FORTRAN_ANYD(*cfp),
                                              BL_TO_FORTRAN_ANYD(*flux[idim]),
                                              dx, &dt, &nc, &idim, &side,
                                              m_ratio.getVect());
                }
            }
            {
                const Box& hibx_is = hibx & cfp->box();
                const int side = 1;
                if (hibx_is.ok())
                {
                    amrex_ya_flux_reg_fineadd(BL_TO_FORTRAN_BOX(hibx_is),
                                              BL_TO_FORTRAN_ANYD(*cfp),
                                              BL_TO_FORTRAN_ANYD(*flux[idim]),
                                              dx, &dt, &nc, &idim, &side,
                                              m_ratio.getVect());
                }
            }
        }
    }
}


void
YAFluxRegister::Reflux (MultiFab& state, int dc)
{
    if (!m_cfp_mask.empty())
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
        {
            for (int i = 0; i < m_ncomp; ++i) {
                m_cfpatch[mfi].mult(m_cfp_mask[mfi],0,i);
            }
        }
    }

    m_crse_data.ParallelCopy(m_cfpatch, m_crse_geom.periodicity(), FabArrayBase::ADD);

    BL_ASSERT(state.nComp() >= dc + m_ncomp);
    MultiFab::Add(state, m_crse_data, 0, dc, m_ncomp, 0);
}

}
