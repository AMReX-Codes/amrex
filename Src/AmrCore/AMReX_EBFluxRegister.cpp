
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBFluxRegister_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_VisMF.H>

namespace amrex {

EBFluxRegister::EBFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                const IntVect& ref_ratio, int fine_lev, int nvar)
{
    define(fba, cba, fdm, cdm, fgeom, cgeom, ref_ratio, fine_lev, nvar);
}

void
EBFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
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
            const Box& bx = fab.box();

            fab.setVal(crse_cell);
            bool has_fine = false;

            for (const auto& iv : pshifts)
            {
                cfba.intersections(bx+iv, isects, false, 1);
                for (const auto& is : isects)
                {
                    Box ibx = is.second - iv;
                    ibx &= cdomain;
                    fab.setVal(crse_fine_boundary_cell, ibx, 0, 1);
                }

                for (const auto& is : isects)
                {
                    Box ibx = is.second - iv;
                    ibx &= cdomain;
                    ibx &= cfba[is.first];
                    fab.setVal(fine_cell, ibx, 0, 1);
                    has_fine = true;
                }
            }

            if (has_fine) {
                m_crse_fab_flag[mfi.LocalIndex()] = fine_cell;
            }
        }
    }
}


void
EBFluxRegister::CrseSetVal (Real val)
{
    m_crse_data.setVal(val);
}

void
EBFluxRegister::CrseAdd (const MFIter& mfi, const std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt)
{
    BL_ASSERT(mfi.nComp() == flux[0].nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    FArrayBox& fab = m_crse_data[mfi];
    const Box& bx = mfi.tilebox();
    const int nc = fab.nComp();

    const IArrayBox& flag = m_crse_flag[mfi];

    amrex_eb_flux_reg_crseadd(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD(fab),
                              BL_TO_FORTRAN_ANYD(flag),
                              AMREX_D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                                           BL_TO_FORTRAN_ANYD(flux[1]),
                                           BL_TO_FORTRAN_ANYD(flux[2])),
                              dx, &dt,&nc);
                                           
}

}
