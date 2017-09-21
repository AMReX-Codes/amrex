
#include <AMReX_BCUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_filcc_f.H>

namespace amrex
{
    void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Array<BCRec>& bc)
    {
        if (Geometry::isAllPeriodic()) return;
        if (phi.nGrow() == 0) return;

        AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

        const Box& domain_box = geom.Domain();
        Box grown_domain_box = domain_box;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (Geometry::isPeriodic(idim)) {
                grown_domain_box.grow(idim,phi.nGrow());
            }
        }
        // Inside grown_domain_box, we have good data.
    
        const Real* dx = geom.CellSizeF();
        const Real* prob_lo = Geometry::ProbLo();

#ifdef AMREX_USE_CUDA
        int n_bc = AMREX_SPACEDIM * 2 * phi.nComp();
        int* bcrs_d = (int*) Device::device_malloc(n_bc * sizeof(int));
        const int* bcrs = bc[0].data();
        Device::device_htod_memcpy_async(bcrs_d, bcrs, n_bc * sizeof(int), -1);
#else
        const int* bcrs_d = bc[0].data();
#endif

        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = phi[mfi];
            const Box& fab_box = fab.box(); // including ghost cells

#ifdef AMREX_USE_DEVICE
            Device::prepare_for_launch(fab_box.loVect(), fab_box.hiVect());
#endif
            
            if (! grown_domain_box.contains(fab_box))
            {
                amrex_fab_filcc(BL_TO_FORTRAN_FAB(fab),
                                BL_TO_FORTRAN_BOX(domain_box),
                                dx, prob_lo,
                                bcrs_d);
            }
        }
    }
}
