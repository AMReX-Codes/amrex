#include "initProb_K.H"

#include "MyTest.H"

using namespace amrex;

void
MyTest::initProb ()
{
    const auto prob_lo = geom.ProbLoArray();
    const auto dx      = geom.CellSizeArray();
    const auto a = alpha;
    const auto b = beta;
    const auto ndhi = geom.Domain().bigEnd()+1;

    enum class CoordID { Cartesian, Cyl1D, Cyl2D, Sph1D };
    auto cid = CoordID::Cartesian;
#if (AMREX_SPACEDIM == 1)
    if (this->coord == 1) {
        cid = CoordID::Cyl1D;
        AMREX_ALWAYS_ASSERT(prob_lo[0] == 0);
    } else if (this->coord == 2) {
        cid = CoordID::Sph1D;
        AMREX_ALWAYS_ASSERT(prob_lo[0] == 0);
    }
#elif (AMREX_SPACEDIM == 2)
    if (this->coord == 1) {
        amrex::Abort("2D Cylindrical support will be added later");
    }
#endif

    amrex::ignore_unused(cid,ndhi);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.tilebox(IntVect(1),IntVect(1));
        GpuArray<Array4<Real>,3> rhsfab{rhs[0].array(mfi),
                                        rhs[1].array(mfi),
                                        rhs[2].array(mfi)};
        GpuArray<Array4<Real>,3> solfab{solution[0].array(mfi),
                                        solution[1].array(mfi),
                                        solution[2].array(mfi)};
        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 1)
            if (cid == CoordID::Sph1D) {
                actual_init_prob_sph1d(i,j,k,rhsfab,solfab,prob_lo,dx,a,b,ndhi);
            } else if (cid == CoordID::Cyl1D) {
                actual_init_prob_cyl1d(i,j,k,rhsfab,solfab,prob_lo,dx,a,b,ndhi);
            } else
#endif
            {
                actual_init_prob(i,j,k,rhsfab,solfab,prob_lo,dx,a,b);
            }
        });
    }
}
