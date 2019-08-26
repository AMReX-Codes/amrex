#include <AMReX_algoim.H>
#include <AMReX_EB2.H>
#include <AMReX_Print.H>
#include <AMReX_algoim_K.H>

namespace amrex { namespace algoim {

void
compute_integrals (MultiFab& intg, int nghost)
{
    compute_integrals(intg, IntVect(nghost));
}

void
compute_integrals (MultiFab& intgmf, IntVect nghost)
{
#if (AMREX_SPACEDIM == 2)
    amrex::Abort("amrex::algoim::compute_integrals is 3D only");
#else

    nghost.min(intgmf.nGrowVect());
    AMREX_ASSERT(intgmf.nComp() >= numIntgs);

    const auto& my_factory = dynamic_cast<EBFArrayBoxFactory const&>(intgmf.Factory());

    // const MultiFab&    vfrac = my_factory.getVolFrac();
    const MultiCutFab& bcent = my_factory.getBndryCent();
    const MultiCutFab& bnorm = my_factory.getBndryNormal();
    const auto&        flags = my_factory.getMultiEBCellFlagFab();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(intgmf,mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        Array4<Real> const& intg = intgmf.array(mfi);

        const auto& flagfab = flags[mfi];
        auto typ = flagfab.getType(bx);

        if (typ == FabType::covered)
        {
            AMREX_HOST_DEVICE_FOR_4D ( bx, numIntgs, i, j, k, n,
            {
                intg(i,j,k,n) = 0.0;
            });
        }
        else if (typ == FabType::regular)
        {
            AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
            {
                set_regular(i,j,k,intg);
            });
        }
        else
        {
            // auto const& vf = vfrac.array(mfi);
            auto const& bc = bcent.array(mfi);
            auto const& bn = bnorm.array(mfi);
            auto const& fg = flagfab.array();

            if (Gpu::inLaunchRegion())
            {
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const auto ebflag = fg(i,j,k);
                    if (ebflag.isRegular()) {
                        set_regular(i,j,k,intg);
                    } else if (ebflag.isCovered()) {
                        for (int n = 0; n < numIntgs; ++n) intg(i,j,k,n) = 0.0;
                    } else {
                        EBPlane phi(bc(i,j,k,0),bc(i,j,k,1),bc(i,j,k,2),
                                    bn(i,j,k,0),bn(i,j,k,1),bn(i,j,k,2));

                        const QuadratureRule q = quadGen(phi);

                        intg(i,j,k,i_S_x    ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x; });
                        intg(i,j,k,i_S_y    ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y; });
                        intg(i,j,k,i_S_z    ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return z; });
                        intg(i,j,k,i_S_x2   ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*x; });
                        intg(i,j,k,i_S_y2   ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y*y; });
                        intg(i,j,k,i_S_z2   ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return z*z; });
                        intg(i,j,k,i_S_x_y  ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*y; });
                        intg(i,j,k,i_S_x_z  ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*z; });
                        intg(i,j,k,i_S_y_z  ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y*z; });
                        intg(i,j,k,i_S_x2_y ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*x*y; });
                        intg(i,j,k,i_S_x2_z ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*x*z; });
                        intg(i,j,k,i_S_x_y2 ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*y*y; });
                        intg(i,j,k,i_S_y2_z ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y*y*z; });
                        intg(i,j,k,i_S_x_z2 ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*z*z; });
                        intg(i,j,k,i_S_y_z2 ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y*z*z; });
                        intg(i,j,k,i_S_x2_y2) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*x*y*y; });
                        intg(i,j,k,i_S_x2_z2) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*x*z*z; });
                        intg(i,j,k,i_S_y2_z2) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return y*y*z*z; });
                        intg(i,j,k,i_S_xyz  ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                   { return x*y*z; });
                    }
                });
            }
            else
            {
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    const auto ebflag = fg(i,j,k);
                    if (ebflag.isRegular()) {
                        set_regular(i,j,k,intg);
                    } else if (ebflag.isCovered()) {
                        for (int n = 0; n < numIntgs; ++n) intg(i,j,k,n) = 0.0;
                    } else {
                        EBPlane phi(bc(i,j,k,0),bc(i,j,k,1),bc(i,j,k,2),
                                    bn(i,j,k,0),bn(i,j,k,1),bn(i,j,k,2));

                        const QuadratureRule q = quadGen(phi);

                        intg(i,j,k,i_S_x    ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x; });
                        intg(i,j,k,i_S_y    ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y; });
                        intg(i,j,k,i_S_z    ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return z; });
                        intg(i,j,k,i_S_x2   ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*x; });
                        intg(i,j,k,i_S_y2   ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y*y; });
                        intg(i,j,k,i_S_z2   ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return z*z; });
                        intg(i,j,k,i_S_x_y  ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*y; });
                        intg(i,j,k,i_S_x_z  ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*z; });
                        intg(i,j,k,i_S_y_z  ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y*z; });
                        intg(i,j,k,i_S_x2_y ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*x*y; });
                        intg(i,j,k,i_S_x2_z ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*x*z; });
                        intg(i,j,k,i_S_x_y2 ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*y*y; });
                        intg(i,j,k,i_S_y2_z ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y*y*z; });
                        intg(i,j,k,i_S_x_z2 ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*z*z; });
                        intg(i,j,k,i_S_y_z2 ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y*z*z; });
                        intg(i,j,k,i_S_x2_y2) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*x*y*y; });
                        intg(i,j,k,i_S_x2_z2) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*x*z*z; });
                        intg(i,j,k,i_S_y2_z2) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return y*y*z*z; });
                        intg(i,j,k,i_S_xyz  ) = q.eval([](Real x, Real y, Real z) noexcept
                                                   { return x*y*z; });
                    }
                }
            }
        }
    }
#endif
}

}}
