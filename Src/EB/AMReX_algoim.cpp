#include <AMReX_algoim.H>
#include <AMReX_EB2.H>
#include <AMReX_Print.H>
#include <AMReX_algoim_K.H>

namespace amrex::algoim {

void
compute_integrals (MultiFab& intg, int nghost)
{
    compute_integrals(intg, IntVect(nghost));
}

void
compute_integrals (MultiFab& intgmf, IntVect nghost)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(intgmf, nghost);
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

#ifdef AMREX_USE_OMP
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

                        intg(i,j,k,i_S_x    ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real /*z*/) noexcept
                                                   { return x; });
                        intg(i,j,k,i_S_y    ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real /*z*/) noexcept
                                                   { return y; });
                        intg(i,j,k,i_S_z    ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real /*y*/, Real z) noexcept
                                                   { return z; });
                        intg(i,j,k,i_S_x2   ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real /*z*/) noexcept
                                                   { return x*x; });
                        intg(i,j,k,i_S_y2   ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real /*z*/) noexcept
                                                   { return y*y; });
                        intg(i,j,k,i_S_z2   ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real /*y*/, Real z) noexcept
                                                   { return z*z; });
                        intg(i,j,k,i_S_x_y  ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real /*z*/) noexcept
                                                   { return x*y; });
                        intg(i,j,k,i_S_x_z  ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real z) noexcept
                                                   { return x*z; });
                        intg(i,j,k,i_S_y_z  ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real z) noexcept
                                                   { return y*z; });
                        intg(i,j,k,i_S_x2_y ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real /*z*/) noexcept
                                                   { return x*x*y; });
                        intg(i,j,k,i_S_x2_z ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real z) noexcept
                                                   { return x*x*z; });
                        intg(i,j,k,i_S_x_y2 ) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real /*z*/) noexcept
                                                   { return x*y*y; });
                        intg(i,j,k,i_S_y2_z ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real z) noexcept
                                                   { return y*y*z; });
                        intg(i,j,k,i_S_x_z2 ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real z) noexcept
                                                   { return x*z*z; });
                        intg(i,j,k,i_S_y_z2 ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real z) noexcept
                                                   { return y*z*z; });
                        intg(i,j,k,i_S_x2_y2) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real /*z*/) noexcept
                                                   { return x*x*y*y; });
                        intg(i,j,k,i_S_x2_z2) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real z) noexcept
                                                   { return x*x*z*z; });
                        intg(i,j,k,i_S_y2_z2) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real z) noexcept
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

                        intg(i,j,k,i_S_x    ) = q.eval([](Real x, Real /*y*/, Real /*z*/) noexcept
                                                   { return x; });
                        intg(i,j,k,i_S_y    ) = q.eval([](Real /*x*/, Real y, Real /*z*/) noexcept
                                                   { return y; });
                        intg(i,j,k,i_S_z    ) = q.eval([](Real /*x*/, Real /*y*/, Real z) noexcept
                                                   { return z; });
                        intg(i,j,k,i_S_x2   ) = q.eval([](Real x, Real /*y*/, Real /*z*/) noexcept
                                                   { return x*x; });
                        intg(i,j,k,i_S_y2   ) = q.eval([](Real /*x*/, Real y, Real /*z*/) noexcept
                                                   { return y*y; });
                        intg(i,j,k,i_S_z2   ) = q.eval([](Real /*x*/, Real /*y*/, Real z) noexcept
                                                   { return z*z; });
                        intg(i,j,k,i_S_x_y  ) = q.eval([](Real x, Real y, Real /*z*/) noexcept
                                                   { return x*y; });
                        intg(i,j,k,i_S_x_z  ) = q.eval([](Real x, Real /*y*/, Real z) noexcept
                                                   { return x*z; });
                        intg(i,j,k,i_S_y_z  ) = q.eval([](Real /*x*/, Real y, Real z) noexcept
                                                   { return y*z; });
                        intg(i,j,k,i_S_x2_y ) = q.eval([](Real x, Real y, Real /*z*/) noexcept
                                                   { return x*x*y; });
                        intg(i,j,k,i_S_x2_z ) = q.eval([](Real x, Real /*y*/, Real z) noexcept
                                                   { return x*x*z; });
                        intg(i,j,k,i_S_x_y2 ) = q.eval([](Real x, Real y, Real /*z*/) noexcept
                                                   { return x*y*y; });
                        intg(i,j,k,i_S_y2_z ) = q.eval([](Real /*x*/, Real y, Real z) noexcept
                                                   { return y*y*z; });
                        intg(i,j,k,i_S_x_z2 ) = q.eval([](Real x, Real /*y*/, Real z) noexcept
                                                   { return x*z*z; });
                        intg(i,j,k,i_S_y_z2 ) = q.eval([](Real /*x*/, Real y, Real z) noexcept
                                                   { return y*z*z; });
                        intg(i,j,k,i_S_x2_y2) = q.eval([](Real x, Real y, Real /*z*/) noexcept
                                                   { return x*x*y*y; });
                        intg(i,j,k,i_S_x2_z2) = q.eval([](Real x, Real /*y*/, Real z) noexcept
                                                   { return x*x*z*z; });
                        intg(i,j,k,i_S_y2_z2) = q.eval([](Real /*x*/, Real y, Real z) noexcept
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

void
compute_surface_integrals (MultiFab& sintg, int nghost)
{
    compute_surface_integrals(sintg, IntVect(nghost));
}

void
compute_surface_integrals (MultiFab& sintgmf, IntVect nghost)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(sintgmf, nghost);
    amrex::Abort("amrex::algoim::compute_surface_integrals is 3D only");
#else

    nghost.min(sintgmf.nGrowVect());
    AMREX_ASSERT(sintgmf.nComp() >= numSurfIntgs);

    const auto& my_factory = dynamic_cast<EBFArrayBoxFactory const&>(sintgmf.Factory());

    const MultiFab&    vfrac = my_factory.getVolFrac();
    const MultiCutFab& bcent = my_factory.getBndryCent();
    const MultiCutFab& bnorm = my_factory.getBndryNormal();
    const auto&        flags = my_factory.getMultiEBCellFlagFab();
    const auto&        area  = my_factory.getAreaFrac();
    const auto&        barea = my_factory.getBndryArea();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sintgmf,mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        Array4<Real> const& sintg = sintgmf.array(mfi);

        const auto& flagfab = flags[mfi];
        auto typ = flagfab.getType(bx);

        if (typ == FabType::covered)
        {
            AMREX_HOST_DEVICE_FOR_4D ( bx, numSurfIntgs, i, j, k, n,
            {
                sintg(i,j,k,n) = 0.0;
            });
        }
        else if (typ == FabType::regular)
        {
            AMREX_HOST_DEVICE_FOR_3D ( bx, i, j, k,
            {
                set_regular_surface(i,j,k,sintg);
            });
        }
        else
        {
            auto const& vf  = vfrac.array(mfi);
            auto const& bc  = bcent.array(mfi);
            auto const& bn  = bnorm.array(mfi);
            auto const& fg  = flagfab.array();
            auto const& apx = area[0]->const_array(mfi);
            auto const& apy = area[1]->const_array(mfi);
            auto const& apz = area[2]->const_array(mfi);
            auto const& ba  = barea.array(mfi);

            if (Gpu::inLaunchRegion())
            {
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const auto ebflag = fg(i,j,k);
                    if (ebflag.isRegular()) {
                        set_regular_surface(i,j,k,sintg);
                    } else if (ebflag.isCovered()) {
                        for (int n = 0; n < numSurfIntgs; ++n) sintg(i,j,k,n) = 0.0;
                    } else {
                        constexpr Real almostone = Real(1.) - Real(100.)*std::numeric_limits<Real>::epsilon();

                        if (vf(i,j,k) >= almostone) {
                            for(int n = 0; n < numSurfIntgs; ++n) sintg(i,j,k,n) = 0.0;

                            Real apxm = apx(i  ,j  ,k  );
                            Real apxp = apx(i+1,j  ,k  );
                            Real apym = apy(i  ,j  ,k  );
                            Real apyp = apy(i  ,j+1,k  );
                            Real apzm = apz(i  ,j  ,k  );
                            Real apzp = apz(i  ,j  ,k+1);

                            if (apxm < Real(1.0)) {
                                sintg(i,j,k,i_B_x) = Real(-0.5)*ba(i,j,k);
                            } else if (apym < Real(1.0)) {
                                sintg(i,j,k,i_B_y) = Real(-0.5)*ba(i,j,k);
                            } else if (apzm < Real(1.0)) {
                                sintg(i,j,k,i_B_z) = Real(-0.5)*ba(i,j,k);
                            } else if (apxp < Real(1.0)) {
                                sintg(i,j,k,i_B_x) = Real( 0.5)*ba(i,j,k);
                            } else if (apyp < Real(1.0)) {
                                sintg(i,j,k,i_B_y) = Real( 0.5)*ba(i,j,k);
                            } else if (apzp < Real(1.0)) {
                                sintg(i,j,k,i_B_z) = Real( 0.5)*ba(i,j,k);
                            } else {
                                 amrex::Abort("amrex::algoim::compute_surface_integrals: we are in trouble");
                            }
                        } else {
                            EBPlane phi(bc(i,j,k,0),bc(i,j,k,1),bc(i,j,k,2),
                                        bn(i,j,k,0),bn(i,j,k,1),bn(i,j,k,2));

                            const QuadratureRule q = quadGenSurf(phi);

                            sintg(i,j,k,i_B_x  ) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real /*z*/) noexcept
                                                       { return x; });
                            sintg(i,j,k,i_B_y  ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real /*z*/) noexcept
                                                       { return y; });
                            sintg(i,j,k,i_B_z  ) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real /*y*/, Real z) noexcept
                                                       { return z; });
                            sintg(i,j,k,i_B_x_y) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real /*z*/) noexcept
                                                       { return x*y; });
                            sintg(i,j,k,i_B_x_z) = q([] AMREX_GPU_DEVICE (Real x, Real /*y*/, Real z) noexcept
                                                       { return x*z; });
                            sintg(i,j,k,i_B_y_z) = q([] AMREX_GPU_DEVICE (Real /*x*/, Real y, Real z) noexcept
                                                       { return y*z; });
                            sintg(i,j,k,i_B_xyz) = q([] AMREX_GPU_DEVICE (Real x, Real y, Real z) noexcept
                                                       { return x*y*z; });
                        }
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
                        set_regular_surface(i,j,k,sintg);
                    } else if (ebflag.isCovered()) {
                        for (int n = 0; n < numSurfIntgs; ++n) sintg(i,j,k,n) = 0.0;
                    } else {
                        constexpr Real almostone = Real(1.) - Real(100.)*std::numeric_limits<Real>::epsilon();

                        if (vf(i,j,k) >= almostone) {
                            for(int n = 0; n < numSurfIntgs; ++n) sintg(i,j,k,n) = 0.0;

                            Real apxm = apx(i  ,j  ,k  );
                            Real apxp = apx(i+1,j  ,k  );
                            Real apym = apy(i  ,j  ,k  );
                            Real apyp = apy(i  ,j+1,k  );
                            Real apzm = apz(i  ,j  ,k  );
                            Real apzp = apz(i  ,j  ,k+1);

                            if (apxm < Real(1.0)) {
                                sintg(i,j,k,i_B_x) = Real(-0.5)*ba(i,j,k);
                            } else if (apym < Real(1.0)) {
                                sintg(i,j,k,i_B_y) = Real(-0.5)*ba(i,j,k);
                            } else if (apzm < Real(1.0)) {
                                sintg(i,j,k,i_B_z) = Real(-0.5)*ba(i,j,k);
                            } else if (apxp < Real(1.0)) {
                                sintg(i,j,k,i_B_x) = Real( 0.5)*ba(i,j,k);
                            } else if (apyp < Real(1.0)) {
                                sintg(i,j,k,i_B_y) = Real( 0.5)*ba(i,j,k);
                            } else if (apzp < Real(1.0)) {
                                sintg(i,j,k,i_B_z) = Real( 0.5)*ba(i,j,k);
                            } else {
                                 amrex::Abort("amrex::algoim::compute_surface_integrals: we are in trouble");
                            }
                        } else {
                            EBPlane phi(bc(i,j,k,0),bc(i,j,k,1),bc(i,j,k,2),
                                        bn(i,j,k,0),bn(i,j,k,1),bn(i,j,k,2));

                            const QuadratureRule q = quadGenSurf(phi);

                            sintg(i,j,k,i_B_x  ) = q.eval([](Real x, Real /*y*/, Real /*z*/) noexcept
                                                       { return x; });
                            sintg(i,j,k,i_B_y  ) = q.eval([](Real /*x*/, Real y, Real /*z*/) noexcept
                                                       { return y; });
                            sintg(i,j,k,i_B_z  ) = q.eval([](Real /*x*/, Real /*y*/, Real z) noexcept
                                                       { return z; });
                            sintg(i,j,k,i_B_x_y) = q.eval([](Real x, Real y, Real /*z*/) noexcept
                                                       { return x*y; });
                            sintg(i,j,k,i_B_x_z) = q.eval([](Real x, Real /*y*/, Real z) noexcept
                                                       { return x*z; });
                            sintg(i,j,k,i_B_y_z) = q.eval([](Real /*x*/, Real y, Real z) noexcept
                                                       { return y*z; });
                            sintg(i,j,k,i_B_xyz) = q.eval([](Real x, Real y, Real z) noexcept
                                                       { return x*y*z; });
                        }
                    }
                }
            }
        }
    }
#endif
}

}
