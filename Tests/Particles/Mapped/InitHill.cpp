#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void
InitHill (MultiFab& a_z_loc, Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(a_z_loc.nComp() == 1);

    const Real hpi = Real(0.5) * amrex::Math::pi<Real>();
    const Real tpi = Real(2.0) * amrex::Math::pi<Real>();

    const int zdir = AMREX_SPACEDIM-1;

    auto domain = geom.Domain();
    auto problo = geom.ProbLo();
    auto probhi = geom.ProbHi();
    const auto dx = geom.CellSizeArray();

    // Define dz_unit so that at k = khi, k*dz_unit = 1
    Real  z_size = (probhi[zdir] - problo[zdir]);
    Real dz_unit = z_size / static_cast<Real>(domain.length()[zdir]);

    for(MFIter mfi(a_z_loc); mfi.isValid(); ++mfi)
    {
        const Box& gtbx = mfi.growntilebox();

        auto loc_arr = a_z_loc.array(mfi);

#if (AMREX_SPACEDIM == 2)
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            Real x = static_cast<Real>(i) * dx[0];
            Real floor = 0.25 * (1.0 + cos(tpi * x));
            loc_arr(i,j,0,0) = floor + 0.5 * (z_size - floor) * sin(hpi * j * dz_unit);
            if (i == 0) amrex::Print() << "INITIAL MAPPING AT i=0 " << j << " " << loc_arr(i,j,0,1) << std::endl;
        });
#elif (AMREX_SPACEDIM == 3)
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x     = static_cast<Real>(i) * dx[0];
            Real floor = cos(tpi * x);
            loc_arr(i,j,k,0) = floor + (z_size - floor) * sin(hpi * j * dz_unit);
            if (i == 0) amrex::Print() << "INITIAL MAPPING AT i=0 " << k << " " << loc_arr(i,j,k,0) << " " << loc_arr(i,j,k,1) << std::endl;
        });
#endif
    }
    a_z_loc.FillBoundary(geom.periodicity());
}
