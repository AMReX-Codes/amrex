#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void
InitAnnulus (MultiFab& a_xyz_loc, Geometry& geom)
{
    const Real tpi = Real(2.0) * amrex::Math::pi<Real>();

#if (AMREX_SPACEDIM == 3)
    const auto dx = geom.CellSizeArray();
#endif

    auto domain = geom.Domain();
    auto probhi = geom.ProbHi();
    auto problo = geom.ProbLo();

    Real ilen = static_cast<Real>(domain.length(0));
    Real jlen = static_cast<Real>(domain.length(1));

    // Center of the annulus
    Real cx = 0.5 * (problo[0]+probhi[0]);
    Real cy = 0.5 * (problo[1]+probhi[1]);

    // Inner radius
    Real ri = 0.1 * (probhi[0]-problo[0]);

    // Outer radius
    Real ro = 0.2 * (probhi[0]-problo[0]);

    // This is "delta r" which is in the j direction
    Real dr = (ro - ri) / jlen;

    // This is "delta theta" which is in the i direction
    Real dtheta = tpi / ilen;

    amrex::Print() << "INNER/OUTER RAD " << ri << " " << ro << std::endl;

    // loc_arr is nodal so no offset
    for(MFIter mfi(a_xyz_loc); mfi.isValid(); ++mfi)
    {
        const Box& gtbx = mfi.growntilebox();

        auto loc_arr = a_xyz_loc.array(mfi);

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real radius = ri + (static_cast<int>(j)) * dr;
            Real theta  =      (static_cast<int>(i)) * dtheta;

            loc_arr(i,j,k,0) = cx + (radius*cos(theta));
            loc_arr(i,j,k,1) = cy + (radius*sin(theta));
#if (AMREX_SPACEDIM == 3)
            loc_arr(i,j,k,2) = problo[2] + k*dx[2];
#endif
            if (i == 0) amrex::Print() << "INITIAL MAPPING AT i=0 " << j << " " << loc_arr(i,j,k,0) << " " << loc_arr(i,j,k,1) << std::endl;
        });

    }
    a_xyz_loc.FillBoundary(geom.periodicity());
}
