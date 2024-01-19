#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void
InitStretched (MultiFab& a_z_loc, Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(a_z_loc.nComp() == 1);

    bool verbose = true;

    const Real hpi = Real(0.5) * amrex::Math::pi<Real>();

    const int zdir = AMREX_SPACEDIM-1;

    auto domain = geom.Domain();
    auto problo = geom.ProbLo();
    auto probhi = geom.ProbHi();
    // const auto dx = geom.CellSizeArray();

    // Define dz_unit so that at k = khi, k*dz_unit = 1
    Real  z_size = (probhi[zdir] - problo[zdir]);
    Real dz_unit = z_size / static_cast<Real>(domain.length()[zdir]);

    for (MFIter mfi(a_z_loc); mfi.isValid(); ++mfi)
    {
        const Box&  tbx = mfi.tilebox();
        const Box& gtbx = mfi.grownnodaltilebox(1);
        // amrex::Print() << "TBX  " <<  tbx << std::endl;
        // amrex::Print() << "GTBX " << gtbx << std::endl;

        auto tlo = lbound(tbx);
        auto thi = ubound(tbx);

        auto loc_arr = a_z_loc.array(mfi);

#if (AMREX_SPACEDIM == 2)
        ParallelFor(makeSlab(gtbx,1,0), [=] AMREX_GPU_DEVICE (int i, int , int ) noexcept
        {
            for (int j = tlo.y; j <= thi.y; j++)
            {
                loc_arr(i,j,0,0) = problo[1] + z_size * sin(hpi * j * dz_unit);
            }
            loc_arr(i,tlo.y-1,0,0) = 2.0 * loc_arr(i,tlo.y,0,0) -  loc_arr(i,tlo.y+1,0,0);
            loc_arr(i,thi.y+1,0,0) = 2.0 * loc_arr(i,thi.y,0,0) -  loc_arr(i,thi.y-1,0,0);

            if (verbose && i == 0) {
                for (int j = tlo.y-1; j <= thi.y+1; j++) {
                    if (j < 0) {
                        amrex::Print() << "INITIAL MAPPING AT " << IntVect(i,j) << " " << loc_arr(i,j,0,0) <<  std::endl;
                    } else {
                        amrex::Print() << "INITIAL MAPPING AT " << IntVect(i,j) << " " << loc_arr(i,j,0,0) <<
                                          " with dz = " << loc_arr(i,j,0,0) - loc_arr(i,j-1,0,0) << std::endl;
                    }
                } // j
            } // i
        });
#elif (AMREX_SPACEDIM == 3)
        ParallelFor(makeSlab(gtbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
        {
            for (int k = tlo.z; k <= thi.z; k++)
            {
                loc_arr(i,j,k,0) = problo[2] + z_size * sin(hpi * k * dz_unit);
            }
            loc_arr(i,j,tlo.z-1,0) = 2.0 * loc_arr(i,j,tlo.z,0) -  loc_arr(i,j,tlo.z+1,0);
            loc_arr(i,j,thi.z+1,0) = 2.0 * loc_arr(i,j,thi.z,0) -  loc_arr(i,j,thi.z-1,0);

            if (verbose && i == 0 && j == 0) {
                for (int k = tlo.z-1; k <= thi.z+1; k++) {
                    if (k < 0) {
                        amrex::Print() << "INITIAL MAPPING AT " << IntVect(i,j,k) << " " << loc_arr(i,j,k,0) <<  std::endl;
                    } else {
                        amrex::Print() << "INITIAL MAPPING AT " << IntVect(i,j,k) << " " << loc_arr(i,j,k,0) <<
                                          " with dz = " << loc_arr(i,j,k,0) - loc_arr(i,j,k-1,0) << std::endl;
                    }
                } // j
            } // i
        });
#endif
    }
    a_z_loc.FillBoundary(geom.periodicity());
}
