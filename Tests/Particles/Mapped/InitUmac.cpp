#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

enum struct ProbType {
    Annulus, Stretched, Hill
};

void
InitUmac (MultiFab* umac, const MultiFab& a_z_loc, Geometry& /*geom*/, ProbType prob_type)
{
    BL_PROFILE("InitUmac");

    // auto probhi = geom.ProbHi();
    // auto problo = geom.ProbLo();

    // For right now we just define a shear flow in x, i.e. in 3D: (u,v,w) = (u(z),0.0,0.0)
    //                                                      in 2D: (u,v)   = (u(y),0.0)

    umac[1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
    umac[2].setVal(0.);
#endif


    for(MFIter mfi(umac[0]); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.growntilebox();
        auto height_arr = a_z_loc.array(mfi);
        auto umac_x_arr = umac[0].array(mfi);

        ParallelFor( tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Physical location of x-face
#if (AMREX_SPACEDIM == 2)
            Real z = Real(0.5)*(height_arr(i,j,k) + height_arr(i,j+1,k));
#elif (AMREX_SPACEDIM == 3)
            Real z = Real(0.25)*(height_arr(i,j,k) + height_arr(i,j+1,k) + height_arr(i,j,k+1) + height_arr(i,j+1,k+1));
#endif

            // Normal velocity on x-face based on height at face center
            umac_x_arr(i,j,k,0) = Real(1.0) + Real(2.0) * z;

#if (AMREX_SPACEDIM == 2)
            if (i == 0) amrex::Print() << "UMAC AT " << IntVect(AMREX_D_DECL(i,j,k)) << " " << umac_x_arr(i,j,k) << std::endl;
#elif (AMREX_SPACEDIM == 3)
            if (i == 0 && j == 0) amrex::Print() << "UMAC AT " << IntVect(AMREX_D_DECL(i,j,k)) << " " << umac_x_arr(i,j,k) << std::endl;
#endif
        });
    }
    int zdir = AMREX_SPACEDIM-1;
    for(MFIter mfi(umac[zdir]); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.growntilebox();
        auto umac_vert_arr = umac[zdir].array(mfi);

        ParallelFor( tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // umac_vert_arr(i,j,k) = 0.05;
            umac_vert_arr(i,j,k) = 0.0;
        });
    }

}
