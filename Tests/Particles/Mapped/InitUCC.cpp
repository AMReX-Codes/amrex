#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

enum struct ProbType {
    Annulus, Stretched, Hill
};

void
InitUCC (MultiFab& ucc, const MultiFab& a_xyz_loc, Geometry& geom, ProbType& prob_type)
{
    BL_PROFILE("InitUCC");

    auto probhi = geom.ProbHi();
    auto problo = geom.ProbLo();

    // Center of the annulus
    Real cx = 0.5 * (problo[0]+probhi[0]);
    Real cy = 0.5 * (problo[1]+probhi[1]);

    amrex::Print() << "UCC NG " << ucc.nGrow() << std::endl;

    if (prob_type == ProbType::Annulus) {

        AMREX_ALWAYS_ASSERT(a_xyz_loc.nComp() == AMREX_SPACEDIM);

        for (MFIter mfi(ucc); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();

            auto loc_arr = a_xyz_loc.const_array(mfi);
            auto ucc_arr = ucc.array(mfi);

            ParallelFor( tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Physical location of cell center
                Real x = 0.125 * (loc_arr(i,j  ,k  ,0) + loc_arr(i+1,j  ,k  ,0) +
                                  loc_arr(i,j+1,k  ,0) + loc_arr(i+1,j+1,k  ,0) +
                                  loc_arr(i,j  ,k+1,0) + loc_arr(i+1,j  ,k+1,0) +
                                  loc_arr(i,j+1,k+1,0) + loc_arr(i+1,j+1,k+1,0) );
                Real y = 0.125 * (loc_arr(i,j  ,k  ,1) + loc_arr(i+1,j  ,k  ,1) +
                                  loc_arr(i,j+1,k  ,1) + loc_arr(i+1,j+1,k  ,1) +
                                  loc_arr(i,j  ,k+1,1) + loc_arr(i+1,j  ,k+1,1) +
                                  loc_arr(i,j+1,k+1,1) + loc_arr(i+1,j+1,k+1,1) );

                Real theta;
                if (x == cx) {
                   theta = Real(0.);
                } else {
                   theta = atan((y-cy)/(x-cx));
                }

                Real    rad = sqrt( x*x + y*y);

                ucc_arr(i,j,k,0) =  rad*sin(theta);
                ucc_arr(i,j,k,1) = -rad*cos(theta);

#if (AMREX_SPACEDIM == 3)
                Real z = loc_arr(i,j,k,2);
                ucc_arr(i,j,k,2) =  0.0;
#endif

                if (i == 0) amrex::Print() << "UCC AT " <<  IntVect(AMREX_D_DECL(i,j,k)) << " "
                                                        << RealVect(AMREX_D_DECL(ucc_arr(i,j,k,0),ucc_arr(i,j,k,1),ucc_arr(i,j,k,2)))
                                                        << std::endl;


            });
        }

    } else if (prob_type == ProbType::Stretched) {

        // Set to 1 in x-direction
        ucc.setVal(1.0,0,1,ucc.nGrow());

        // Set to 0 in other directions
        ucc.setVal(0.0,1,AMREX_SPACEDIM-1,ucc.nGrow());

    } else if (prob_type == ProbType::Hill) {

        // Set to 1 in x-direction
        ucc.setVal(1.0,0,1,ucc.nGrow());

        // Set to 0 in other directions
        ucc.setVal(0.0,1,AMREX_SPACEDIM-1,ucc.nGrow());
    }

}
