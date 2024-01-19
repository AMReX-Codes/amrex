#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

enum struct ProbType {
    Annulus, Stretched, Hill
};

void
InitUND (MultiFab& und, const MultiFab& a_xyz_loc, Geometry& geom, ProbType& prob_type)
{
    BL_PROFILE("InitUND");

    auto probhi = geom.ProbHi();
    auto problo = geom.ProbLo();

    // Center of the annulus
    Real cx = 0.5 * (problo[0]+probhi[0]);
    Real cy = 0.5 * (problo[1]+probhi[1]);

    amrex::Print() << "UND NG " << und.nGrow() << std::endl;

    if (prob_type == ProbType::Annulus) {

        // We only do this problem with fully mapped coordinates
        AMREX_ALWAYS_ASSERT(a_xyz_loc.nComp() == AMREX_SPACEDIM);

        for (MFIter mfi(und); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();

            auto loc_arr = a_xyz_loc.const_array(mfi);
            auto und_arr = und.array(mfi);

            ParallelFor( tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Physical location of cell center
                Real x = loc_arr(i,j,k,1);
                Real y = loc_arr(i,j,k,1);

                Real theta;
                if (x == cx) {
                   theta = Real(0.);
                } else {
                   theta = atan((y-cy)/(x-cx));
                }

                Real    rad = sqrt( x*x + y*y);

                und_arr(i,j,k,0) =  rad*sin(theta);
                und_arr(i,j,k,1) = -rad*cos(theta);

#if (AMREX_SPACEDIM == 3)
                Real z = loc_arr(i,j,k,2);
                und_arr(i,j,k,2) =  0.0;
#endif

                if (i == 0) amrex::Print() << "UND AT " <<  IntVect(AMREX_D_DECL(i,j,k)) << " "
                                                        << RealVect(AMREX_D_DECL(und_arr(i,j,k,0),und_arr(i,j,k,1),und_arr(i,j,k,2)))
                                                        << std::endl;


            });
        }

    } else if (prob_type == ProbType::Stretched) {

        // Set to 1 in x-direction
        und.setVal(1.0,0,1,und.nGrow());

        // Set to 0 in other directions
        und.setVal(0.0,1,AMREX_SPACEDIM-1,und.nGrow());

    } else if (prob_type == ProbType::Hill) {

        // Set to 1 in x-direction
        und.setVal(1.0,0,1,und.nGrow());

        // Set to 0 in other directions
        und.setVal(0.0,1,AMREX_SPACEDIM-1,und.nGrow());
    }

}
