#include <cmath>

#include <WarpXUtil.H>
#include <WarpXConst.H>
#include <AMReX_ParmParse.H>
#include <WarpX.H>

using namespace amrex;

void ReadBoostedFrameParameters(Real& gamma_boost, Real& beta_boost,
                                Vector<int>& boost_direction)
{
    ParmParse pp("warpx");
    pp.query("gamma_boost", gamma_boost);
    if( gamma_boost > 1. ) {
        beta_boost = std::sqrt(1.-1./pow(gamma_boost,2));
        std::string s;
        pp.get("boost_direction", s);
        if (s == "x" || s == "X") {
            boost_direction[0] = 1;
        }
#if (AMREX_SPACEDIM == 3)
        else if (s == "y" || s == "Y") {
            boost_direction[1] = 1;
        }
#endif
        else if (s == "z" || s == "Z") {
            boost_direction[2] = 1;
        }
        else {
            const std::string msg = "Unknown boost_dir: "+s;
            Abort(msg.c_str());
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( s == "z" || s == "Z" ,
                                          "The boost must be in the z direction.");
    }
}

void ConvertLabParamsToBoost()
{
    Real gamma_boost = 1., beta_boost = 0.;
    int max_level = 0;
    Vector<int> boost_direction {0,0,0};

    ReadBoostedFrameParameters(gamma_boost, beta_boost, boost_direction);

    if (gamma_boost <= 1.) return;

    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);
    Vector<Real> fine_tag_lo(AMREX_SPACEDIM);
    Vector<Real> fine_tag_hi(AMREX_SPACEDIM);

    ParmParse pp_geom("geometry");
    ParmParse pp_wpx("warpx");
    ParmParse pp_amr("amr");

    pp_geom.getarr("prob_lo",prob_lo,0,AMREX_SPACEDIM);
    BL_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
    pp_geom.getarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);
    BL_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_amr.query("max_level", max_level);
    if (max_level > 0){
      pp_wpx.getarr("fine_tag_lo", fine_tag_lo);
      pp_wpx.getarr("fine_tag_hi", fine_tag_hi);
    }


#if (AMREX_SPACEDIM == 3)
    Vector<int> dim_map {0, 1, 2};
#else
    Vector<int> dim_map {0, 2};
#endif

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (boost_direction[dim_map[idim]]) {
            amrex::Real convert_factor;
            // Assume that the window travels with speed +c
            convert_factor = 1./( gamma_boost * ( 1 - beta_boost ) );
            prob_lo[idim] *= convert_factor;
            prob_hi[idim] *= convert_factor;
            if (max_level > 0){
              fine_tag_lo[idim] *= convert_factor;
              fine_tag_hi[idim] *= convert_factor;
            }
            break;
        }
    }
    pp_geom.addarr("prob_lo", prob_lo);
    pp_geom.addarr("prob_hi", prob_hi);
    if (max_level > 0){
      pp_wpx.addarr("fine_tag_lo", fine_tag_lo);
      pp_wpx.addarr("fine_tag_hi", fine_tag_hi);
    }
}

/* \brief Function that sets the value of MultiFab MF to zero for z between 
 * zmin and zmax.
 */
void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin, amrex::Real zmax){
    BL_PROFILE("WarpX::NullifyMF()");
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();
        // Get box lower and upper physical z bound, and dz
        const amrex::Real zmin_box = WarpX::LowerCorner(bx, lev)[2];
        const amrex::Real zmax_box = WarpX::UpperCorner(bx, lev)[2];
        amrex::Real dz  = WarpX::CellSize(lev)[2];
        // Get box lower index in the z direction
#if (AMREX_SPACEDIM==3)
        const int lo_ind = bx.loVect()[2];
#else
        const int lo_ind = bx.loVect()[1];
#endif
        // Check if box intersect with [zmin, zmax]
        if ( (zmax>zmin_box && zmin<=zmax_box) ){
            Array4<Real> arr = mf[mfi].array();
            // Set field to 0 between zmin and zmax
            ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept{
#if (AMREX_SPACEDIM==3)
                    const Real z_gridpoint = zmin_box+(k-lo_ind)*dz;
#else
                    const Real z_gridpoint = zmin_box+(j-lo_ind)*dz;
#endif 
                    if ( (z_gridpoint >= zmin) && (z_gridpoint < zmax) ) {
                        arr(i,j,k) = 0.;
                    };
                }
            );
        }
    }
}
