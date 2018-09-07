#include <cmath>

#include <WarpXUtil.H>
#include <WarpXConst.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void ReadBoostedFrameParameters(Real& gamma_boost, Real& beta_boost,
                                Vector<int>& boost_direction)
{
    ParmParse pp("warpx");
    pp.query("gamma_boost", gamma_boost);
    beta_boost = std::sqrt(1.-1./pow(gamma_boost,2));
    if( gamma_boost > 1. ) {
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

    Real gamma_boost, beta_boost;
    int max_level;
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
      pp_wpx.addarr("fine_tag_hi", fine_tag_lo);
    }
}
