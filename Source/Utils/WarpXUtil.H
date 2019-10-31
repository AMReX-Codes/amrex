#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <string>

void ReadBoostedFrameParameters(amrex::Real& gamma_boost, amrex::Real& beta_boost,
                                amrex::Vector<int>& boost_direction);

void ConvertLabParamsToBoost();

void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin,
               amrex::Real zmax);

namespace WarpXUtilIO{
    /**
     * A helper function to write binary data on disk.
     * @param[in] filename where to write
     * @param[in] data Vector containing binary data to write on disk
     * return true if it succeeds, false otherwise
     */
    bool WriteBinaryDataOnFile(std::string filename, const amrex::Vector<char>& data);
}
