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

namespace WarpXUtilAlgo{

template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
const T* upper_bound(const T* first, const T* last, const T& val)
{
    const T* it;
    size_t count, step;
    count = last-first;
    while(count>0){
        it = first;
        step = count/2;
        it += step;
        if (!(val<*it)){
            first = ++it;
                count -= step + 1;
        }
        else{
            count = step;
        }
    }
        return first;
}

template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T bilinear_interp(T x0, T x1, T y0, T y1, T f00, T f01, T f10, T f11, T x, T y)
{
    const T one_over_dx = static_cast<T>(1.0)/(x1-x0);
    const T fx0 = ((x1-x)*f00 + (x-x0)*f10)*one_over_dx;
    const T fx1 = ((x1-x)*f01 + (x-x0)*f11)*one_over_dx;
    return ((y1-y)*fx0 + (y-y0)*fx1)/(y1-y0);
}

template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T trilinear_interp(T x0, T x1,T y0, T y1, T z0, T z1,
    T f000, T f001, T f010, T f011, T f100, T f101, T f110, T f111,
    T x, T y, T z)
{
    const T fxy0 = bilinear_interp(
        x0, x1, y0, y1, f000, f010, f100, f110, x, y);
    const T fxy1 = bilinear_interp(
        x0, x1, y0, y1, f001, f011, f101, f111, x, y);
    );
    return ((z1-z)*fxy0 + (z-z0)*fxy1)/(z1-z0);
}

}