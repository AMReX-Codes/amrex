#ifndef AMREX_FFT_HELPER_H_
#define AMREX_FFT_HELPER_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Math.H>

namespace amrex::FFT
{

enum struct Scaling { full, symmetric, none };
enum struct Direction { forward, backward };

struct Info
{
    //! Supported only in 3D. When batch_mode is true, FFT is performed on
    //! the first two dimensions only and the third dimension size is the
    //! batch size.
    bool batch_mode = false;

    Info& setBatchMode (bool x) { batch_mode = x; return *this; }
};

template <typename T>
struct PoissonSpectral
{
    PoissonSpectral (Geometry const& geom)
        : fac({AMREX_D_DECL(T(2)*Math::pi<T>()/T(geom.ProbLength(0)),
                            T(2)*Math::pi<T>()/T(geom.ProbLength(1)),
                            T(2)*Math::pi<T>()/T(geom.ProbLength(2)))}),
          dx({AMREX_D_DECL(T(geom.CellSize(0)),
                           T(geom.CellSize(1)),
                           T(geom.CellSize(2)))}),
          scale(T(1.0/geom.Domain().d_numPts())),
          len(geom.Domain().length())
    {
        static_assert(std::is_floating_point_v<T>);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (int i, int j, int k, GpuComplex<T>& spectral_data) const
    {
        amrex::ignore_unused(i,j,k);
        // the values in the upper-half of the spectral array in y and z
        // are here interpreted as negative wavenumbers
        AMREX_D_TERM(T a = fac[0]*i;,
                     T b = (j < len[1]/2) ? fac[1]*j : fac[1]*(len[1]-j);,
                     T c = (k < len[2]/2) ? fac[2]*k : fac[2]*(len[2]-k));
        T k2 = AMREX_D_TERM(T(2)*(std::cos(a*dx[0])-T(1))/(dx[0]*dx[0]),
                           +T(2)*(std::cos(b*dx[1])-T(1))/(dx[1]*dx[1]),
                           +T(2)*(std::cos(c*dx[2])-T(1))/(dx[2]*dx[2]));
        if (k2 != T(0)) {
            spectral_data /= k2;
        } else {
            // interpretation here is that the average value of the
            // solution is zero
            spectral_data = 0;
        }
        spectral_data *= scale;
    };

    GpuArray<T,AMREX_SPACEDIM> fac;
    GpuArray<T,AMREX_SPACEDIM> dx;
    T scale;
    IntVect len;
};

}

#endif
