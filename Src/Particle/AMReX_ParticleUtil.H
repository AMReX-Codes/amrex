#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>

namespace amrex
{
  namespace utility
  {
    template<class T> 
    AMREX_GPU_HOST_DEVICE
    const T& min(const T& a, const T& b)
    {
      return (b < a) ? b : a;
    }

    template<class T>
    AMREX_GPU_HOST_DEVICE
    const T& max(const T& a, const T& b)
    {
      return (a < b) ? b : a;
    }    

    class IndexConverter
    {
    public:
      Box domain;
      
      IndexConverter(const Box& a_domain) : domain(a_domain) {};

      IntVect To3DIndex(const int in);
      int To1DIndex(const IntVect& in) { return domain.index(in); };
    };
  }

  AMREX_GPU_HOST_DEVICE
  int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
		    const IntVect& a_tile_size, Box& tbx);
  
}

#endif // include guard
