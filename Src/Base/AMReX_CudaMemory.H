#ifndef AMREX_CUDA_MEMORY_H_
#define AMREX_CUDA_MEMORY_H_

#include <AMReX_CudaQualifiers.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_CudaControl.H>
#include <cstdlib>

namespace amrex {

struct Managed {

#ifdef AMREX_USE_CUDA

  void *operator new(size_t len)
  {
    void *ptr;
    cudaMallocManaged(&ptr, len);
    cudaDeviceSynchronize();
    return ptr;
  }

  void operator delete(void *ptr)
  {
    cudaDeviceSynchronize();
    cudaFree(ptr);
  }

#endif
};

struct Pinned {

#ifdef AMREX_USE_CUDA

  void *operator new(size_t len)
  {
    void *ptr;
    cudaMallocHost(&ptr, len);
    cudaDeviceSynchronize();
    return ptr;
  }

  void operator delete(void *ptr)
  {
    cudaDeviceSynchronize();
    cudaFreeHost(ptr);
  }

#endif
};

template <class T, typename std::enable_if<std::is_pod<T>::value,int>::type = 0 > 
struct ManagedData
{
    ManagedData ()
    {
        cudaMalloc(&d_d, std::size_t(sizeof(T)));
    }

    ManagedData (T const & h_d)
    : ManagedData()
    {
        *d_d = h_d;
    }

    ~ManagedData ()
    {
        cudaFree(d_d);
    }

    T* devicePtr() &
    {
        return d_d;
    }

    T const * devicePtr() const&
    {
        return d_d;
    }

    T hostValue () const
    {
        T t;
//        cudaMemcpy(&t, d_d, sizeof(T), cudaMemcpyDeviceToHost);
        return t; 
    }

    void updateDevice(const T &t)
    {
//        cudaMemcpy(d_d, &t, sizeof(T), cudaMemcpyHostToDevice);
    }

    T* data() && = delete; 
    ManagedData(ManagedData const &) = delete;
    ManagedData(ManagedData &&) = delete;
    void operator = (ManagedData const &) = delete;
    void operator = (ManagedData &&) = delete; 

    private:
    T* d_d = nullptr;
};


template <class T, typename = amrex::EnableIf_t<std::is_pod<T>::value> >
struct DeviceScalar
{
    DeviceScalar (DeviceScalar const&) = delete;
    DeviceScalar (DeviceScalar &&) = delete;
    void operator= (DeviceScalar const&) = delete;
    void operator= (DeviceScalar &&) = delete;

#if AMREX_USE_CUDA

    DeviceScalar (T init_val) {
        if (Cuda::inLaunchRegion()) {
            cudaMalloc(&dp, sizeof(T));
            cudaMemcpy(dp, &init_val, sizeof(T), cudaMemcpyHostToDevice);
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
            *dp = init_val;
        }
    }

    ~DeviceScalar () {
        if (Cuda::inLaunchRegion()) {
            cudaFree(dp);
        } else {
            std::free(dp);
        }
    }

    T* dataPtr () { return dp; }
    T const* dataPtr () const { return dp; }
    T dataValue () const {
        if (Cuda::inLaunchRegion()) {
            T r;
            cudaMemcpy(&r, dp, sizeof(T), cudaMemcpyDeviceToHost);
            return r;
        } else {
            return *dp;
        }
    }

private:
    T* dp;

#else

    DeviceScalar (T init_val) : d(init_val) {}
    ~DeviceScalar () {}

    T* dataPtr () { return &d; }
    T const* dataPtr () const { return &d; }
    T dataValue () const { return d; }

private:
    T d;

#endif
};

}  //namespace


#endif
