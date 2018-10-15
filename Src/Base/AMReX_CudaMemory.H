#ifndef AMREX_CUDA_MEMORY_H_
#define AMREX_CUDA_MEMORY_H_

#include <AMReX_CudaQualifiers.H>

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

template <class T, typename std::enable_if<std::is_pod<T>::value,int>::type = 0 >
struct HostDeviceScalar
{
    HostDeviceScalar (T& h) : hp(&h) {
#if AMREX_USE_CUDA
        cudaMalloc(&dp, sizeof(T));
        cudaMemcpy(dp, hp, sizeof(T), cudaMemcpyHostToDevice);
#else
        dp = hp;
#endif
    }
    
    ~HostDeviceScalar () {
#if AMREX_USE_CUDA
        cudaMemcpy(hp, dp, sizeof(T), cudaMemcpyDeviceToHost);
        cudaFree(dp);
#endif
    }

    HostDeviceScalar (HostDeviceScalar const&) = delete;
    HostDeviceScalar (HostDeviceScalar &&) = delete;
    void operator= (HostDeviceScalar const&) = delete;
    void operator= (HostDeviceScalar &&) = delete;

    T* devicePtr () { return dp; }
    T const* devicePtr () const { return dp; }

    void DeviceToHost () {
#if AMREX_USE_CUDA
        cudaMemcpy(hp, dp, sizeof(T), cudaMemcpyDeviceToHost);
#endif
    }

    void HostToDevice () {
#if AMREX_USE_CUDA
        cudaMemcpy(dp, hp, sizeof(T), cudaMemcpyHostToDevice);
#endif
    }

private:
    T* hp;
    T* dp;
};

}  //namespace


#endif
