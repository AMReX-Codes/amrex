#ifndef AMREX_CUDA_CONTAINERS_H_
#define AMREX_CUDA_CONTAINERS_H_

#include <AMReX_Vector.H>
#include <AMReX_CudaAllocators.H>

namespace amrex {

namespace Cuda {

#ifdef AMREX_USE_CUDA

template <class T>
class DeviceVector
    :
        public thrust::device_vector<T>
{
public:

    using thrust::device_vector<T>::device_vector;
    
    //! get access to the underlying data pointer
    T* dataPtr () { return thrust::raw_pointer_cast(this->data()); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return thrust::raw_pointer_cast(this->data()); }

    void assign(typename thrust::device_vector<T>::size_type count, const T& value )
    {
        thrust::fill(thrust::device, dataPtr(), dataPtr() + count, value);
    }        
};

template <class T>
class HostVector
    :
        public thrust::host_vector<T>
{
public:

    using thrust::host_vector<T>::host_vector;
    
    //! get access to the underlying data pointer
    T* dataPtr () { return thrust::raw_pointer_cast(this->data()); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return thrust::raw_pointer_cast(this->data()); }

    void assign(typename thrust::host_vector<T>::size_type count, const T& value )
    {
        thrust::fill(thrust::host, dataPtr(), dataPtr() + count, value);
    }
};

template <class T>
class ManagedVector
    :
        public Vector<T, CudaManagedAllocator<T> >
{
public:

    using Vector<T, CudaManagedAllocator<T> >::Vector;
    
    void assign(typename thrust::host_vector<T>::size_type count, const T& value )
    {
        if (Cuda::inLaunchRegion())
        {
            thrust::fill(thrust::device,
                         thrust::raw_pointer_cast(this->data()),
                         thrust::raw_pointer_cast(this->data() + count), value);
        }
        else
        {
            this->Vector<T, CudaManagedAllocator<T> >::assign(count, value);
        }
    }
};
        
#else
    template <class T>
    using DeviceVector = Vector<T>;

    template <class T>
    using HostVector = Vector<T>;

    template <class T>
    using ManagedVector = Vector<T>;        
#endif
    
}

    
}
    

#endif
