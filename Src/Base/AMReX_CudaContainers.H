#ifndef AMREX_CUDA_CONTAINERS_H_
#define AMREX_CUDA_CONTAINERS_H_

#include <AMReX_Vector.H>
#include <AMReX_PODVector.H>
#include <AMReX_CudaAllocators.H>

namespace amrex {

namespace Cuda {

#ifdef AMREX_USE_CUDA

    template <class T>
    using DeviceVector = PODVector<T, ThrustDeviceAllocator<T> >;
	
    template <class T>
    using HostVector = PODVector<T>;
    
    template <class T>
    using ManagedDeviceVector = PODVector<T, ManagedAllocator<T> >;
    
    template <class T>
    using ManagedVector = PODVector<T, ManagedAllocator<T> >;
	
/**

   This is thrust::device_vector allocated using pinned memory.

 **/
template <class T>
class PinnedDeviceVector
    :
    public thrust::device_vector<T, ThrustPinnedAllocator<T> >
{
public:

    using thrust::device_vector<T, ThrustPinnedAllocator<T> >::device_vector;

    //! get access to the underlying data pointer
    T* data ()
    { 
      return thrust::raw_pointer_cast(
       this->thrust::device_vector<T, ThrustPinnedAllocator<T>>::data());
    }

    //! get access to the underlying data pointer
    const T* data () const
    { 
      return thrust::raw_pointer_cast(
       this->thrust::device_vector<T, ThrustPinnedAllocator<T>>::data());
    }

    //! get access to the underlying data pointer
    T* dataPtr () { return this->data(); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return this->data(); }

    T& operator[] (std::size_t i)
    {
      return *thrust::raw_pointer_cast(
      &(this->thrust::device_vector<T, ThrustPinnedAllocator<T>>::operator[](i)));
    }

    const T& operator[] (std::size_t i) const
    {
      return *thrust::raw_pointer_cast(
      &(this->thrust::device_vector<T, ThrustPinnedAllocator<T>>::operator[](i)));
    }
    
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
	    this->thrust::device_vector<T, ThrustPinnedAllocator<T> >::assign(count, value);
        }
    }
};

/**

   This is thrust::device_vector allocated using either pinned or device memory, based on a
   runtime switch.

 **/
template <class T>
class PolymorphicDeviceVector
    :
    public thrust::device_vector<T, ThrustPolymorphicAllocator<T> >
{
public:

    using thrust::device_vector<T, ThrustPolymorphicAllocator<T> >::device_vector;

    //! get access to the underlying data pointer
    T* data ()
    { 
      return thrust::raw_pointer_cast(
       this->thrust::device_vector<T, ThrustPolymorphicAllocator<T>>::data());
    }

    //! get access to the underlying data pointer
    const T* data () const
    { 
		return thrust::raw_pointer_cast(
		    this->thrust::device_vector<T, ThrustPolymorphicAllocator<T>>::data());
    }

    //! get access to the underlying data pointer
    T* dataPtr () { return this->data(); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return this->data(); }

    T& operator[] (std::size_t i)
    {
      return *thrust::raw_pointer_cast(
      &(this->thrust::device_vector<T, ThrustPolymorphicAllocator<T>>::operator[](i)));
    }

    const T& operator[] (std::size_t i) const
    {
      return *thrust::raw_pointer_cast(
      &(this->thrust::device_vector<T, ThrustPolymorphicAllocator<T>>::operator[](i)));
    }
    
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
	    this->thrust::device_vector<T, ThrustPolymorphicAllocator<T> >::assign(count, value);
        }
    }
};

#else
    //! When Cuda is off, all these containers revert to amrex::Vector. 
    template <class T>
    using DeviceVector = PODVector<T>;

    template <class T>
    using HostVector = PODVector<T>;

    template <class T>
    using ManagedVector = PODVector<T>;

    template <class T>
    using ManagedDeviceVector = PODVector<T>;

    template <class T>
    using PinnedDeviceVector = PODVector<T>;

    template <class T>
    using PolymorphicDeviceVector = PODVector<T>;
#endif

    template<class InIter, class OutIter>
    OutIter thrust_copy (InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_CUDA
        return thrust::copy(begin, end, result);
#else
        return std::copy(begin, end, result);
#endif
    }

}


}


#endif
