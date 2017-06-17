
#ifndef BL_ARENA_H
#define BL_ARENA_H

#include <cstddef>

namespace amrex {

class Arena;

Arena* The_Arena ();

/**
* \brief 
* A virtual base class for objects that manage their own dynamic
* memory allocation.  
*/

class Arena
{
public:

    virtual ~Arena ();
    /**
    * Allocate a dynamic memory arena
    * \param sz size of the memory request
    * \return a pointer to the allocated memory
    */
    virtual void* alloc (std::size_t sz) = 0;
    /**
    * \brief A pure virtual function for deleting the arena pointed to by pt
    */
    virtual void free (void* pt) = 0;
    /**
    * Allocate a dynamic device memory arena
    * \param sz size of the memory request
    * \return a pointer to the allocated memory
    */
    virtual void* alloc_device (std::size_t sz) = 0;
    /**
    * \brief A pure virtual function for deleting the device arena pointed to by pt
    */
    virtual void free_device (void* pt) = 0;
    /**
    * \brief Given a minimum required arena size of sz bytes, this returns
    * the next largest arena size that will align to align_size bytes
    */
    static std::size_t align (std::size_t sz);

#ifdef CUDA
    void SetDeviceMemory() { device_use_managed_memory = false; }
    void SetReadOnly() { device_set_readonly = true; }
    void SetPreferred() { device_set_preferred = true; }
#endif

protected:

#if 0
    union Word
    {
        void*  p;
        double d;
        long   l;
        void (*f) ();
    };
    static const unsigned int align_size = sizeof(Word);
#endif

    static const unsigned int align_size = 16;

#ifdef CUDA
    bool device_use_managed_memory = true;
    bool device_set_readonly = false;
    bool device_set_preferred = false;
#endif

};

}

#endif /*BL_ARENA_H*/
