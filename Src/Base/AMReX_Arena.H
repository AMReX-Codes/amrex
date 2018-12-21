
#ifndef BL_ARENA_H
#define BL_ARENA_H

#include <AMReX_BLassert.H>
#include <cstddef>

namespace amrex {

class Arena;

Arena* The_Arena ();
Arena* The_Device_Arena ();
Arena* The_Managed_Arena ();
Arena* The_Pinned_Arena ();

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
    * \brief Given a minimum required arena size of sz bytes, this returns
    * the next largest arena size that will align to align_size bytes
    */
    static std::size_t align (std::size_t sz);

#ifdef AMREX_USE_GPU
    void SetDeviceMemory() { device_use_managed_memory = false; device_use_hostalloc = false; return; }
    void SetReadOnly() { BL_ASSERT(device_use_managed_memory); device_set_readonly = true; return; }
    void SetPreferred() { BL_ASSERT(device_use_managed_memory); device_set_preferred = true; return; }
    void SetHostAlloc() { device_use_hostalloc = true; device_use_managed_memory = false; return; }
#endif

    static void Initialize ();
    static void PrintUsage ();
    static void Finalize ();

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

#ifdef AMREX_USE_GPU
    bool device_use_managed_memory = true;
    bool device_set_readonly = false;
    bool device_set_preferred = false;
    bool device_use_hostalloc = false;
#endif

};

}

#endif /*BL_ARENA_H*/
