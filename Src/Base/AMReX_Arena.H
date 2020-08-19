#ifndef AMREX_ARENA_H_
#define AMREX_ARENA_H_

#include <AMReX_BLassert.H>
#include <cstddef>
#include <cstdlib>

namespace amrex {

/**
 * \brief Given a minimum required size of size bytes, this returns
 * the next largest arena size that will align to align_requirement bytes
 */
std::size_t aligned_size (std::size_t align_requirement, std::size_t size);

class Arena;

Arena* The_Arena ();
Arena* The_Device_Arena ();
Arena* The_Managed_Arena ();
Arena* The_Pinned_Arena ();
Arena* The_Cpu_Arena ();

struct ArenaInfo
{
    bool use_cpu_memory = false;
    bool device_use_managed_memory = true;
    bool device_set_readonly = false;
    bool device_set_preferred = false;
    bool device_use_hostalloc = false;
    ArenaInfo& SetDeviceMemory () noexcept {
        device_use_managed_memory = false;
        device_use_hostalloc = false;
        return *this;
    }
    ArenaInfo& SetReadOnly () noexcept {
        BL_ASSERT(device_use_managed_memory);
        device_set_readonly = true;
        return *this;
    }
    ArenaInfo& SetPreferred () noexcept {
        BL_ASSERT(device_use_managed_memory);
        device_set_preferred = true;
        return *this;
    }
    ArenaInfo& SetHostAlloc () noexcept {
        device_use_hostalloc = true;
        device_use_managed_memory = false;
        return *this;
    }
    ArenaInfo& SetCpuMemory () noexcept {
        use_cpu_memory = true;
        device_use_managed_memory = false;
        device_set_readonly = false;
        device_set_preferred = false;
        device_use_hostalloc = false;
        return *this; 
    }
};

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

    static void Initialize ();
    static void PrintUsage ();
    static void Finalize ();

#if 0
    union Word
    {
        void*  p;
        long long ll;
        long double ld;
        void (*f) ();
    };
    static const std::size_t align_size = sizeof(Word);
#endif

    static const std::size_t align_size = 16;

protected:

    ArenaInfo arena_info;

    void* allocate_system (std::size_t nbytes);
    void deallocate_system (void* p, std::size_t nbytes);
};

}

#endif /*BL_ARENA_H*/
