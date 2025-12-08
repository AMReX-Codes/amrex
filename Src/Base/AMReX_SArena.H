#ifndef AMREX_SARENA_H_
#define AMREX_SARENA_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>


namespace amrex {
/**
* \brief Doc TODO
*/
class SArena
    :
    public Arena
{
public:
    /**
    * \brief Allocates a dynamic memory arena of size sz.
    * Returns a pointer to this memory.
    * Doc TODO
    */
    [[nodiscard]] void* alloc (std::size_t sz) final;
    /**
    * \brief Deletes the arena pointed to by pt.
    * Doc TODO
    */
    void free (void* pt) final;

    [[nodiscard]] bool isDeviceAccessible () const final;
    [[nodiscard]] bool isHostAccessible () const final;

    [[nodiscard]] bool isManaged () const final;
    [[nodiscard]] bool isDevice () const final;
    [[nodiscard]] bool isPinned () const final;

#ifdef AMREX_USE_GPU
    //! Is this CUDA stream ordered memory allocator?
    [[nodiscard]] bool isStreamOrderedArena () const final { return true; }
#endif
};

}

#endif
