#ifndef AMREX_BARENA_H_
#define AMREX_BARENA_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>


namespace amrex {
/**
* \brief A Concrete Class for Dynamic Memory Management
* This is the simplest dynamic memory management class derived from Arena.
* Makes calls to std::malloc and std::free.
*/
class BArena
    :
    public Arena
{
public:
    /**
    * \brief Allocates a dynamic memory arena of size sz.
    * Returns a pointer to this memory.
    */
    [[nodiscard]] void* alloc (std::size_t sz) final;
    /**
    * \brief Deletes the arena pointed to by pt.
    */
    void free (void* pt) final;

    [[nodiscard]] bool isDeviceAccessible () const final;
    [[nodiscard]] bool isHostAccessible () const final;

    [[nodiscard]] bool isManaged () const final;
    [[nodiscard]] bool isDevice () const final;
    [[nodiscard]] bool isPinned () const final;
};

}

#endif
