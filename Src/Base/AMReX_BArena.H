
#ifndef BL_BARENA_H
#define BL_BARENA_H

#include <AMReX_Arena.H>


namespace amrex {
/**
* \brief A Concrete Class for Dynamic Memory Management
* This is the simplest dynamic memory management class derived from Arena.
* Makes calls to ::operator new() and ::operator delete().
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
    virtual void* alloc (std::size_t sz) override;
    /**
    * \brief Deletes the arena pointed to by pt.
    */
    virtual void free (void* pt) override;
#ifdef AMREX_USE_CUDA
    /**
    * \brief Allocates a dynamic device memory arena of size sz.
    * Returns a pointer to this memory.
    */
    virtual void* alloc_device (std::size_t sz) override;
    /**
    * \brief Deletes the device arena pointed to by pt.
    */
    virtual void free_device (void* pt) override;
#endif
};

}

#endif /*BL_BARENA_H*/
