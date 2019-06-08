#ifndef AMREX_BARENA_H_
#define AMREX_BARENA_H_

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
    virtual void* alloc (std::size_t sz) override;
    /**
    * \brief Deletes the arena pointed to by pt.
    */
    virtual void free (void* pt) override;
};

}

#endif
