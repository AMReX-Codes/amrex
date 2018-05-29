
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
    * \brief Version for managed CUDA memory. 
    * Allocates a dynamic memory arena of size sz.
    * Memory is allocated based on attributes of parent object.
    * Returns a pointer to this memory.
    */
    virtual void* alloc (void* parent_ptr, std::size_t sz) override;
    /**
    * \brief Deletes the arena pointed to by pt.
    */
    virtual void free (void* pt) override;
};

}

#endif /*BL_BARENA_H*/
