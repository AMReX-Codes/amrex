
#ifndef BL_CUDARENA_H
#define BL_CUDARENA_H

#include <AMReX_Arena.H>


namespace amrex {
/**
* \brief A Concrete Class for Dynamic Memory Management
* This class allocates memory using CUDA managed memory.
*/

class CUDArena
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

#endif /*BL_CUDARENA_H*/
