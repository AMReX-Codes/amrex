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
    virtual void* alloc (std::size_t sz) override final;
    /**
    * \brief Deletes the arena pointed to by pt.
    */
    virtual void free (void* pt) override final;

    virtual bool isDeviceAccessible () const override final;
    virtual bool isHostAccessible () const override final;

    virtual bool isManaged () const override final;
    virtual bool isDevice () const override final;
    virtual bool isPinned () const override final;
};

}

#endif
