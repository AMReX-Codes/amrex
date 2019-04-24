
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>
#include <AMReX_Print.H>
#include <AMReX_MakeType.H>
#include <AMReX_Vector.H>

namespace amrex
{

class FArrayBox;

enum class FabType : int {
    covered = -1,
    regular = 0,
    singlevalued = 1,
    multivalued = 2,
    undefined = 100
};

struct FabInfo
{
    bool alloc = true;
    bool shared = false;

    FabInfo& SetAlloc (bool a) noexcept {
        alloc = a;
        return *this;
    }

    FabInfo& SetShared (bool s) noexcept {
        shared = s;
        return *this;
    }
};

template <class FAB>
class FabFactory
{
public:
    virtual ~FabFactory () {}
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
    virtual void destroy (FAB* fab) const = 0;
    virtual FabFactory<FAB>* clone () const = 0;
#ifdef AMREX_USE_GPU
    virtual FAB* createDeviceAlias (FAB const& src) const = 0;
    virtual void destroyDeviceAlias (FAB* fab) const = 0;
#endif
};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
    {
        return ::new FAB(box, ncomps, info.alloc, info.shared);
    }

    virtual void destroy (FAB* fab) const override
    {
        ::delete fab;
    }

#ifdef AMREX_USE_GPU
    virtual FAB* createDeviceAlias (FAB const& src) const override
    {
        return new FAB(src, amrex::make_alias, 0, src.nComp());
    }

    virtual void destroyDeviceAlias (FAB* fab) const override
    {
        if (fab) delete fab;
    }
#endif

    virtual DefaultFabFactory<FAB>* clone () const override {
        return new DefaultFabFactory<FAB>();
    }
};

}

#endif
