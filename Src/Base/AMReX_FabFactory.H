
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>
#include <AMReX_Print.H>
#include <AMReX_MakeType.H>
#include <AMReX_Vector.H>
#include <AMReX_Arena.H>

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
    Arena* arena = nullptr;

    FabInfo& SetAlloc (bool a) noexcept {
        alloc = a;
        return *this;
    }

    FabInfo& SetShared (bool s) noexcept {
        shared = s;
        return *this;
    }

    FabInfo& SetArena (Arena* ar) noexcept {
        arena = ar;
        return *this;
    }
};

template <class FAB>
class FabFactory
{
public:
    virtual ~FabFactory () {}
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
    virtual FAB* create_alias (FAB const& /*rhs*/, int /*scomp*/, int /*ncomp*/) const { return nullptr; }
    virtual void destroy (FAB* fab) const = 0;
    virtual FabFactory<FAB>* clone () const = 0;
};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int /*box_index*/) const override
    {
        return new FAB(box, ncomps, info.alloc, info.shared, info.arena);
    }

    virtual FAB* create_alias (FAB const& rhs, int scomp, int ncomp) const override
    {
        return new FAB(rhs, amrex::make_alias, scomp, ncomp);
    }

    virtual void destroy (FAB* fab) const override
    {
        delete fab;
    }

    virtual DefaultFabFactory<FAB>* clone () const override {
        return new DefaultFabFactory<FAB>();
    }
};

}

#endif
