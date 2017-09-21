
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>
#include <AMReX_Print.H>

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
#ifdef AMREX_USE_DEVICE
    bool device = false;
#endif
    bool shared = false;

    FabInfo& SetAlloc (bool a) {
        alloc = a;
        return *this;
    }

#ifdef AMREX_USE_DEVICE
    FabInfo& SetDevice (bool d) {
	device = d;
	return *this;
    }
#endif

    FabInfo& SetShared (bool s) {
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
    virtual FabFactory<FAB>* clone () const = 0;
};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
    {
        FAB* fab = new FAB(box, ncomps, info.alloc, info.shared);
#ifdef AMREX_USE_DEVICE
	if (info.device) fab->allocate_device();
#endif
	return fab;
    }

    virtual DefaultFabFactory<FAB>* clone () const override {
        return new DefaultFabFactory<FAB>();
    }
};

}

#endif
