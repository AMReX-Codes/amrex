
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>

namespace amrex
{

struct FabInfo
{
    bool alloc = true;
    bool device = false;
    bool shared = false;

    FabInfo& SetAlloc (bool a) {
        alloc = a;
        return *this;
    }

    FabInfo& SetDevice (bool d) {
	device = d;
	return *this;
    }

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
};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
    {
        FAB* fab = new FAB(box, ncomps, info.alloc, info.shared);
	if (info.device) fab->allocate_device();
	return fab;
    }
};

}

#endif
