
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>
#include <AMReX_Print.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBISLayout.H>
#endif

namespace amrex
{

class FArrayBox;

struct FabInfo
{
    bool alloc = true;
    bool shared = false;

    FabInfo& SetAlloc (bool a) {
        alloc = a;
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
        return new FAB(box, ncomps, info.alloc, info.shared);
    }
};

#ifdef AMREX_USE_EB

class FArrayBoxFactory
    : public FabFactory<FArrayBox>
{
public:

    FArrayBoxFactory () = default;
    FArrayBoxFactory (const EBISLayout& a_ebisl);
    virtual ~FArrayBoxFactory () = default;

    FArrayBoxFactory (const FArrayBoxFactory&) = delete;
    FArrayBoxFactory (FArrayBoxFactory&&) = delete;
    FArrayBoxFactory& operator= (const FArrayBoxFactory&) = delete;
    FArrayBoxFactory& operator= (FArrayBoxFactory&&) = delete;


    virtual FArrayBox* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

private:
    EBISLayout m_ebisl;
};

#else

using FArrayBoxFactory = DefaultFabFactory<FArrayBox>;

#endif

}

#endif
