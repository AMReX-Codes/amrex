#ifndef AMREX_MULTICUTFAB_H_
#define AMREX_MULTICUTFAB_H_

#include <AMReX_FabArray.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class CutFab
    : public FArrayBox
{
public:

    CutFab () {}
    CutFab (const Box& b, int ncomps, bool alloc=true, bool shared=false)
        : FArrayBox(b, ncomps, alloc, shared) {}

#ifdef AMREX_USE_GPU
    CutFab (CutFab const& rhs, MakeType make_type, int scomp, int ncomp)
        : FArrayBox(rhs, make_type, scomp, ncomp) {}
#endif

    ~CutFab () = default;

    CutFab (CutFab&& rhs) noexcept = default;

    CutFab (const CutFab&) = delete;
    CutFab& operator= (const CutFab&) = delete;
    CutFab& operator= (CutFab&&) = delete;

    AMREX_GPU_HOST_DEVICE
    std::size_t copyFromMem (const void* src) {
        return copyFromMem (box(), 0, nComp(), src);
    }

    AMREX_GPU_HOST_DEVICE
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src)
    {
        if (dptr != nullptr) {
            return FArrayBox::copyFromMem(dstbox, dstcomp, numcomp, src);
        } else {
            return sizeof(CutFab::value_type)*static_cast<std::size_t>(dstbox.numPts()*numcomp);
        }
    }

    AMREX_GPU_HOST_DEVICE
    CutFab& copy (const CutFab & src,
                  const Box&     srcbox,
                  int            srccomp,
                  const Box&     destbox,
                  int            destcomp,
                  int            numcomp)
    {
        if (dptr != nullptr) {
            FArrayBox::copy(src,srcbox,srccomp,destbox,destcomp,numcomp);
        }
        return *this;
    }

    AMREX_GPU_HOST_DEVICE
    CutFab& copy (const CutFab & src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp)
    {
        if (dptr != nullptr) {
            FArrayBox::copy(src, bx, scomp, dcomp, ncomp);
        }
        return *this;
    }
};

class MultiCutFab
{
public:

    MultiCutFab ();

    MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                 int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags);

    ~MultiCutFab ();

    MultiCutFab (MultiCutFab&& rhs) noexcept = default;
    
    MultiCutFab (const MultiCutFab& rhs) = delete;
    MultiCutFab& operator= (const MultiCutFab& rhs) = delete;
    MultiCutFab& operator= (MultiCutFab&& rhs) = delete;

    void define (const BoxArray& ba, const DistributionMapping& dm,
                 int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags);

    const CutFab& operator[] (const MFIter& mfi) const noexcept;
    CutFab& operator[] (const MFIter& mfi) noexcept;

    CutFab      * fabPtr (const MFIter& mfi) noexcept;
    CutFab const* fabPtr (const MFIter& mfi) const noexcept;

    Array4<Real      > array (const MFIter& mfi) noexcept;
    Array4<Real const> array (const MFIter& mfi) const noexcept;

    bool ok (const MFIter& mfi) const noexcept;

    void setVal (Real val);

    FabArray<CutFab>& data () noexcept { return m_data; }
    const FabArray<CutFab>& data () const noexcept { return m_data; }

    const BoxArray& boxArray () const noexcept { return m_data.boxArray(); }
    const DistributionMapping& DistributionMap () const noexcept { return m_data.DistributionMap(); }
    int nComp () const noexcept { return m_data.nComp(); }
    int nGrow () const noexcept { return m_data.nGrow(); }

    void ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng,
                       const Periodicity& period = Periodicity::NonPeriodic());

    MultiFab ToMultiFab (Real regular_value, Real covered_value) const;

private:

    FabArray<CutFab> m_data;
    const FabArray<EBCellFlagFab>* m_cellflags;

    void remove ();
};

}

#endif
