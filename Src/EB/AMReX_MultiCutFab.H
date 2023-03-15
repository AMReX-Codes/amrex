#ifndef AMREX_MULTICUTFAB_H_
#define AMREX_MULTICUTFAB_H_
#include <AMReX_Config.H>

#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class CutFab final
    : public FArrayBox
{
public:

    CutFab () = default;

    explicit CutFab (Arena* ar) : FArrayBox(ar) {}

    CutFab (const Box& b, int n, Arena* ar)
        : FArrayBox(b, n, ar) {}

    CutFab (const Box& b, int ncomps, bool alloc=true,
            bool shared=false, Arena* ar = nullptr)
        : FArrayBox(b, ncomps, alloc, shared, ar) {}

    CutFab (CutFab const& rhs, MakeType make_type, int scomp, int ncomp)
        : FArrayBox(rhs, make_type, scomp, ncomp) {}

    ~CutFab () noexcept override = default;

    CutFab (CutFab&& rhs) noexcept = default;

    CutFab (const CutFab&) = delete;
    CutFab& operator= (const CutFab&) = delete;
    CutFab& operator= (CutFab&&) = delete;

    template <RunOn run_on, typename BUF = CutFab::value_type>
    std::size_t copyFromMem (const void* src) {
        return copyFromMem<run_on, BUF>(box(), 0, nComp(), src);
    }

    template <RunOn run_on, typename BUF = CutFab::value_type>
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src)
    {
        if (dptr != nullptr) {
            return FArrayBox::copyFromMem<run_on, BUF>(dstbox, dstcomp, numcomp, src);
        } else {
            return sizeof(BUF)*static_cast<std::size_t>(dstbox.numPts()*numcomp);
        }
    }

    template <RunOn run_on>
    CutFab& copy (const CutFab & src,
                  const Box&     srcbox,
                  int            srccomp,
                  const Box&     destbox,
                  int            destcomp,
                  int            numcomp)
    {
        if (dptr != nullptr) {
            FArrayBox::copy<run_on>(src,srcbox,srccomp,destbox,destcomp,numcomp);
        }
        return *this;
    }

    template <RunOn run_on>
    CutFab& copy (const CutFab & src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp)
    {
        if (dptr != nullptr) {
            FArrayBox::copy<run_on>(src, bx, scomp, dcomp, ncomp);
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

    const CutFab& operator[] (int global_box_index) const noexcept;
    CutFab& operator[] (int global_box_index) noexcept;

    Array4<Real      > array (const MFIter& mfi) noexcept;
    Array4<Real const> array (const MFIter& mfi) const noexcept;
    Array4<Real const> const_array (const MFIter& mfi) const noexcept;

    //! Is it OK to call operator[] with this MFIter?
    bool ok (const MFIter& mfi) const noexcept;

    //! Is it OK to call operator[] with this global box index?
    bool ok (int global_box_index) const noexcept;

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
