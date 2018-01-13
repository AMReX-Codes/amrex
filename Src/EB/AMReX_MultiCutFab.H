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

    virtual ~CutFab () = default;

    CutFab (CutFab&& rhs) noexcept = default;

    CutFab (const CutFab&) = delete;
    CutFab& operator= (const CutFab&) = delete;
    CutFab& operator= (CutFab&&) = delete;

    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    CutFab& copy (const CutFab & src,
                  const Box&     srcbox,
                  int            srccomp,
                  const Box&     destbox,
                  int            destcomp,
                  int            numcomp);
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

    const CutFab& operator[] (const MFIter& mfi) const;
    CutFab& operator[] (const MFIter& mfi);

    bool ok (const MFIter& mfi) const;

    void setVal (Real val);

    FabArray<CutFab>& data () { return m_data; }
    const FabArray<CutFab>& data () const { return m_data; }

    int nComp () const { return m_data.nComp(); }
    int nGrow () const { return m_data.nGrow(); }

    void ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng);
    void ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng, const Periodicity& period);
    void FillBoundary (const Periodicity& period);

private:

    FabArray<CutFab> m_data;
    const FabArray<EBCellFlagFab>* m_cellflags;

    void remove ();
};

}

#endif
