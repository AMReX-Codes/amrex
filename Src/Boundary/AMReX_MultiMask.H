#ifndef AMREX_MULTIMASK_H_
#define AMREX_MULTIMASK_H_

#include <AMReX_Mask.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

class MultiMaskIter;

class MultiMask
{
public:
    MultiMask () { ; }
    MultiMask (const BoxArray& ba, const DistributionMapping& dm, int ncomp);
    MultiMask (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
	       Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval);

    ~MultiMask () = default;

    MultiMask (MultiMask&& rhs) noexcept = default;

    MultiMask (const MultiMask& rhs) = delete;
    MultiMask& operator= (const MultiMask& rhs) = delete;
    MultiMask& operator= (MultiMask&& rhs) = delete;

    void define (const BoxArray& ba, const DistributionMapping& dm, int ncomp);
    void define (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
		 Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval);

    Mask& operator[] (const MFIter& mfi) noexcept { return m_fa[mfi]; }
    const Mask& operator[] (const MFIter& mfi) const  noexcept { return m_fa[mfi]; }

    Array4<int const> array (const MFIter& mfi) const noexcept { return m_fa.array(mfi); }
    Array4<int      > array (const MFIter& mfi)       noexcept { return m_fa.array(mfi); }

    int nComp () const noexcept { return m_fa.nComp(); }

    const BoxArray& boxArray () const noexcept { return m_fa.boxArray(); }

    const DistributionMapping& DistributionMap () const noexcept { return m_fa.DistributionMap(); }

    static void Copy (MultiMask& dst, const MultiMask& src);

    friend class MultiMaskIter;

private:

    FabArray<Mask> m_fa;
};

class MultiMaskIter
    : public MFIter
{
public:
    explicit MultiMaskIter (const MultiMask& mm) 
	: MFIter(mm.m_fa) { }
};

}

#endif

