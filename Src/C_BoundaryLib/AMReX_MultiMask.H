#ifndef _MULTIMASK_H_
#define _MULTIMASK_H_

#include <Mask.H>
#include <FabArray.H>
#include <Geometry.H>

#ifdef _OPENMP
#include <omp.h>
#endif

class MultiMaskIter;

class MultiMask
{
public:
    MultiMask () { ; }
    MultiMask (const BoxArray& ba, const DistributionMapping& dm, int ncomp);
    MultiMask (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
	       Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval);

    Mask& operator[] (const MFIter& mfi) { return m_fa[mfi]; }
    const Mask& operator[] (const MFIter& mfi) const { return m_fa[mfi]; }

    int nComp () const { return m_fa.nComp(); }

    const BoxArray& boxArray () const { return m_fa.boxArray(); }

    const DistributionMapping& DistributionMap () const { return m_fa.DistributionMap(); }

    static void Copy (MultiMask& dst, const MultiMask& src);

    friend class MultiMaskIter;

private:

    FabArray<Mask> m_fa;

    //
    // These are disabled.
    //
    MultiMask (const MultiMask& rhs);
    MultiMask& operator= (const MultiMask& rhs);
};

class MultiMaskIter
    : public MFIter
{
public:
    explicit MultiMaskIter (const MultiMask& mm) 
	: MFIter(mm.m_fa) { }
};

#endif

