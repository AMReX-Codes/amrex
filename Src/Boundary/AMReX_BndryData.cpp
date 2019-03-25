
#include <AMReX_BndryData.H>
#include <AMReX_Utility.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {

//
// Mask info required for this many cells past grid edge
//  (here, e.g. ref=4, crse stencil width=3, and let stencil slide 2 past grid edge)
//
int BndryData::NTangHalfWidth = 5;  // ref_ratio + 1, so won't work if ref_ratio > 4

BndryData::BndryData () noexcept
    :
    m_ncomp(-1), m_defined(false) {}

BndryData::BndryData (const BoxArray& _grids,
		      const DistributionMapping& _dmap,
                      int             _ncomp, 
                      const Geometry& _geom)
    :
    geom(_geom),
    m_ncomp(_ncomp),
    m_defined(false)
{
    define(_grids,_dmap,_ncomp,_geom);
}

void
BndryData::setBoundCond (Orientation     _face,
                         int              _n,
                         int              _comp,
                         const BoundCond& _bcn) noexcept
{
    bcond[_n][_face][_comp] = _bcn;
}

void
BndryData::setBoundCond (Orientation     _face,
                         const MFIter&   mfi,
                         int              _comp,
                         const BoundCond& _bcn) noexcept
{
    bcond[mfi][_face][_comp] = _bcn;
}

void
BndryData::setBoundLoc (Orientation _face,
                        int         _n,
                        Real        _val) noexcept
{
    bcloc[_n][_face] = _val;
}

void
BndryData::setBoundLoc (Orientation _face,
                        const MFIter& mfi,
                        Real        _val) noexcept
{
    bcloc[mfi][_face] = _val;
}

const Vector< Vector<BoundCond> >&
BndryData::bndryConds (int igrid) const noexcept
{
    return bcond[igrid];
}

const Vector< Vector<BoundCond> >&
BndryData::bndryConds (const MFIter& mfi) const noexcept
{
    return bcond[mfi];
}

const BndryData::RealTuple&
BndryData::bndryLocs (int igrid) const noexcept
{
    return bcloc[igrid];
}

const BndryData::RealTuple&
BndryData::bndryLocs (const MFIter& mfi) const noexcept
{
    return bcloc[mfi];
}

void
BndryData::init (const BndryData& src)
{
    geom      = src.geom;
    m_ncomp   = src.m_ncomp;
    m_defined = src.m_defined;
    bcloc     = src.bcloc;
    bcond     = src.bcond;

    masks.clear();
    masks.resize(2*AMREX_SPACEDIM);
    for (int i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
	const MultiMask& smasks = src.masks[i];
	masks[i].define(smasks.boxArray(), smasks.DistributionMap(), smasks.nComp());
	MultiMask::Copy(masks[i], smasks);
    }
}

BndryData::BndryData (const BndryData& src)
    :
    BndryRegister(src),
    m_ncomp(src.m_ncomp)
{
    init(src);
}

BndryData&
BndryData::operator= (const BndryData& src)
{
    if (this != &src)
    {
        BndryRegister::operator=(src);
        init(src);
    }
    return *this;
}

BndryData::~BndryData ()
{
}

void
BndryData::define (const BoxArray& _grids,
		   const DistributionMapping& _dmap,
                   int             _ncomp,
                   const Geometry& _geom)
{
    BL_PROFILE("BndryData::define()");

    if (m_defined)
    {
        if (_grids == boxes() && m_ncomp == _ncomp && _geom.Domain() == geom.Domain())
            //
            // We want to allow reuse of BndryData objects that were define()d exactly as a previous call.
            //
            return;
        //
        // Otherwise we'll just abort.  We could make this work but it's just as easy to start with a fresh Bndrydata object.
        //
        amrex::Abort("BndryData::define(): object already built");
    }
    geom    = _geom;
    m_ncomp = _ncomp;

    BndryRegister::setBoxes(_grids);

    masks.clear();
    masks.resize(2*AMREX_SPACEDIM);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();
        BndryRegister::define(face,IndexType::TheCellType(),0,1,1,_ncomp,_dmap);
	masks[face].define(grids, _dmap, geom, face, 0, 2, NTangHalfWidth, 1, true);
    }
    //
    // Define "bcond" and "bcloc".
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    //

    bcloc.define(grids, _dmap);
    bcond.define(grids, _dmap);

    for (FabSetIter bfsi(bndry[Orientation(0,Orientation::low)]);
         bfsi.isValid();
         ++bfsi)
    {
        Vector< Vector<BoundCond> >& abc = bcond[bfsi];

        abc.resize(2*AMREX_SPACEDIM);

        for (OrientationIter fi; fi; ++fi)
        {
            abc[fi()].resize(_ncomp);
        }
    }

    m_defined = true;
}

void
BndryData::setValue (Orientation face, int n, Real val) noexcept
{
    auto& fab = bndry[face][n]; 
    auto arr = fab.array();
    const Box& bx = fab.box();
    const int ncomp = m_ncomp;
    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, m,
    {
        arr(i,j,k,m) = val;
    });
}

}
