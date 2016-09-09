
#include <winstd.H>
#include <BndryData.H>
#include <Utility.H>
#include <LO_BCTYPES.H>
#include <ParallelDescriptor.H>

//
// Mask info required for this many cells past grid edge
//  (here, e.g. ref=4, crse stencil width=3, and let stencil slide 2 past grid edge)
//
int BndryData::NTangHalfWidth = 5;  // ref_ratio + 1, so won't work if ref_ratio > 4

BndryData::BndryData ()
    :
    m_ncomp(-1), m_defined(false) {}

BndryData::BndryData (const BoxArray& _grids,
                      int             _ncomp, 
                      const Geometry& _geom,
		      ParallelDescriptor::Color color)
    :
    geom(_geom),
    m_ncomp(_ncomp),
    m_defined(false)
{
    define(_grids,_ncomp,_geom,color);
}

void
BndryData::setBoundCond (Orientation     _face,
                         int              _n,
                         int              _comp,
                         const BoundCond& _bcn)
{
    bcond[_n][_face][_comp] = _bcn;
}

void
BndryData::setBoundLoc (Orientation _face,
                        int         _n,
                        Real        _val)
{
    bcloc[_n][_face] = _val;
}

const Array< Array<BoundCond> >&
BndryData::bndryConds (int igrid) const
{
    std::map< int, Array< Array<BoundCond> > >::const_iterator it = bcond.find(igrid);
    BL_ASSERT(it != bcond.end());
    return it->second;
}

const BndryData::RealTuple&
BndryData::bndryLocs (int igrid) const
{
    std::map<int,RealTuple>::const_iterator it = bcloc.find(igrid);
    BL_ASSERT(it != bcloc.end());
    return it->second;
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
    masks.resize(2*BL_SPACEDIM, PArrayManage);
    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
	const MultiMask& smasks = src.masks[i];
	masks.set(i, new MultiMask(smasks.boxArray(), smasks.DistributionMap(), smasks.nComp()));
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
	for (int i = 0; i < 2*BL_SPACEDIM; i++) {
	    bndry[i].clear();
	}
        init(src);
    }
    return *this;
}

BndryData::~BndryData ()
{
}

void
BndryData::define (const BoxArray& _grids,
                   int             _ncomp,
                   const Geometry& _geom,
		   ParallelDescriptor::Color color)
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
        BoxLib::Abort("BndryData::define(): object already built");
    }
    geom    = _geom;
    m_ncomp = _ncomp;

    BndryRegister::setBoxes(_grids);

    masks.clear();
    masks.resize(2*BL_SPACEDIM, PArrayManage);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();

        BndryRegister::define(face,IndexType::TheCellType(),0,1,1,_ncomp,color);
	
	masks.set(face, new MultiMask(grids, bndry[face].DistributionMap(), geom,
				      face, 0, 2, NTangHalfWidth, 1, true));
    }
    //
    // Define "bcond" and "bcloc".
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    //
    for (FabSetIter bfsi(bndry[Orientation(0,Orientation::low)]);
         bfsi.isValid();
         ++bfsi)
    {
        const int idx = bfsi.index();
        //
        // Insert with a hint since we know the indices are increasing.
        //
        bcloc.insert(bcloc.end(),std::map<int,RealTuple>::value_type(idx,RealTuple()));

        std::map< int, Array< Array<BoundCond> > >::value_type v(idx,Array< Array<BoundCond> >());

        Array< Array<BoundCond> >& abc = bcond.insert(bcond.end(),v)->second;

        abc.resize(2*BL_SPACEDIM);

        for (OrientationIter fi; fi; ++fi)
        {
            abc[fi()].resize(_ncomp);
        }
    }

    m_defined = true;
}

#if 0
//xxxxx
std::ostream&
operator<< (std::ostream&    os,
            const BndryData& bd)
{
    if (ParallelDescriptor::NProcs() != 1)
	BoxLib::Abort("BndryData::operator<<(): not implemented in parallel");

    const BoxArray& grds  = bd.boxes();
    const int       ngrds = grds.size();
    const int       ncomp = bd.nComp();

    os << "[BndryData with " << ngrds << " grids and " << ncomp << " comps:\n";

    for (int grd = 0; grd < ngrds; grd++)
    {
        const BndryData::RealTuple& bdl = bd.bndryLocs(grd);

        const Array< Array<BoundCond> > & bcs = bd.bndryConds(grd);

        const BndryData::MaskTuple& msk = bd.bndryMasks(grd);

        for (OrientationIter face; face; ++face)
        {
            const Orientation f = face();

            os << "::: face " << (int)(f) << " of grid " << grds[grd] << "\nBC = ";

            const Array<BoundCond>& bc = bcs[f];

            for (int i = 0; i < ncomp; ++i)
                os << bc[i] << ' ';
            os << " LOC = " << bdl[f] << '\n';
            os << msk[f];
            os << bd.bndry[f][grd];
        }
        os << "------------------------------------------------" << '\n';
    }

    return os;
}

void 
BndryData::writeOn (std::ostream& os) const
{
    if (ParallelDescriptor::NProcs() != 1)
	BoxLib::Abort("BndryData::writeOn(): not implemented in parallel");

    const int ngrds = grids.size();
    const int ncomp = nComp();

    os << ngrds << " " << ncomp << '\n';

    for (int grd = 0; grd < ngrds; grd++)
    {
        os << grids[grd] << '\n';
    }

    os << geom << '\n';

    for (int grd = 0; grd < ngrds; grd++)
    {
        const BndryData::RealTuple& bdl = bndryLocs(grd);

        const Array< Array<BoundCond> >& bcs = bndryConds(grd);

        for (OrientationIter face; face; ++face)
        {
            const Orientation f = face();

            const Array<BoundCond>& bc = bcs[f];

            for (int cmp = 0; cmp < ncomp; cmp++)
                os << bc[cmp] << ' ';
            os << '\n';

            os << bdl[f] << '\n';
        }
    }

    std::map<int,MaskTuple>::const_iterator it;

    for (OrientationIter face; face; ++face)
    {
        const Orientation f = face();

        for (int grd = 0; grd < ngrds; grd++)
        {
            it = masks.find(grd);
            BL_ASSERT(it != masks.end());
            it->second[f]->writeOn(os);
            bndry[f][grd].writeOn(os);
        }
    }
}

void 
BndryData::readFrom (std::istream& is)
{
    if (ParallelDescriptor::NProcs() != 1)
	BoxLib::Abort("BndryData::readFrom(): not implemented in parallel");

    int tmpNgrids, tmpNcomp;
    is >> tmpNgrids >> tmpNcomp;

    BoxArray tmpBa(tmpNgrids);
    for (int grd = 0; grd < tmpNgrids; grd++)
    {
        Box readBox;
        is >> readBox;
        tmpBa.set(grd, readBox);
    }

    Geometry tmpGeom;
    is >> tmpGeom;

    define(tmpBa, tmpNcomp, tmpGeom);

    for (int grd = 0; grd < tmpNgrids; grd++)
    {
        for (OrientationIter face; face; ++face)
        {
            Orientation f = face();

            int intBc;
            for (int cmp = 0; cmp < tmpNcomp; cmp++)
            {
                is >> intBc;
                setBoundCond(f, grd, cmp, BoundCond(intBc));
            }

            Real realBcLoc;
            is >> realBcLoc;
            setBoundLoc(f, grd, realBcLoc);
        }
    }

    std::map<int,MaskTuple>::const_iterator it;

    for (OrientationIter face; face; ++face)
    {
        const Orientation f = face();

        for (int grd = 0; grd < tmpNgrids; grd++)
        {
            it = masks.find(grd);
            BL_ASSERT(it != masks.end());
            it->second[f]->readFrom(is);
            bndry[f][grd].readFrom(is);
        }
    }
}
#endif
