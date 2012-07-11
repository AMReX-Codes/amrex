
#include <winstd.H>
#include <BndryData.H>
#include <Utility.H>
#include <LO_BCTYPES.H>
#include <ParallelDescriptor.H>

//
// Mask info required for this many cells past grid edge
//  (here, e.g. ref=4, crse stencil width=3, and let stencil slide 2 past grid edge)
//
int BndryData::NTangHalfWidth = 5;

BndryData::BndryData()
    :
    m_ncomp(-1) {}

BndryData::BndryData (const BoxArray& _grids,
                      int             _ncomp, 
                      const Geometry& _geom)
    :
    geom(_geom),
    m_ncomp(_ncomp)
{
    define(_grids,_ncomp,_geom);
}

const Array<BoundCond>&
BndryData::bndryConds (const Orientation& _face, int igrid) const
{
    std::map< int,Array<BoundCond> >::const_iterator it = bcond[_face].find(igrid);

    BL_ASSERT(it != bcond[_face].end());

    return it->second;
}

void
BndryData::init (const BndryData& src)
{
    //
    // Got to save the geometric info.
    //
    m_ncomp = src.m_ncomp;

    geom = src.geom;
    //
    // Redefine grids and bndry array.
    //
    const int ngrd  = grids.size();

    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face = fi();

        bcond[face] = src.bcond[face];

        bcloc[face] = src.bcloc[face];

        masks[face].resize(ngrd);

        for (FabSetIter bfsi(bndry[face]); bfsi.isValid(); ++bfsi)
        {
            const int grd        = bfsi.index();
            const Mask& src_mask = src.masks[face][grd];
            Mask* m = new Mask(src_mask.box(),src_mask.nComp());
            m->copy(src_mask);
            masks[face].set(grd,m);
        }
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
        clear_masks();
        init(src);
    }
    return *this;
}

BndryData::~BndryData ()
{
    //
    // Masks was not allocated with PArrayManage -- manually dealloc.
    //
    clear_masks();
}

void
BndryData::clear_masks ()
{
    for (OrientationIter oitr; oitr; oitr++)
    {
        const Orientation face = oitr();

        for (int k = 0, N = masks[face].size(); k < N; k++)
        {
            if (masks[face].defined(k))
            {
                delete masks[face].remove(k);
            }
        }
    }
}

void
BndryData::define (const BoxArray& _grids,
                   int             _ncomp,
                   const Geometry& _geom)
{
    m_ncomp = _ncomp;

    geom = _geom;

    BndryRegister::setBoxes(_grids);

    const int ngrd = grids.size();

    BL_ASSERT(ngrd > 0);

    Array<IntVect> pshifts(27);

    std::vector< std::pair<int,Box> > isects;

    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face      = fi();
        const int         coord_dir = face.coordDir();

        masks[face].resize(ngrd);
        bcloc[face].resize(ngrd);

        BndryRegister::define(face,IndexType::TheCellType(),0,1,0,_ncomp);
        //
        // Alloc mask and set to quad_interp value.
        //
        for (FabSetIter bfsi(bndry[face]); bfsi.isValid(); ++bfsi)
        {
            const int igrid = bfsi.index();

            bcond[face][igrid].resize(_ncomp);

            Box face_box = BoxLib::adjCell(grids[igrid], face, 1);
            //
            // Extend box in directions orthogonal to face normal.
            //
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                if (dir == coord_dir) continue;
                face_box.grow(dir,NTangHalfWidth);
            }
            Mask* m = new Mask(face_box);
            m->setVal(outside_domain,0);
            const Box dbox = geom.Domain() & face_box;
            m->setVal(not_covered,dbox,0);
            //
            // Now have to set as not_covered the periodic translates as well.
            //
            if (geom.isAnyPeriodic())
            {
                geom.periodicShift(geom.Domain(), face_box, pshifts);

                for (int iiv = 0, N = pshifts.size(); iiv < N; iiv++)
                {
                    m->shift(pshifts[iiv]);
                    const Box target = geom.Domain() & m->box();
                    m->setVal(not_covered,target,0);
                    m->shift(-pshifts[iiv]);
                }
            }
            masks[face].set(igrid,m);
            //
            // Turn mask off on intersection with grids at this level.
            //
            grids.intersections(face_box,isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                m->setVal(covered, isects[ii].second, 0);
            }
            //
            // Handle special cases if is periodic.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(face_box))
            {
                geom.periodicShift(geom.Domain(), face_box, pshifts);

                for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
                {
                    m->shift(pshifts[iiv]);

                    grids.intersections(m->box(),isects);

                    for (int ii = 0, N = isects.size(); ii < N; ii++)
                        m->setVal(covered, isects[ii].second, 0);

                    m->shift(-pshifts[iiv]);
                }
            }
        }
    }
}

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
        for (OrientationIter face; face; ++face)
        {
            Orientation f = face();

            os << "::: face " << (int)(f) << " of grid " << grds[grd] << "\nBC = ";

            const Array<BoundCond>& bc = bd.bndryConds(f,grd);

            for (int i = 0; i < ncomp; ++i)
                os << bc[i] << ' ';

            os << " LOC = " << bd.bcloc[f][grd] << '\n';
            os << bd.masks[f][grd];
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

    int ngrds = grids.size();
    int ncomp = nComp();

    os << ngrds << " " << ncomp << '\n';

    for (int grd = 0; grd < ngrds; grd++)
    {
        os << grids[grd] << '\n';
    }

    os << geom << '\n';

    for (int grd = 0; grd < ngrds; grd++)
    {
        for (OrientationIter face; face; ++face)
        {
            Orientation f = face();

            const Array<BoundCond>& bc = bndryConds(f,grd);

            for (int cmp = 0; cmp < ncomp; cmp++)
                os << bc[cmp] << ' ';
            os << '\n';

            os << bcloc[f][grd] << '\n';
        }
    }

    for (OrientationIter face; face; ++face)
    {
        Orientation f = face();

        for (int grd = 0; grd < ngrds; grd++)
        {
            masks[f][grd].writeOn(os);
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

    for (OrientationIter face; face; ++face)
    {
        Orientation f = face();

        for (int grd = 0; grd < tmpNgrids; grd++)
        {
            masks[f][grd].readFrom(is);
            bndry[f][grd].readFrom(is);
        }
    }
}
