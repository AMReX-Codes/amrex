//BL_COPYRIGHT_NOTICE

//
// $Id: BndryData.cpp,v 1.4 1998-07-07 17:25:12 lijewski Exp $
//

#include <BndryData.H>
#include <Utility.H>

BndryData::BndryData (const BoxArray&      _grids,
                      int                  _ncomp, 
                      const ProxyGeometry& _geom)
    :
    geom(_geom)
{
    define(_grids,_ncomp,_geom);
}

BndryData::BndryData (const BndryData& src)
{
    (*this) = src;
}

BndryData::~BndryData ()
{
    //
    // Masks was not allocated with PArrayManage, must manually dealloc.
    //
    clear_masks();
}

void
BndryData::clear_masks ()
{
    for (OrientationIter oitr; oitr; oitr++)
    {
        Orientation face = oitr();
        for (int k = 0; k < masks[face].length(); k++)
        {
            if (masks[face].defined(k))
            {
                delete masks[face].remove(k);
            }
        }
    }
}

ostream&
operator<< (ostream&         os,
            const BndryData& bd)
{
    ParallelDescriptor::Abort("BndryData::operator<< not yet implemented in parallel");
    const BoxArray& grds = bd.boxes();
    int ngrds = grds.length();
    int ncomp = bd.bcond[0][0].length();
    os << "[BndryData with " << ngrds << " grids and "<<ncomp<<" comps:\n";
    for (int grd = 0; grd < ngrds; grd++)
    {
    for (OrientationIter face; face; ++face)
    {
        Orientation f = face();
        os << "::: face " << f << " of grid " << grds[grd] << "\n";
        os << "BC = " ;
        for( int i=0; i<ncomp; ++i)
        {
        os << bd.bcond[f][grd][i] << " ";
        }
        os << " LOC = " << bd.bcloc[f][grd] << "\n";
        os << bd.masks[f][grd];
        os << bd.bndry[f][grd];
    }
    os << "--------------------------------------------------" << endl;
    }
    return os;
}

BndryData&
BndryData::operator= (const BndryData& src)
{
    //
    // Got to save the geometric info.
    //
    geom = src.geom;
    //
    // Redefine grids and bndry array.
    //
    BndryRegister::operator= ( (BndryRegister) src);
    int ngrd = grids.length();
    int ncomp = src.bcond[0][0].length();
    clear_masks();
    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();
        bcond[face].resize(ngrd);
        for (int grd = 0; grd < ngrd; ++grd)
        {
            bcond[face][grd].resize(ncomp);
        }
        bcloc[face].resize(ngrd);
        masks[face].resize(ngrd);
        for (ConstFabSetIterator bfsi(bndry[face]); bfsi.isValid(false); ++bfsi)
        {
            const int grd = bfsi.index();
            bcond[face][grd] = src.bcond[face][grd];
            bcloc[face][grd] = src.bcloc[face][grd];
            const Mask& src_mask = src.masks[face][grd];
            Mask* m = new Mask(src_mask.box(),src_mask.nComp());
            m->copy(src_mask);
            masks[face].set(grd,m);
        }
    }
    return *this;
}

void
BndryData::define (const BoxArray&      _grids,
                   int                  _ncomp,
                   const ProxyGeometry& _geom)
{
    geom = _geom;
    BndryRegister::setBoxes(_grids);
    const int ngrd = grids.length();
    assert(ngrd > 0);

    Array<IntVect> pshifts(27);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();
        const int coord_dir = face.coordDir();
        masks[face].resize(ngrd);
        bcloc[face].resize(ngrd);
        bcond[face].resize(ngrd);
        for (int ig = 0; ig < ngrd; ++ig)
        {
            bcond[face][ig].resize(_ncomp);
        }
        BndryRegister::define(face,IndexType::TheCellType(),0,1,0,_ncomp);
        //
        // Alloc mask and set to quad_interp value.
        //
        for (ConstFabSetIterator bfsi(bndry[face]); bfsi.isValid(false); ++bfsi)
        {
            Box face_box = ::adjCell(grids[bfsi.index()], face, 1);
            //
            // Extend box in directions orthogonal to face normal.
            //
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                if (dir == coord_dir)
                    continue;
                face_box.grow(dir,1);
            }
            Mask* m = new Mask(face_box);
            m->setVal(outside_domain,0);
            Box dbox = geom.Domain() & face_box;
            m->setVal(not_covered,dbox,0);
            //
            // Now have to set as not_covered the periodic translates as well.
            //
            if (geom.isAnyPeriodic())
            {
                geom.periodicShift(geom.Domain(), face_box, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    m->shift(pshifts[iiv]);
                    Box target = geom.Domain() & m->box();
                    m->setVal(not_covered,target,0);
                    m->shift(-pshifts[iiv]);
                }
            }
            masks[face].set(bfsi.index(),m);
            //
            // Turn mask off on intersection with grids at this level.
            //
            for (int g = 0; g < ngrd; g++)
            {
                if (grids[g].intersects(face_box))
                {
                    Box ovlp = grids[g] & face_box;
                    m->setVal(covered,ovlp,0);
                }
            }
            //
            // Handle special cases if is periodic.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(face_box))
            {
                geom.periodicShift(geom.Domain(), face_box, pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    m->shift(pshifts[iiv]);
                    for (int g = 0; g < ngrd; g++)
                    {
                        if (grids[g].intersects(m->box()))
                        {
                            Box ovlp = grids[g] & m->box();
                            m->setVal(covered,ovlp,0);
                        }
                    }
                    m->shift(-pshifts[iiv]);
                }
            }
        }
    }
}
