//BL_COPYRIGHT_NOTICE

//
// $Id: BndryRegister.cpp,v 1.6 1998-03-30 20:24:11 lijewski Exp $
//

#include <BndryRegister.H>

BndryRegister::BndryRegister () {}

BndryRegister::~BndryRegister () {}

BndryRegister::BndryRegister (const BndryRegister& src)
{
    grids.define(src.grids);
    int ngrd = grids.length();
    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        bndry[i].resize(ngrd);
        const FabSet& srcfs = src.bndry[i];
        bndry[i].DefineGrids(grids);
        bndry[i].DefineDistributionMap(grids);
        for (ConstFabSetIterator mfi(srcfs); mfi.isValid(); ++mfi)
        {
            FArrayBox* fab = new FArrayBox(mfi().box(),mfi().nComp());
            fab->copy(mfi());
            bndry[i].setFab(mfi.index(),fab);
        }
    }
}

BndryRegister::BndryRegister (const BoxArray& _grids,
                              int             _in_rad,
                              int             _out_rad,
                              int             _extent_rad,
                              int             _ncomp)
    :
    grids(_grids)
{
    assert(grids.ready());
    assert(grids[0].cellCentered());
    assert(_ncomp > 0);

    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),_in_rad,
               _out_rad,_extent_rad,_ncomp);
    }
}

BndryRegister&
BndryRegister::operator= (const BndryRegister& src)
{
    if (grids.ready())
    {
        grids.clear();
        for (int i = 0; i < 2*BL_SPACEDIM; i++)
            bndry[i].clear();
    }
    grids.define(src.grids);
    int ngrd = grids.length();
    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        bndry[i].resize(ngrd);
        bndry[i].DefineGrids(grids);
        bndry[i].DefineDistributionMap(grids);
        const FabSet& srcfs = src.bndry[i];
        for (ConstFabSetIterator mfi(srcfs); mfi.isValid(); ++mfi)
        {
            FArrayBox* fab = new FArrayBox(mfi().box(), mfi().nComp());
            fab->copy(mfi());
            bndry[i].setFab(mfi.index(),fab);
        }
    }
    return *this;
}

void
BndryRegister::setBoxes (const BoxArray& _grids)
{
    assert(!grids.ready());
    assert(_grids.ready());
    assert(_grids[0].cellCentered());
    grids.define(_grids);
    //
    // Insure bndry regions are not allocated.
    //
    for (int k = 0; k < 2*BL_SPACEDIM; k++)
    {
        if (bndry[k].ready())
            bndry[k].clear();
    }
}

void
BndryRegister::define (const Orientation& _face,
                       const IndexType&   _typ,
                       int                _in_rad,
                       int                _out_rad,
                       int                _extent_rad,
                       int                _ncomp)
{
    int myproc = ParallelDescriptor::MyProc();
    assert(grids.ready());
    int ngrd = grids.length();
    FabSet &fabs = bndry[_face];
    assert( !fabs.ready());
    fabs.resize(ngrd);
    fabs.DefineGrids(grids);
    fabs.DefineDistributionMap(grids);
    int coord_dir = _face.coordDir();
    int lo_side = _face.isLow();
    //
    // Don't use a FabSetIterator here.
    //
    for (int mfiindex = 0; mfiindex < grids.length(); ++mfiindex)
    {
        //
        // TODO -- get rid of mfi.index() here.
        //
        Box b;
        //
        // First construct proper box for direction normal to face.
        //
        if (_out_rad > 0)
        {
            if (_typ.ixType(coord_dir) == IndexType::CELL)
            {
                b = adjCell(grids[mfiindex], _face, _out_rad);
            }
            else
            {
                b = bdryNode(grids[mfiindex], _face, _out_rad);
            }
            if (_in_rad > 0)
            {
                //
                // Grow in opposite direction to face.
                //
                Orientation opposite = _face.flip();
                b.grow(opposite, _in_rad);
            }
        }
        else
        {
            if (_in_rad > 0)
            {
                //
                // adjCells in opposite direction to face.
                //
                if (_typ.ixType(coord_dir) == IndexType::CELL)
                {
                    b = adjCell(grids[mfiindex], _face, _in_rad);
                }
                else
                {
                    b = bdryNode(grids[mfiindex], _face, _in_rad);
                }
                b.shift(coord_dir, lo_side?_in_rad:-_in_rad);
            }
            else
            {
                BoxLib::Error("strange values for in_rad, out_rad");
            }
        }
        //
        // Now alter box in all other index directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (dir == coord_dir)
                continue;
            if (_typ.ixType(dir) == IndexType::NODE)
            {
                b.surroundingNodes(dir);
            }
            if (_extent_rad > 0)
                b.grow(dir,_extent_rad);
        }
        assert( b.ok() );
        fabs.setBox(mfiindex, b);
        if (fabs.DistributionMap()[mfiindex] == myproc)
        {
            //
            // Local.
            //
            assert( ! fabs.defined(mfiindex) );
            fabs.clear(mfiindex);
            FArrayBox* fab = new FArrayBox(b,_ncomp);
            fabs.setFab(mfiindex,fab);
        }
    }
}

void BndryRegister::setVal(Real v)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].setVal(v);
    }
}

BndryRegister&
BndryRegister::linComb (Real            a,
                        const MultiFab& mfa,
                        int             a_comp,
                        Real            b,
                        const MultiFab& mfb,
                        int             b_comp,
                        int             dest_comp,
                        int             num_comp,
                        int             n_ghost)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].linComb(a,mfa,a_comp,b,mfb,b_comp,
                              dest_comp,num_comp,n_ghost);
    }
    return *this;
}

BndryRegister&
BndryRegister::copyFrom (const MultiFab& src,
                         int             nghost,
                         int             src_comp,
                         int             dest_comp,
                         int             num_comp)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyFrom(src,nghost,src_comp,dest_comp,num_comp);
    }
    return *this;
}

BndryRegister&
BndryRegister::plusFrom (const MultiFab& src,
                         int             nghost,
                         int             src_comp,
                         int             dest_comp,
                         int             num_comp)
{
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(src,nghost,src_comp,dest_comp,num_comp);
    }
    return *this;
}
