//BL_COPYRIGHT_NOTICE

//
// $Id: BndryRegister.cpp,v 1.10 1999-05-10 17:18:29 car Exp $
//

#include <BndryRegister.H>

const Real BL_SAFE_BOGUS = -666.e30;

BndryRegister::BndryRegister () {}

BndryRegister::~BndryRegister () {}

BndryRegister::BndryRegister (const BoxArray& _grids,
                              int             _in_rad,
                              int             _out_rad,
                              int             _extent_rad,
                              int             _ncomp)
    :
    grids(_grids)
{
    BLassert(_ncomp > 0);
    BLassert(grids.ready());
    BLassert(grids[0].cellCentered());

    for (OrientationIter face; face; ++face)
    {
        define(face(),
               IndexType::TheCellType(),
               _in_rad,
               _out_rad,
               _extent_rad,
               _ncomp);
    }
}

void
BndryRegister::init (const BndryRegister& src)
{
    grids.define(src.grids);

    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        FabSet& fs = bndry[i];

        fs.define(src.bndry[i].boxArray(), src.bndry[i].nComp());

        for (ConstFabSetIterator mfi(src.bndry[i]); mfi.isValid(); ++mfi)
        {
            fs[mfi.index()].copy(mfi());
        }
    }
}

BndryRegister::BndryRegister (const BndryRegister& src)
{
    init(src);
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

    init(src);

    return *this;
}

void
BndryRegister::define (const Orientation& _face,
                       const IndexType&   _typ,
                       int                _in_rad,
                       int                _out_rad,
                       int                _extent_rad,
                       int                _ncomp)
{
    BLassert(grids.ready());

    FabSet& fabs = bndry[_face];

    BLassert(!fabs.ready());

    const int coord_dir = _face.coordDir();
    const int lo_side   = _face.isLow();
    //
    // Build the BoxArray on which to define the FabSet on this face.
    //
    BoxArray fsBA(grids.length());

    for (int idx = 0; idx < grids.length(); ++idx)
    {
        Box b;
        //
        // First construct proper box for direction normal to face.
        //
        if (_out_rad > 0)
        {
            if (_typ.ixType(coord_dir) == IndexType::CELL)
                b = ::adjCell(grids[idx], _face, _out_rad);
            else
                b = ::bdryNode(grids[idx], _face, _out_rad);

            if (_in_rad > 0)
                b.grow(_face.flip(), _in_rad);
        }
        else
        {
            if (_in_rad > 0)
            {
                if (_typ.ixType(coord_dir) == IndexType::CELL)
                    b = ::adjCell(grids[idx], _face, _in_rad);
                else
                    b = ::bdryNode(grids[idx], _face, _in_rad);

                b.shift(coord_dir, lo_side ? _in_rad : -_in_rad);
            }
            else
                BoxLib::Error("BndryRegister::define(): strange values for in_rad, out_rad");
        }
        //
        // Now alter box in all other index directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (dir == coord_dir)
                continue;
            if (_typ.ixType(dir) == IndexType::NODE)
                b.surroundingNodes(dir);
            if (_extent_rad > 0)
                b.grow(dir,_extent_rad);
        }

        BLassert(b.ok());

        fsBA.set(idx,b);
    }

    BLassert(fsBA.ok());

    fabs.define(fsBA,_ncomp);

    // 
    // Go ahead and assign values to the boundary register fabs
    // since in some places APPLYBC (specifically in the tensor
    // operator) the boundary registers are used for a few calculations
    // before the masks are tested to see if you need them.
    //
    fabs.setVal(BL_SAFE_BOGUS);
}

void
BndryRegister::setBoxes (const BoxArray& _grids)
{
    BLassert(!grids.ready());
    BLassert(_grids.ready());
    BLassert(_grids[0].cellCentered());

    grids.define(_grids);
    //
    // Check that bndry regions are not allocated.
    //
    for (int k = 0; k < 2*BL_SPACEDIM; k++)
        BLassert(!bndry[k].ready());
}

void BndryRegister::setVal (Real v)
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
        bndry[face()].linComb(a,
                              mfa,
                              a_comp,
                              b,
                              mfb,
                              b_comp,
                              dest_comp,
                              num_comp,
                              n_ghost);
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
