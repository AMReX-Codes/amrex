//
// $Id: BndryRegister.cpp,v 1.17 2005-01-07 17:38:55 lijewski Exp $
//
#include <winstd.H>

#include <BndryRegister.H>
#include <Orientation.H>

const Real BL_SAFE_BOGUS = -666.e200;

BndryRegister::BndryRegister () {}

BndryRegister::~BndryRegister () {}

BndryRegister::BndryRegister (const BoxArray& grids,
                              int             in_rad,
                              int             out_rad,
                              int             extent_rad,
                              int             ncomp)
    :
    grids(grids)
{
    BL_ASSERT(ncomp > 0);
    BL_ASSERT(grids[0].cellCentered());

    for (OrientationIter face; face; ++face)
    {
        define(face(),IndexType::TheCellType(),in_rad,out_rad,extent_rad,ncomp);
    }
}

void
BndryRegister::init (const BndryRegister& src)
{
    grids.define(src.grids);

    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        bndry[i].define(src.bndry[i].boxArray(), src.bndry[i].nComp());

        for (FabSetIter mfi(src.bndry[i]); mfi.isValid(); ++mfi)
        {
            bndry[i][mfi.index()].copy(src.bndry[i][mfi]);
        }
    }
}

BndryRegister::BndryRegister (const BndryRegister& src)
{
    init(src);
}

const BoxArray&
BndryRegister::boxes () const
{
    return grids;
}

int
BndryRegister::size () const
{
    return grids.size();
}

const FabSet&
BndryRegister::operator[] (const Orientation& _face) const
{
    return bndry[_face];
}

FabSet&
BndryRegister::operator[] (const Orientation& _face)
{
    return bndry[_face];
}

BndryRegister&
BndryRegister::operator= (const BndryRegister& src)
{
    if (this != &src)
    {
        if (grids.size() > 0)
        {
            grids.clear();

            for (int i = 0; i < 2*BL_SPACEDIM; i++)
                bndry[i].clear();
        }

        init(src);
    }
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
    BL_ASSERT(grids.size() > 0);

    FabSet& fabs = bndry[_face];

    BL_ASSERT(fabs.size() == 0);

    const int coord_dir = _face.coordDir();
    const int lo_side   = _face.isLow();
    //
    // Build the BoxArray on which to define the FabSet on this face.
    //
    BoxArray fsBA(grids.size());

    for (int idx = 0; idx < grids.size(); ++idx)
    {
        Box b;
        //
        // First construct proper box for direction normal to face.
        //
        if (_out_rad > 0)
        {
            if (_typ.ixType(coord_dir) == IndexType::CELL)
                b = BoxLib::adjCell(grids[idx], _face, _out_rad);
            else
                b = BoxLib::bdryNode(grids[idx], _face, _out_rad);

            if (_in_rad > 0)
                b.grow(_face.flip(), _in_rad);
        }
        else
        {
            if (_in_rad > 0)
            {
                if (_typ.ixType(coord_dir) == IndexType::CELL)
                    b = BoxLib::adjCell(grids[idx], _face, _in_rad);
                else
                    b = BoxLib::bdryNode(grids[idx], _face, _in_rad);

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

        BL_ASSERT(b.ok());

        fsBA.set(idx,b);
    }

    BL_ASSERT(fsBA.ok());

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
    BL_ASSERT(grids.size() == 0);
    BL_ASSERT(_grids.size() > 0);
    BL_ASSERT(_grids[0].cellCentered());

    grids.define(_grids);
    //
    // Check that bndry regions are not allocated.
    //
    for (int k = 0; k < 2*BL_SPACEDIM; k++)
        BL_ASSERT(bndry[k].size() == 0);
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
