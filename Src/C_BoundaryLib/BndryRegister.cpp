
#include <winstd.H>
#include <BndryRegister.H>
#include <Orientation.H>
#include <Utility.H>

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
            bndry[i][mfi].copy(src.bndry[i][mfi]);
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
BndryRegister::defineDoit (Orientation _face,
                           IndexType   _typ,
                           int         _in_rad,
                           int         _out_rad,
                           int         _extent_rad,
                           BoxArray&   fsBA)
{
    BL_PROFILE("BndryRegister::defineDoit()");

    BL_ASSERT(grids.size() > 0);

    const int coord_dir = _face.coordDir();
    const int lo_side   = _face.isLow();
    //
    // Build the BoxArray on which to define the FabSet on this face.
    //
    const int N = grids.size();

    fsBA.resize(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int idx = 0; idx < N; ++idx)
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
}

void
BndryRegister::define (Orientation _face,
                       IndexType   _typ,
                       int         _in_rad,
                       int         _out_rad,
                       int         _extent_rad,
                       int         _ncomp)
{
    BoxArray fsBA;

    defineDoit(_face,_typ,_in_rad,_out_rad,_extent_rad,fsBA);

    FabSet& fabs = bndry[_face];

    BL_ASSERT(fabs.size() == 0);

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
BndryRegister::define (Orientation                _face,
                       IndexType                  _typ,
                       int                        _in_rad,
                       int                        _out_rad,
                       int                        _extent_rad,
                       int                        _ncomp,
                       const DistributionMapping& _dm)
{
    BoxArray fsBA;

    defineDoit(_face,_typ,_in_rad,_out_rad,_extent_rad,fsBA);

    FabSet& fabs = bndry[_face];

    BL_ASSERT(fabs.size() == 0);

    fabs.define(fsBA,_ncomp,_dm);
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
BndryRegister::operator+= (const BndryRegister& rhs)
{
    BL_ASSERT(grids == rhs.grids);
    for (OrientationIter face; face; ++face) {
	for (FabSetIter bfsi(rhs[face()]); bfsi.isValid(); ++bfsi) {
	    bndry[face()][bfsi] += rhs[face()][bfsi];
	}
    }
    return *this;
}

BndryRegister&
BndryRegister::plus (const BndryRegister& rhs)
{
    return operator+=(rhs);
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

void
BndryRegister::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        grids.writeOn(os);
        os << '\n';
    }

    for (OrientationIter face; face; ++face)
    {
        //
        // Take name here and make a "new" name unique to each face.
        // Simplest thing would probably to append "_n" to the name,
        // where n is the integer value of face().
        //
        const int i = face();
        BL_ASSERT(i >= 0 && i <= 7);

        std::string facename = BoxLib::Concatenate(name + '_', i, 1);

        bndry[face()].write(facename);
    }
}

void
BndryRegister::read (const std::string& name, std::istream& is)
{
    grids.readFrom(is);

    for (OrientationIter face; face; ++face)
    {
        //
        // Take name here and make a "new" name unique to each face.
        // Simplest thing would probably to append "_n" to the name,
        // where n is the integer value of face().
        //
        const int i = face();
        BL_ASSERT(i >= 0 && i <= 7);

        std::string facename = BoxLib::Concatenate(name + '_', i, 1);

        bndry[face()].read(facename);
    }
}

void
BndryRegister::MakeSidecarsSmaller(int ioProcNumSCS, int ioProcNumAll,
                                   int scsMyId, MPI_Comm scsComm)
{
  // ---- BoxArrays
  BoxLib::BroadcastBoxArray(grids, scsMyId, ioProcNumSCS, scsComm);

  // ---- FabSet
  for(int i(0); i < (2 * BL_SPACEDIM); ++i) {
    bndry[i].MakeSidecarsSmaller(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
  }
}



