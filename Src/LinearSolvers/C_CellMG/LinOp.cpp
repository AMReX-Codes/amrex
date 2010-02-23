
//
// $Id: LinOp.cpp,v 1.39 2010-02-23 21:08:02 lijewski Exp $
//
#include <winstd.H>

#include <cstdlib>

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <LO_BCTYPES.H>
#include <LO_F.H>
#include <LinOp.H>
#include <Profiler.H>



bool LinOp::initialized = false;
int LinOp::def_harmavg  = 0;
int LinOp::def_verbose  = 0;
int LinOp::def_maxorder = 2;

#ifndef NDEBUG
//
// LinOp::applyBC fills LinOp_grow ghost cells with data expected in
// LinOp::apply therefore, the incoming MultiFab to LinOp::applyBC better
// have this many ghost allocated.
//
const int LinOp_grow = 1;
#endif

void
LinOp::initialize ()
{
    ParmParse pp("Lp");

    pp.query("harmavg", def_harmavg);
    pp.query("v", def_verbose);
    pp.query("maxorder", def_maxorder);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
      {
        std::cout << "def_harmavg = " << def_harmavg << '\n';
        std::cout << "def_maxorder = " << def_maxorder << '\n';
      }

    initialized = true;
}

LinOp::LinOp (const BndryData& _bgb,
              const Real       _h)
    :
    bgb(_bgb)
{
    Real __h[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        __h[i] = _h;
    }
    initConstruct(__h);
}

LinOp::LinOp (const BndryData& _bgb,
              const Real*      _h)
    :
    bgb(_bgb)
{
    initConstruct(_h);
}

LinOp::~LinOp ()
{


    for (int i = 0; i < maskvals.size(); ++i)
    {
        for (int j = 0; j < maskvals[i].size(); ++j)
        {
            for (int k = 0; k < maskvals[i][j].size(); ++k)
            {
                delete maskvals[i][j][k];
            }
        }
    }

    for (int i = 0; i < lmaskvals.size(); ++i)
    {
        for (int j = 0; j < lmaskvals[i].size(); ++j)
        {
            for (int k = 0; k < lmaskvals[i][j].size(); ++k)
            {
                delete lmaskvals[i][j][k];
            }
        }
    }

}

void
LinOp::initConstruct (const Real* _h)
{
    if (!initialized)
        initialize();
    //
    // We'll reserve() space to cut down on copying during resize()s.
    //
    const int N = 10;

    h.reserve(N);
    gbox.reserve(N);
    undrrelxr.reserve(N);
    maskvals.reserve(N);
    lmaskvals.reserve(N);
    geomarray.reserve(N);

    harmavg = def_harmavg;
    verbose = def_verbose;
    gbox.resize(1);
    const int level = 0;
    gbox[level] = bgb.boxes();
    geomarray.resize(1);
    geomarray[level] = bgb.getGeom();
    h.resize(1);
    maxorder = def_maxorder;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        h[level][i] = _h[i];
    }
    undrrelxr.resize(1);
    undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, 1);

    maskvals.resize(1);
    maskvals[level].resize(gbox[level].size());
    //
    // For each orientation, build NULL masks, then use distributed allocation.
    //
    for (int i = 0; i < gbox[level].size(); i++)
    {
        maskvals[level][i].resize(2*BL_SPACEDIM, 0);
    }

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face = oitr();
        const FabSet& bndry = bgb[face];
        for (int i = 0; i < gbox[level].size(); i++)
        {
            if (bndry.DistributionMap()[i] == ParallelDescriptor::MyProc())
            {
                const PArray<Mask>& pam = bgb.bndryMasks(face);
                BL_ASSERT(maskvals[level][i][face] == 0);
                maskvals[level][i][face] = new Mask(pam[i].box(), 1);
                maskvals[level][i][face]->copy(pam[i]);
            }
        }
    }

    lmaskvals.resize(1);
    lmaskvals[level].resize(gbox[level].size());
    //
    // For each orientation, build NULL masks, then use distributed allocation.
    //
    for (int i = 0; i < gbox[level].size(); i++)
    {
        lmaskvals[level][i].resize(2*BL_SPACEDIM, 0);
    }

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face = oitr();
        const FabSet& bndry = bgb[face];
        for (int i = 0; i < gbox[level].size(); i++)
        {
            if (bndry.DistributionMap()[i] == ParallelDescriptor::MyProc())
            {
                const PArray<Mask>& pam = bgb.bndryMasks(face);
                BL_ASSERT(lmaskvals[level][i][face] == 0);
                lmaskvals[level][i][face] = new Mask(pam[i].box(), 1);
                lmaskvals[level][i][face]->copy(pam[i]);
            }
        }
    }
}

void
LinOp::apply (MultiFab&      out,
              MultiFab&      in,
              int            level,
              LinOp::BC_Mode bc_mode,
              bool           local)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "apply()");

    applyBC(in,0,1,level,bc_mode,local);
    Fapply(out,in,level);
}

void
LinOp::applyBC (MultiFab&      inout,
                int            src_comp,
                int            num_comp,
                int            level,
                LinOp::BC_Mode bc_mode,
                bool           local)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::applyBC()");
    //
    // The inout MultiFab needs at least LinOp_grow ghost cells for applyBC.
    //
    BL_ASSERT(inout.nGrow() >= LinOp_grow);
    //
    // The inout MultiFab must have at least Periodic_BC_grow cells for the
    // algorithms taking care of periodic boundary conditions.
    //
    BL_ASSERT(inout.nGrow() >= LinOp_grow);
    //
    // No coarsened boundary values, cannot apply inhomog at lev>0.
    //
    BL_ASSERT(!(level > 0 && bc_mode == Inhomogeneous_BC));

    int flagden = 1; // Fill in undrrelxr.
    int flagbc  = 1; // Fill boundary data.

    if (bc_mode == LinOp::Homogeneous_BC)
        //
        // No data if homogeneous.
        //
        flagbc = 0;
    //
    // Only single-component solves supported (verified) by this class.
    //
    BL_ASSERT(num_comp == 1);

    inout.FillBoundary(src_comp,num_comp,local);

    prepareForLevel(level);
    //
    // Do periodic fixup.
    //
    geomarray[level].FillPeriodicBoundary(inout,src_comp,num_comp,false,local);
    //
    // Fill boundary cells.
    //
    const int comp = 0;

    const int N = inout.IndexMap().size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
    {
        const int gn = inout.IndexMap()[i];

        BL_ASSERT(gbox[level][gn] == inout.box(gn));

        for (OrientationIter oitr; oitr; ++oitr)
        {
            Orientation o(oitr());

            const Array< Array<BoundCond> >& b   = bgb.bndryConds(o);
            const Array<Real>&               r   = bgb.bndryLocs(o);
            FabSet&                          f   = (*undrrelxr[level])[o];
            int                              cdr = o;
            const FabSet&                    fs  = bgb.bndryValues(o);
            const Mask&                      m   = local ? (*lmaskvals[level][gn][o]) : (*maskvals[level][gn][o]);
            Real                             bcl = r[gn];
            int                              bct = b[gn][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         inout[gn].dataPtr(src_comp),
                         ARLIM(inout[gn].loVect()), ARLIM(inout[gn].hiVect()),
                         &cdr, &bct, &bcl,
                         fs[gn].dataPtr(), 
                         ARLIM(fs[gn].loVect()), ARLIM(fs[gn].hiVect()),
                         m.dataPtr(),
                         ARLIM(m.loVect()), ARLIM(m.hiVect()),
                         f[gn].dataPtr(),
                         ARLIM(f[gn].loVect()), ARLIM(f[gn].hiVect()),
                         inout.box(gn).loVect(),
                         inout.box(gn).hiVect(), &num_comp, h[level]);
        }
    }
}

void
LinOp::residual (MultiFab&       residL,
                 const MultiFab& rhsL,
                 MultiFab&       solnL,
                 int             level,
                 LinOp::BC_Mode  bc_mode,
                 bool            local)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "residual()");

    apply(residL, solnL, level, bc_mode, local);

    for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
        int nc = residL.nComp();
        //
        // Only single-component solves supported (verified) by this class.
        //
        BL_ASSERT(nc == 1);
        BL_ASSERT(gbox[level][solnLmfi.index()] == solnLmfi.validbox());
        FORT_RESIDL(
            residL[solnLmfi].dataPtr(), 
            ARLIM(residL[solnLmfi].loVect()), ARLIM(residL[solnLmfi].hiVect()),
            rhsL[solnLmfi].dataPtr(), 
            ARLIM(rhsL[solnLmfi].loVect()), ARLIM(rhsL[solnLmfi].hiVect()),
            residL[solnLmfi].dataPtr(), 
            ARLIM(residL[solnLmfi].loVect()), ARLIM(residL[solnLmfi].hiVect()),
            solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(), &nc);
    }
}

void
LinOp::smooth (MultiFab&       solnL,
               const MultiFab& rhsL,
               int             level,
               LinOp::BC_Mode  bc_mode)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::smooth()");

    for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++)
    {
        applyBC(solnL, 0, 1, level, bc_mode);
        Fsmooth(solnL, rhsL, level, redBlackFlag);
    }
}

void
LinOp::jacobi_smooth (MultiFab&       solnL,
                      const MultiFab& rhsL,
                      int             level,
                      LinOp::BC_Mode  bc_mode)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::jacobi_smooth()");
    Fsmooth_jacobi(solnL, rhsL, level);
}

Real
LinOp::norm (int nm, int level, const bool local)
{
    BoxLib::Error("LinOp::norm: Placeholder for pure virtual function");
    return 0;
}

void
LinOp::prepareForLevel (int level)
{
    if (level == 0) return;

    LinOp::prepareForLevel(level-1);

    if (h.size() > level) return;
    //
    // Assume from here down that this is a new level one coarser than existing
    //
    BL_ASSERT(h.size() == level);
    h.resize(level+1);
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        h[level][i] = h[level-1][i]*2.0;
    }
    geomarray.resize(level+1);
    geomarray[level].define(BoxLib::coarsen(geomarray[level-1].Domain(),2));
    const Box& curdomain = geomarray[level].Domain();
    //
    // Add a box to the new coarser level (assign removes old BoxArray).
    //
    gbox.resize(level+1);
    gbox[level] = gbox[level-1];
    gbox[level].coarsen(2);
    //
    // Add the BndryRegister of relax values to the new coarser level.
    //
    BL_ASSERT(undrrelxr.size() == level);
    undrrelxr.resize(level+1);
    undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, 1);
    //
    // Add an Array of Array of maskvals to the new coarser level
    // For each orientation, build NULL masks, then use distributed allocation
    // Initial masks for coarse levels, ignore outside_domain possibility since
    // we always solve homogeneous equation on coarse levels.
    //
    BL_ASSERT(maskvals.size() == level);
    maskvals.resize(level+1);
    maskvals[level].resize(gbox[level].size());
    for (int i = 0; i < gbox[level].size(); i++)
        maskvals[level][i].resize(2*BL_SPACEDIM, 0);

    Array<IntVect> pshifts(27);

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face = oitr();
        //
        // Use bgb's distribution map for masks.
        //
        for (FabSetIter bndryfsi(bgb[face]); bndryfsi.isValid(); ++bndryfsi)
        {
            const int gn = bndryfsi.index();
            Box bx_k = BoxLib::adjCell(gbox[level][gn], face, 1);
            BL_ASSERT(maskvals[level][gn][face] == 0);
            maskvals[level][gn][face] = new Mask(bx_k, 1);
            Mask& curmask = *(maskvals[level][gn][face]);
            curmask.setVal(BndryData::not_covered);
            std::vector< std::pair<int,Box> > isects = gbox[level].intersections(bx_k);
            for (int ii = 0; ii < isects.size(); ii++)
                if (isects[ii].first != gn)
                    curmask.setVal(BndryData::covered, isects[ii].second, 0);
            //
            // Now take care of periodic wraparounds.
            //
            Geometry& curgeom = geomarray[level];

            if (curgeom.isAnyPeriodic() && !curdomain.contains(bx_k))
            {
                curgeom.periodicShift(curdomain, bx_k, pshifts);

                for (int iiv = 0; iiv < pshifts.size(); iiv++)
                {
                    curmask.shift(pshifts[iiv]);
                    std::vector< std::pair<int,Box> > isects = gbox[level].intersections(curmask.box());
                    for (int ii = 0; ii < isects.size(); ii++)
                        curmask.setVal(BndryData::covered, isects[ii].second, 0);
                    curmask.shift(-pshifts[iiv]);
                }
            }
        }
    }
    //
    // Ditto for lmaskvals
    //
    BL_ASSERT(lmaskvals.size() == level);
    lmaskvals.resize(level+1);
    lmaskvals[level].resize(gbox[level].size());
    for (int i = 0; i < gbox[level].size(); i++)
        lmaskvals[level][i].resize(2*BL_SPACEDIM, 0);

    const int MyProc = ParallelDescriptor::MyProc();

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face = oitr();
        //
        // Use bgb's distribution map for masks.
        //
        for (FabSetIter bndryfsi(bgb[face]); bndryfsi.isValid(); ++bndryfsi)
        {
            const int gn = bndryfsi.index();
            Box bx_k = BoxLib::adjCell(gbox[level][gn], face, 1);
            BL_ASSERT(lmaskvals[level][gn][face] == 0);
            lmaskvals[level][gn][face] = new Mask(bx_k, 1);
            Mask& curmask = *(lmaskvals[level][gn][face]);
            curmask.setVal(BndryData::not_covered);
            std::vector< std::pair<int,Box> > isects = gbox[level].intersections(bx_k);
            for (int ii = 0; ii < isects.size(); ii++)
                if (isects[ii].first != gn && bgb[face].DistributionMap()[isects[ii].first] == MyProc)
                    curmask.setVal(BndryData::covered, isects[ii].second, 0);
            //
            // Now take care of periodic wraparounds.
            //
            Geometry& curgeom = geomarray[level];

            if (curgeom.isAnyPeriodic() && !curdomain.contains(bx_k))
            {
                curgeom.periodicShift(curdomain, bx_k, pshifts);

                for (int iiv = 0; iiv < pshifts.size(); iiv++)
                {
                    curmask.shift(pshifts[iiv]);
                    std::vector< std::pair<int,Box> > isects = gbox[level].intersections(curmask.box());
                    for (int ii = 0; ii < isects.size(); ii++)
                        if (bgb[face].DistributionMap()[isects[ii].first] == MyProc)
                            curmask.setVal(BndryData::covered, isects[ii].second, 0);
                    curmask.shift(-pshifts[iiv]);
                }
            }
        }
    }
}

void
LinOp::makeCoefficients (MultiFab&       cs,
                         const MultiFab& fn,
                         int             level)
{
    int nc = 1;
    //
    // Determine index type of incoming MultiFab.
    //
    const IndexType iType(fn.boxArray()[0].ixType());
    const IndexType cType(D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL));
    const IndexType xType(D_DECL(IndexType::NODE, IndexType::CELL, IndexType::CELL));
    const IndexType yType(D_DECL(IndexType::CELL, IndexType::NODE, IndexType::CELL));
#if (BL_SPACEDIM == 3)    
    const IndexType zType(D_DECL(IndexType::CELL, IndexType::CELL, IndexType::NODE));
#endif

    int cdir;
    if (iType == cType)
    {
        cdir = -1;
    }
    else if (iType == xType)
    {
        cdir = 0;
    }
    else if (iType == yType)
    {
        cdir = 1;
#if (BL_SPACEDIM == 3)
    }
    else if (iType == zType)
    {
        cdir = 2;
#endif    
    }
    else
    {
        BoxLib::Error("LinOp::makeCoeffients: Bad index type");
    }

    BoxArray d(gbox[level]);
    if (cdir >= 0)
        d.surroundingNodes(cdir);
    //
    // Only single-component solves supported (verified) by this class.
    //
    const int nComp=1;
    const int nGrow=0;
    cs.define(d, nComp, nGrow, Fab_allocate);

    const BoxArray& grids = gbox[level];

    MFIter csmfi(cs);

    switch (cdir)
    {
    case -1:
        for ( ; csmfi.isValid(); ++csmfi)
        {
            FORT_AVERAGECC(cs[csmfi].dataPtr(), ARLIM(cs[csmfi].loVect()),
                           ARLIM(cs[csmfi].hiVect()),fn[csmfi].dataPtr(),
                           ARLIM(fn[csmfi].loVect()),ARLIM(fn[csmfi].hiVect()),
                           grids[csmfi.index()].loVect(),
                           grids[csmfi.index()].hiVect(), &nc);
        }
        break;
    case 0:
    case 1:
    case 2:
        if (harmavg)
        {
            for ( ; csmfi.isValid(); ++csmfi)
            {
                FORT_HARMONIC_AVERAGEEC(cs[csmfi].dataPtr(),
                                        ARLIM(cs[csmfi].loVect()),
                                        ARLIM(cs[csmfi].hiVect()),
                                        fn[csmfi].dataPtr(),
                                        ARLIM(fn[csmfi].loVect()),
                                        ARLIM(fn[csmfi].hiVect()),
                                        grids[csmfi.index()].loVect(),
                                        grids[csmfi.index()].hiVect(),
                                        &nc,&cdir);
            }
        }
        else
        {
            for ( ; csmfi.isValid(); ++csmfi)
            {
                FORT_AVERAGEEC(cs[csmfi].dataPtr(),ARLIM(cs[csmfi].loVect()),
                               ARLIM(cs[csmfi].hiVect()),fn[csmfi].dataPtr(), 
                               ARLIM(fn[csmfi].loVect()),ARLIM(fn[csmfi].hiVect()),
                               grids[csmfi.index()].loVect(),
                               grids[csmfi.index()].hiVect(),&nc, &cdir);
            }
        }
        break;
    default:
        BoxLib::Error("LinOp:: bad coefficient coarsening direction!");
    }
}

std::ostream&
operator<< (std::ostream& os,
            const LinOp&  lp)
{
    if (ParallelDescriptor::IOProcessor())
    {
        os << "LinOp" << '\n';
        os << "Grids: " << '\n';
        for (int level = 0; level < lp.h.size(); ++level)
        {
            os << " level = " << level << ": " << lp.gbox[level] << '\n';
        }
        os << "Grid Spacing: " << '\n';
        for (int level = 0; level < lp.h.size(); ++level)
        {
            os << " level = " << level << ", dx = ";
            for (int d =0; d < BL_SPACEDIM; ++d)
            {
                os << lp.h[level][d] << "  ";
            }
            os << '\n';
        }
        os << "Harmonic average? " << (lp.harmavg == 1 ? "yes" : "no") << '\n';
        os << "Verbosity: " << lp.verbose << '\n';
        os << "Max Order: " << lp.maxorder << '\n';
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << "Masks:" << '\n';
    }
    for (int level = 0; level < lp.h.size(); ++level)
    {
        if (ParallelDescriptor::IOProcessor())
            os << "level = " << level << '\n';

        for (int nproc = 0; nproc < ParallelDescriptor::NProcs(); ++nproc)
        {
            if (nproc == ParallelDescriptor::MyProc())
            {
                os << "Processor " << nproc << '\n';
                for (OrientationIter oitr; oitr; ++oitr)
                {
                    Orientation face = oitr();
                    for (int i=0; i<lp.boxArray().size(); ++i)
                    {
                        if (lp.maskvals[level][i][face])
                        {
                            os << *lp.maskvals[level][i][face];
                        }
                    }
                }
            }
        }
    }
    
    return os;
}

const Geometry&
LinOp::getGeom (int level)
{
    return geomarray[level];
}

const Real * 
LinOp::getDx (int level)
{
    return h[level];
}

Real
LinOp::get_alpha () const
{
    BoxLib::Abort("LinOp::get_alpha");
    return 0;
}   
    
Real
LinOp::get_beta () const 
{   
    BoxLib::Abort("LinOp::get_beta");
    return 0; 
}

const MultiFab&
LinOp::aCoefficients (int level)
{
    static MultiFab junk;
    BoxLib::Abort("LinOp::aCoefficients");
    return junk;
}

const MultiFab&
LinOp::bCoefficients (int dir,int level)
{
    static MultiFab junk;
    BoxLib::Abort("LinOp::bCoefficients");
    return junk;
}

