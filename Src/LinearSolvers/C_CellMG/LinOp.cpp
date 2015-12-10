
#include <winstd.H>
#include <cstdlib>

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <LO_BCTYPES.H>
#include <LO_F.H>
#include <LinOp.H>

namespace
{
    bool initialized = false;
}
//
// Set default values for these in Initialize()!!!
//
int LinOp::def_harmavg;
int LinOp::def_verbose;
int LinOp::def_maxorder;

#ifndef NDEBUG
//
// LinOp::applyBC fills LinOp_grow ghost cells with data expected in
// LinOp::apply therefore, the incoming MultiFab to LinOp::applyBC better
// have this many ghost allocated.
//
const int LinOp_grow = 1;
#endif

void
LinOp::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    LinOp::def_harmavg  = 0;
    LinOp::def_verbose  = 0;
    LinOp::def_maxorder = 2;

    ParmParse pp("Lp");

    pp.query("harmavg",  def_harmavg);
    pp.query("v",        def_verbose);
    pp.query("maxorder", def_maxorder);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        std::cout << "def_harmavg = "  << def_harmavg  << '\n';
        std::cout << "def_maxorder = " << def_maxorder << '\n';
    }

    BoxLib::ExecOnFinalize(LinOp::Finalize);

    initialized = true;
}

void
LinOp::Finalize ()
{
    ;
}

void
LinOp::bndryData (const BndryData& bd)
{
    BL_ASSERT(gbox[0] == bd.boxes());
    *bgb = bd;
}

LinOp::LinOp (const BndryData& _bgb,
              const Real       _h)
    :
    bgb(new BndryData(_bgb))
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
    bgb(new BndryData(_bgb))
{
    initConstruct(_h);
}

LinOp::LinOp (BndryData*  _bgb,
              const Real* _h)
    :
    bgb(_bgb)
{
    initConstruct(_h);
}

LinOp::~LinOp ()
{
    delete bgb;

    for (int i = 0, N = maskvals.size(); i < N; ++i)
    {
	for (int j = 0, M = maskvals[i].size(); j < M; ++j)
	{
            MaskTuple& a = maskvals[i].at_local(j);
            for (int k = 0; k < 2*BL_SPACEDIM; ++k)
                delete a[k];
        }
    }

    for (int i = 0, N = lmaskvals.size(); i < N; ++i)
    {
	for (int j = 0, M = lmaskvals[i].size(); j < M; ++j)
	{
            MaskTuple& a = lmaskvals[i].at_local(j);
            for (int k = 0; k < 2*BL_SPACEDIM; ++k)
                delete a[k];
        }
    }
}

void
LinOp::clearToLevel (int level) {}

void
LinOp::initConstruct (const Real* _h)
{
    Initialize();
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
    gbox[level] = bgb->boxes();
    geomarray.resize(1);
    geomarray[level] = bgb->getGeom();
    h.resize(1);
    maxorder = def_maxorder;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        h[level][i] = _h[i];
    }
    undrrelxr.resize(1);
    undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, 1);

     maskvals.resize(1);
    lmaskvals.resize(1);
    //
    // For each orientation, build NULL masks, then use distributed allocation.
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    maskvals[0].reserve((*bgb)[Orientation(0,Orientation::low)].local_size());
    lmaskvals[0].reserve((*bgb)[Orientation(0,Orientation::low)].local_size());
    for (FabSetIter bndryfsi((*bgb)[Orientation(0,Orientation::low)]);
         bndryfsi.isValid();
         ++bndryfsi)
    {
        const int        i   = bndryfsi.index();
        MaskTuple&       ma  =  maskvals[level][i];
        MaskTuple&       lma = lmaskvals[level][i];
        const MaskTuple& bdm = bgb->bndryMasks(i);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face = oitr();
            const Mask*       m    = bdm[face];
             ma[face] = new Mask(m->box(),1);
            lma[face] = new Mask(m->box(),1);
             ma[face]->copy(*m);
            lma[face]->copy(*m);
        }
    }
}

void
LinOp::apply (MultiFab&      out,
              MultiFab&      in,
              int            level,
              LinOp::BC_Mode bc_mode,
              bool           local,
	      int            src_comp,
	      int            dst_comp,
	      int            num_comp,
	      int            bndry_comp)
{
    applyBC(in,src_comp,num_comp,level,bc_mode,local,bndry_comp);
    Fapply(out,dst_comp,in,src_comp,num_comp,level);
}

void
LinOp::applyBC (MultiFab&      inout,
                int            src_comp,
                int            num_comp,
                int            level,
                LinOp::BC_Mode bc_mode,
                bool           local,
		int            bndry_comp)
{
    BL_PROFILE("LinOp::applyBC()");
    //
    // The inout MultiFab needs at least LinOp_grow ghost cells for applyBC.
    //
    BL_ASSERT(inout.nGrow() >= LinOp_grow);
    //
    // No coarsened boundary values, cannot apply inhomog at lev>0.
    //
    BL_ASSERT(level < numLevels());
    BL_ASSERT(!(level > 0 && bc_mode == Inhomogeneous_BC));

    int flagden = 1; // Fill in undrrelxr.
    int flagbc  = 1; // Fill boundary data.

    if (bc_mode == LinOp::Homogeneous_BC)
        //
        // No data if homogeneous.
        //
        flagbc = 0;

    const bool cross = true;

    inout.FillBoundary(src_comp,num_comp,local,cross);

    prepareForLevel(level);
    //
    // Do periodic fixup.
    //
    BL_ASSERT(level<geomarray.size());
    geomarray[level].FillPeriodicBoundary(inout,src_comp,num_comp,false,local);
    //
    // Fill boundary cells.
    //
    // OMP over boxes
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(inout); mfi.isValid(); ++mfi)
    {
        const int gn = mfi.index();

        BL_ASSERT(gbox[level][gn] == inout.box(gn));

        BL_ASSERT(level<maskvals.size() && maskvals[level].local_index(gn)>=0);
        BL_ASSERT(level<lmaskvals.size() && lmaskvals[level].local_index(gn)>=0);
        BL_ASSERT(level<undrrelxr.size());

        const MaskTuple&                 ma  =  maskvals[level][gn];
        const MaskTuple&                 lma = lmaskvals[level][gn];
        const BndryData::RealTuple&      bdl = bgb->bndryLocs(gn);
        const Array< Array<BoundCond> >& bdc = bgb->bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation o = oitr();

            FabSet&       f   = (*undrrelxr[level])[o];
            int           cdr = o;
            const FabSet& fs  = bgb->bndryValues(o);
            const Mask&   m   = local ? (*lma[o]) : (*ma[o]);
            Real          bcl = bdl[o];
            BL_ASSERT(bdc[o].size()>bndry_comp);
            int           bct = bdc[o][bndry_comp];

            const Box&       vbx   = inout.box(gn);
            FArrayBox&       iofab = inout[mfi];
            BL_ASSERT(f.size()>gn);
            BL_ASSERT(fs.size()>gn);

            FArrayBox&       ffab  = f[mfi];
            const FArrayBox& fsfab = fs[mfi];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         iofab.dataPtr(src_comp),
                         ARLIM(iofab.loVect()), ARLIM(iofab.hiVect()),
                         &cdr, &bct, &bcl,
                         fsfab.dataPtr(bndry_comp), 
                         ARLIM(fsfab.loVect()), ARLIM(fsfab.hiVect()),
                         m.dataPtr(),
                         ARLIM(m.loVect()), ARLIM(m.hiVect()),
                         ffab.dataPtr(),
                         ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                         vbx.loVect(),
                         vbx.hiVect(), &num_comp, h[level]);
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
    BL_PROFILE("LinOp::residual()");

    apply(residL, solnL, level, bc_mode, local);

    const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
    {
        const int nc = residL.nComp();
        //
        // Only single-component solves supported (verified) by this class.
        //
        BL_ASSERT(nc == 1);

        const Box& tbx = solnLmfi.tilebox();

        BL_ASSERT(gbox[level][solnLmfi.index()] == solnLmfi.validbox());

        FArrayBox& residfab = residL[solnLmfi];
        const FArrayBox& rhsfab = rhsL[solnLmfi];

        FORT_RESIDL(
            residfab.dataPtr(), 
            ARLIM(residfab.loVect()), ARLIM(residfab.hiVect()),
            rhsfab.dataPtr(), 
            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
            residfab.dataPtr(), 
            ARLIM(residfab.loVect()), ARLIM(residfab.hiVect()),
            tbx.loVect(), tbx.hiVect(), &nc);
    }
}

void
LinOp::smooth (MultiFab&       solnL,
               const MultiFab& rhsL,
               int             level,
               LinOp::BC_Mode  bc_mode)
{
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
    applyBC(solnL, 0, 1, level, bc_mode);
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
    BL_PROFILE("LinOp::prepareForLevel()");

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
    BL_ASSERT( maskvals.size() == level);
    BL_ASSERT(lmaskvals.size() == level);
     maskvals.resize(level+1);
    lmaskvals.resize(level+1);

    Array<IntVect> pshifts(27);

    std::vector< std::pair<int,Box> > isects;
    //
    // Use bgb's distribution map for masks.
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    maskvals[level].reserve((*bgb)[Orientation(0,Orientation::low)].local_size());
    lmaskvals[level].reserve((*bgb)[Orientation(0,Orientation::low)].local_size());
    for (FabSetIter bndryfsi((*bgb)[Orientation(0,Orientation::low)]);
         bndryfsi.isValid();
         ++bndryfsi)
    {
        const int   gn = bndryfsi.index();
        MaskTuple&  ma =  maskvals[level][gn];
        MaskTuple& lma = lmaskvals[level][gn];

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face = oitr();
            const Box&        bx_k = BoxLib::adjCell(gbox[level][gn], face, 1);
             ma[face] = new Mask(bx_k,1);
            lma[face] = new Mask(bx_k,1);
            Mask&  curmask = *( ma[face]);
            Mask& lcurmask = *(lma[face]);
             curmask.setVal(BndryData::not_covered);
            lcurmask.setVal(BndryData::not_covered);
            gbox[level].intersections(bx_k,isects);
            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                if (isects[ii].first != gn)
                {
                     curmask.setVal(BndryData::covered, isects[ii].second, 0);
                    lcurmask.setVal(BndryData::covered, isects[ii].second, 0);
                }
            }
            //
            // Now take care of periodic wraparounds.
            //
            Geometry& curgeom = geomarray[level];

            if (curgeom.isAnyPeriodic() && !curdomain.contains(bx_k))
            {
                curgeom.periodicShift(curdomain, bx_k, pshifts);

                for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
                {
                     curmask.shift(pshifts[iiv]);
                    lcurmask.shift(pshifts[iiv]);
                    BL_ASSERT(curmask.box() == lcurmask.box());
                    gbox[level].intersections(curmask.box(),isects);
                    for (int ii = 0, N = isects.size(); ii < N; ii++)
                    {
                         curmask.setVal(BndryData::covered, isects[ii].second, 0);
                        lcurmask.setVal(BndryData::covered, isects[ii].second, 0);
                    }
                     curmask.shift(-pshifts[iiv]);
                    lcurmask.shift(-pshifts[iiv]);
                }
            }
        }
    }

    gbox[level].clear_hash_bin();
}

void
LinOp::makeCoefficients (MultiFab&       cs,
                         const MultiFab& fn,
                         int             level)
{
    BL_PROFILE("LinOp::makeCoefficients()");

    int nc = 1;
    //
    // Determine index type of incoming MultiFab.
    //
    const IndexType iType(fn.boxArray().ixType());
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

    const bool tiling = true;

    switch (cdir)
    {
    case -1:
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter csmfi(cs,tiling); csmfi.isValid(); ++csmfi)
        {
            const Box& tbx = csmfi.tilebox();
            FArrayBox&       csfab = cs[csmfi];
            const FArrayBox& fnfab = fn[csmfi];

            FORT_AVERAGECC(csfab.dataPtr(), ARLIM(csfab.loVect()),
                           ARLIM(csfab.hiVect()),fnfab.dataPtr(),
                           ARLIM(fnfab.loVect()),ARLIM(fnfab.hiVect()),
                           tbx.loVect(),tbx.hiVect(), &nc);
        }
        break;
    case 0:
    case 1:
    case 2:
        if (harmavg)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
  	    for (MFIter csmfi(cs,tiling); csmfi.isValid(); ++csmfi)
            {
	        const Box& tbx = csmfi.tilebox();
                FArrayBox&       csfab = cs[csmfi];
                const FArrayBox& fnfab = fn[csmfi];

                FORT_HARMONIC_AVERAGEEC(csfab.dataPtr(),
                                        ARLIM(csfab.loVect()),
                                        ARLIM(csfab.hiVect()),
                                        fnfab.dataPtr(),
                                        ARLIM(fnfab.loVect()),
                                        ARLIM(fnfab.hiVect()),
                                        tbx.loVect(),tbx.hiVect(),
                                        &nc,&cdir);
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter csmfi(cs,tiling); csmfi.isValid(); ++csmfi)
            {
                const Box& tbx = csmfi.tilebox();
                FArrayBox&       csfab = cs[csmfi];
                const FArrayBox& fnfab = fn[csmfi];

                FORT_AVERAGEEC(csfab.dataPtr(),ARLIM(csfab.loVect()),
                               ARLIM(csfab.hiVect()),fnfab.dataPtr(), 
                               ARLIM(fnfab.loVect()),ARLIM(fnfab.hiVect()),
	                       tbx.loVect(),tbx.hiVect(),
                               &nc, &cdir);
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

        const BLMap<LinOp::MaskTuple>& m = lp.maskvals[level];

        for (int nproc = 0; nproc < ParallelDescriptor::NProcs(); ++nproc)
        {
            if (nproc == ParallelDescriptor::MyProc())
            {
                os << "Processor " << nproc << '\n';

                for (OrientationIter oitr; oitr; ++oitr)
                {
                    const Orientation face = oitr();

		    for (int i=0, N=m.size(); i<N; ++i)
                    {
                        os << m.at_local(i)[face];
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

int
LinOp::maxOrder (int maxorder_)
{
    BL_ASSERT(maxorder_ >= 2);
    maxorder_ = (maxorder_ < 2 ? 2 : maxorder_ );
    int omaxorder = maxorder;
    maxorder = maxorder_;
    return omaxorder;
}
