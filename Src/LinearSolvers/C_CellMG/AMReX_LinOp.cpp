
#include <cstdlib>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LO_F.H>
#include <AMReX_LinOp.H>

namespace amrex {

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
int LinOp::LinOp_grow;

// Important:
// LinOp::applyBC fills LinOp_grow ghost cells with data expected in
// LinOp::apply therefore, the incoming MultiFab to LinOp::applyBC better
// have LinOp_grow many ghost allocated.
//

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
    LinOp::LinOp_grow   = 1; // Must be consistent with expectations of apply/applyBC, not parm-parsed

    ParmParse pp("Lp");

    pp.query("harmavg",  def_harmavg);
    pp.query("v",        def_verbose);
    pp.query("maxorder", def_maxorder);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        amrex::Print() << "def_harmavg = "  << def_harmavg  << '\n'
                       << "def_maxorder = " << def_maxorder << '\n';
    }

    amrex::ExecOnFinalize(LinOp::Finalize);

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
    Real _hh[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        _hh[i] = _h;
    }
    initConstruct(_hh);
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
    undrrelxr[0].define(gbox[level], bgb->DistributionMap(), 1, 0, 0, 1);

    maskvals.resize(1);
    maskvals[0].resize(2*BL_SPACEDIM);

    lmaskvals.resize(1);
    lmaskvals[0].resize(2*BL_SPACEDIM);

    for (OrientationIter oitr; oitr; ++oitr)
    {
	const Orientation face = oitr();
	const MultiMask& m = bgb->bndryMasks(face);
	maskvals[0][face].define(m.boxArray(), m.DistributionMap(), 1);
	lmaskvals[0][face].define(m.boxArray(), m.DistributionMap(), 1);
	MultiMask::Copy(maskvals[0][face], m);
	MultiMask::Copy(lmaskvals[0][face], m);
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

    prepareForLevel(level);

    const bool cross = true;
    inout.FillBoundary(src_comp,num_comp,geomarray[level].periodicity(),cross);

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
        BL_ASSERT(level<undrrelxr.size());

        const BndryData::RealTuple&      bdl = bgb->bndryLocs(gn);
        const Vector< Vector<BoundCond> >& bdc = bgb->bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation o = oitr();

            FabSet&       f   = undrrelxr[level][o];
            int           cdr = o;
            const FabSet& fs  = bgb->bndryValues(o);
            const Mask&   m   = local ? lmaskvals[level][o][mfi] : maskvals[level][o][mfi];
            Real          bcl = bdl[o];
            BL_ASSERT(bdc[o].size()>bndry_comp);
            int           bct = bdc[o][bndry_comp];

            const Box&       vbx   = inout.box(gn);
            FArrayBox&       iofab = inout[mfi];
            BL_ASSERT(f.size()>gn);
            BL_ASSERT(fs.size()>gn);

            FArrayBox&       ffab  = f[mfi];
            const FArrayBox& fsfab = fs[mfi];

            amrex_lo_applybc(&flagden, &flagbc, &maxorder,
                         iofab.dataPtr(src_comp),
                         AMREX_ARLIM(iofab.loVect()), AMREX_ARLIM(iofab.hiVect()),
                         &cdr, &bct, &bcl,
                         fsfab.dataPtr(bndry_comp), 
                         AMREX_ARLIM(fsfab.loVect()), AMREX_ARLIM(fsfab.hiVect()),
                         m.dataPtr(),
                         AMREX_ARLIM(m.loVect()), AMREX_ARLIM(m.hiVect()),
                         ffab.dataPtr(),
                         AMREX_ARLIM(ffab.loVect()), AMREX_ARLIM(ffab.hiVect()),
                         vbx.loVect(),
                         vbx.hiVect(), &num_comp, h[level].data());
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
    MultiFab::Xpay(residL, -1.0, rhsL, 0, 0, residL.nComp(), 0);
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
    amrex::Error("LinOp::norm: Placeholder for pure virtual function");
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
    geomarray[level].define(amrex::coarsen(geomarray[level-1].Domain(),2));
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
    undrrelxr[level].define(gbox[level], bgb->DistributionMap(), 1, 0, 0, 1);
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

    maskvals[level].resize(2*BL_SPACEDIM);
    lmaskvals[level].resize(2*BL_SPACEDIM);

    int nGrow = NumGrow(level);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();
	maskvals[level][face].define(gbox[level], bgb->DistributionMap(), geomarray[level],
				     face, 0, nGrow, 0, 1, true);
	lmaskvals[level][face].define(maskvals[level][face].boxArray(),
				      maskvals[level][face].DistributionMap(), 1);
	MultiMask::Copy(lmaskvals[level][face], maskvals[level][face]);
    }
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
    const IndexType cType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL));
    const IndexType xType(AMREX_D_DECL(IndexType::NODE, IndexType::CELL, IndexType::CELL));
    const IndexType yType(AMREX_D_DECL(IndexType::CELL, IndexType::NODE, IndexType::CELL));
#if (BL_SPACEDIM == 3)    
    const IndexType zType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::NODE));
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
        cdir = -100;
        amrex::Error("LinOp::makeCoeffients: Bad index type");
    }

    BoxArray d(gbox[level]);
    if (cdir >= 0)
        d.surroundingNodes(cdir);
    //
    // Only single-component solves supported (verified) by this class.
    //
    const int nComp=1;
    const int nGrow=0;
    cs.define(d, fn.DistributionMap(), nComp, nGrow, MFInfo(), FArrayBoxFactory());

    const bool tiling = true;

    switch (cdir)
    {
    case -1:
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter csmfi(cs,tiling); csmfi.isValid(); ++csmfi)
        {
            const Box& tbx = csmfi.tilebox();
            FArrayBox&       csfab = cs[csmfi];
            const FArrayBox& fnfab = fn[csmfi];

            amrex_lo_averagecc(csfab.dataPtr(), AMREX_ARLIM(csfab.loVect()),
                           AMREX_ARLIM(csfab.hiVect()),fnfab.dataPtr(),
                           AMREX_ARLIM(fnfab.loVect()),AMREX_ARLIM(fnfab.hiVect()),
                           tbx.loVect(),tbx.hiVect(), &nc);
        }
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

                amrex_lo_harmonic_averageec(csfab.dataPtr(),
                                        AMREX_ARLIM(csfab.loVect()),
                                        AMREX_ARLIM(csfab.hiVect()),
                                        fnfab.dataPtr(),
                                        AMREX_ARLIM(fnfab.loVect()),
                                        AMREX_ARLIM(fnfab.hiVect()),
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

                amrex_lo_averageec(csfab.dataPtr(),AMREX_ARLIM(csfab.loVect()),
                               AMREX_ARLIM(csfab.hiVect()),fnfab.dataPtr(), 
                               AMREX_ARLIM(fnfab.loVect()),AMREX_ARLIM(fnfab.hiVect()),
	                       tbx.loVect(),tbx.hiVect(),
                               &nc, &cdir);
            }
        }
        break;
    default:
        amrex::Error("LinOp:: bad coefficient coarsening direction!");
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
                    const Orientation face = oitr();

		    for (MultiMaskIter mmi(lp.maskvals[level][face]); mmi.isValid(); ++mmi)
		    {
                        os << lp.maskvals[level][face][mmi];
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
    return h[level].data();
}

Real
LinOp::get_alpha () const
{
    amrex::Abort("LinOp::get_alpha");
    return 0;
}   
    
Real
LinOp::get_beta () const 
{   
    amrex::Abort("LinOp::get_beta");
    return 0; 
}

const MultiFab&
LinOp::aCoefficients (int level)
{
    static MultiFab junk;
    amrex::Abort("LinOp::aCoefficients");
    return junk;
}

const MultiFab&
LinOp::bCoefficients (int dir,int level)
{
    static MultiFab junk;
    amrex::Abort("LinOp::bCoefficients");
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

}
