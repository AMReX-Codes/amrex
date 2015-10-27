#include <ParmParse.H>
#include "Diffusion.H"
#include "ADR.H"
#include "ADR_F.H"

#include <FMultiGrid.H>

#define MAX_LEV 15

int  Diffusion::verbose      = 0;
int  Diffusion::stencil_type = CC_CROSS_STENCIL;
 
Diffusion::Diffusion(Amr* Parent, BCRec* _phys_bc)
  : 
    parent(Parent),
    LevelData(MAX_LEV),
    phi_flux_reg(MAX_LEV, PArrayManage),
    grids(MAX_LEV),
    volume(MAX_LEV),
    area(MAX_LEV),
    phys_bc(_phys_bc)
{
    read_params();
    make_mg_bc();
}

Diffusion::~Diffusion() {}

void
Diffusion::read_params ()
{
    static bool done = false;

    if (!done)
    {
        ParmParse pp("diffusion");
        pp.query("v", verbose);
        done = true;
    }
}

void
Diffusion::install_level (int                   level,
                          AmrLevel*             level_data,
                          MultiFab&             _volume,
                          MultiFab*             _area)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing Diffusion level " << level << '\n';

    LevelData.clear(level);
    LevelData.set(level, level_data);

    volume.clear(level);
    volume.set(level, &_volume);

    area.set(level, _area);

    BoxArray ba(LevelData[level].boxArray());
    grids[level] = ba;

    if (level > 0) {
       phi_flux_reg.clear(level);
       IntVect crse_ratio = parent->refRatio(level-1);
       phi_flux_reg.set(level,new FluxRegister(grids[level],crse_ratio,level,1));
    }
}

void
Diffusion::zeroPhiFluxReg (int level)
{
  phi_flux_reg[level].setVal(0.);
}

void
Diffusion::applyop (int level, MultiFab& Species, 
                    MultiFab& Crse, MultiFab& DiffTerm, 
                    PArray<MultiFab>& diff_coef)
{
    if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... compute diffusive term at level " << level << '\n';
    }

    PArray<MultiFab> coeffs_curv;
#if (BL_SPACEDIM < 3)
    // NOTE: we just pass DiffTerm here to use in the MFIter loop...
    if (Geometry::IsRZ() || Geometry::IsSPHERICAL())
    {
	coeffs_curv.resize(BL_SPACEDIM, PArrayManage);

	for (int i = 0; i< BL_SPACEDIM; ++i) {
	    coeffs_curv.set(i, new MultiFab(diff_coef[i].boxArray(), 1, 0, Fab_allocate));
	    MultiFab::Copy(coeffs_curv[i], diff_coef[i], 0, 0, 1, 0);
	}

	applyMetricTerms(0, DiffTerm, coeffs_curv);
    }
#endif

    PArray<MultiFab> & coeffs = (coeffs_curv.size() > 0) ? coeffs_curv : diff_coef;

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                   : IntVect::TheZeroVector();

    FMultiGrid fmg(parent->Geom(level), level, crse_ratio);

    if (level == 0) {
	fmg.set_bc(mg_bc, Species);
    } else {
	fmg.set_bc(mg_bc, Crse, Species);
    }

    fmg.set_diffusion_coeffs(coeffs);

    fmg.applyop(Species, DiffTerm);

#if (BL_SPACEDIM < 3)
    // Do this to unweight Res
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
	unweight_cc(level, DiffTerm);
#endif
}

#if (BL_SPACEDIM < 3)
void
Diffusion::applyMetricTerms(int level, MultiFab& Rhs, PArray<MultiFab>& coeffs)
{
    const Real* dx = parent->Geom(level).CellSize();
    int coord_type = Geometry::Coord();
    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        // Modify Rhs and coeffs with the appropriate metric terms.
        BL_FORT_PROC_CALL(APPLY_METRIC,apply_metric)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(Rhs[mfi]),
             D_DECL(BL_TO_FORTRAN(coeffs[0][mfi]),
                    BL_TO_FORTRAN(coeffs[1][mfi]),
                    BL_TO_FORTRAN(coeffs[2][mfi])),
                    dx,&coord_type);
    }
}
#endif

#if (BL_SPACEDIM < 3)
void
Diffusion::unweight_cc(int level, MultiFab& cc)
{
    const Real* dx = parent->Geom(level).CellSize();
    int coord_type = Geometry::Coord();
    for (MFIter mfi(cc); mfi.isValid(); ++mfi)
    {
        int index = mfi.index();
        const Box& bx = grids[level][index];
        BL_FORT_PROC_CALL(UNWEIGHT_CC,unweight_cc)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(cc[mfi]),dx,&coord_type);
    }
}
#endif

void
Diffusion::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);
    for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
    {
        if ( geom.isPeriodic(dir) )
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
        else
        {
            if (phys_bc->lo(dir) == Symmetry   || 
                phys_bc->lo(dir) == SlipWall   || 
                phys_bc->lo(dir) == NoSlipWall || 
                phys_bc->lo(dir) == Outflow)   {
              mg_bc[2*dir + 0] = MGT_BC_NEU;
            }
            if (phys_bc->hi(dir) == Symmetry   || 
                phys_bc->hi(dir) == SlipWall   || 
                phys_bc->hi(dir) == NoSlipWall || 
                phys_bc->hi(dir) == Outflow)   {
              mg_bc[2*dir + 1] = MGT_BC_NEU;
            }
        }
    }

    // Set Neumann bc at r=0.
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
        mg_bc[0] = MGT_BC_NEU;
}

