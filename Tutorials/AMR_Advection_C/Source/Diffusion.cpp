#include <ParmParse.H>
#include "Diffusion.H"
#include "ADR.H"
#include "ADR_F.H"

#include <MacBndry.H>
#include <MGT_Solver.H>
#include <stencil_types.H>
#include <mg_cpp_f.h>

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
                    MultiFab& DiffTerm, PArray<MultiFab>& diff_coef)
{
    if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... compute diffusive term at level " << level << '\n';
    }

    int nlevs = 1;

    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

    if (level == 0) {
      for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        xa[0][i] = 0.;
        xb[0][i] = 0.;
      }
    } else {
      const Real* dx_crse   = parent->Geom(level-1).CellSize();
      for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        xa[0][i] = 0.5 * dx_crse[i];
        xb[0][i] = 0.5 * dx_crse[i];
      }
    }

    //
    // Store the Dirichlet boundary condition for phi in bndry.
    //
    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level],1,geom);

    // Build the homogeneous boundary conditions.  One could setVal
    // the bndry fabsets directly, but we instead do things as if
    // we had a fill-patched mf with grows--in that case the bndry
    // object knows how to grab grow data from the mf on physical
    // boundarys.  Here we creat an mf, setVal, and pass that to
    // the bndry object.
    //
    if (level == 0)
    {
        bndry.setBndryValues(Species,0,0,1,*phys_bc);
    }

    std::vector<BoxArray> bav(nlevs);
    std::vector<DistributionMapping> dmv(nlevs);

    ADR* cs = dynamic_cast<ADR*>(&parent->getLevel(level));
    bav[0] = grids[level];
    MultiFab& S_new = cs->get_new_data(State_Type);
    dmv[0] = S_new.DistributionMap();
    std::vector<Geometry> fgeom(1);
    fgeom[0] = parent->Geom(level);

    MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);
    
    MultiFab* phi_p[1];
    MultiFab* Res_p[1];

    Array< PArray<MultiFab> > coeffs(1);

    DiffTerm.setVal(0.);

    phi_p[0] = &Species;
    Res_p[0] = &DiffTerm;

    // Need to do this even if Cartesian because the array is needed in set_coefficients
    coeffs[0].resize(BL_SPACEDIM,PArrayManage);
    Geometry g = cs->Geom();
    for (int i = 0; i < BL_SPACEDIM ; i++) {
        coeffs[0].set(i, new MultiFab);
        g.GetFaceArea(coeffs[0][i],grids[level],i,0);
        MultiFab::Copy(coeffs[0][i],diff_coef[i],0,0,1,0);
    }

#if (BL_SPACEDIM < 3)
    // NOTE: we just pass Res here to use in the MFIter loop...
    if (Geometry::IsRZ() || Geometry::IsSPHERICAL())
       applyMetricTerms(level,(*Res_p[0]),coeffs[0]);
#endif

    mgt_solver.set_gravity_coefficients(coeffs,xa,xb,0);
 
    mgt_solver.applyop(phi_p, Res_p, bndry);

#if (BL_SPACEDIM < 3)
    // Do this to unweight Res
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
       unweight_cc(level,(*Res_p[0]));
#endif
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

    int nlevs = 1;

    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

    if (level == 0) {
      for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        xa[0][i] = 0.;
        xb[0][i] = 0.;
      }
    } else {
      const Real* dx_crse   = parent->Geom(level-1).CellSize();
      for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        xa[0][i] = 0.5 * dx_crse[i];
        xb[0][i] = 0.5 * dx_crse[i];
      }
    }

    //
    // Store the Dirichlet boundary condition for phi in bndry.
    //
    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level],1,geom);
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                   : IntVect::TheZeroVector();

    // Build the homogeneous boundary conditions.  One could setVal
    // the bndry fabsets directly, but we instead do things as if
    // we had a fill-patched mf with grows--in that case the bndry
    // object knows how to grab grow data from the mf on physical
    // boundarys.  Here we creat an mf, setVal, and pass that to
    // the bndry object.
    //
    BoxArray crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
    const int in_rad     = 0;
    const int out_rad    = 1; 
    const int extent_rad = 2;
    BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
    crse_br.copyFrom(Crse,1,src_comp,dest_comp,num_comp);
    bndry.setBndryValues(crse_br,src_comp,Species,src_comp,
                         dest_comp,num_comp,crse_ratio,*phys_bc);

    std::vector<BoxArray> bav(nlevs);
    std::vector<DistributionMapping> dmv(nlevs);

    ADR* cs = dynamic_cast<ADR*>(&parent->getLevel(level));
    bav[0] = grids[level];
    MultiFab& S_new = cs->get_new_data(State_Type);
    dmv[0] = S_new.DistributionMap();
    std::vector<Geometry> fgeom(1);
    fgeom[0] = parent->Geom(level);

    MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);
    
    MultiFab* phi_p[1];
    MultiFab* Res_p[1];

    Array< PArray<MultiFab> > coeffs(1);

    DiffTerm.setVal(0.);

    phi_p[0] = &Species;
    Res_p[0] = &DiffTerm;

    // Need to do this even if Cartesian because the array is needed in set_coefficients
    coeffs[0].resize(BL_SPACEDIM,PArrayManage);
    Geometry g = cs->Geom();
    for (int i = 0; i < BL_SPACEDIM ; i++) {
        coeffs[0].set(i, new MultiFab);
        g.GetFaceArea(coeffs[0][i],grids[level],i,0);
        MultiFab::Copy(coeffs[0][i],diff_coef[i],0,0,1,0);
    }

#if (BL_SPACEDIM < 3)
    // NOTE: we just pass Res here to use in the MFIter loop...
    if (Geometry::IsRZ() || Geometry::IsSPHERICAL())
       applyMetricTerms(level,(*Res_p[0]),coeffs[0]);
#endif

    mgt_solver.set_gravity_coefficients(coeffs,xa,xb,0);
 
    mgt_solver.applyop(phi_p, Res_p, bndry);

#if (BL_SPACEDIM < 3)
    // Do this to unweight Res
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
       unweight_cc(level,(*Res_p[0]));
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

