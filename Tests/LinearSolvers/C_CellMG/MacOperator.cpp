

#include <AMReX_MacBndry.H>
#include <MacOperator.H>
#include <MacOpMacDrivers.H>
#include <MACOPERATOR_F.H>
#include <AMReX_CGSolver.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();
//
// This is the default maxorder for the linear operator.
//
// To change this ParmParse in a new value for "macop.max_order"
//
static int max_order;

namespace
{
    bool initialized = false;
}

void
MacOperator::Initialize ()
{
    if (initialized) return;
    //
    // This is the default maxorder for the linear operator.
    //
    // To change this ParmParse in a new value for "macop.max_order"
    //
    max_order = 4;

    ParmParse pp("macop");

    pp.query("max_order", max_order);

    amrex::ExecOnFinalize(MacOperator::Finalize);

    initialized = true;
}

void
MacOperator::Finalize ()
{
    initialized = false;
}

MacOperator::MacOperator (AmrCore*             Parent,
                          const BndryData& mgb,
                          const Real*      h)
    :
    ABecLaplacian(mgb,h),
    parent(Parent)
{
    Initialize();
}

MacOperator::~MacOperator () {}

//
// Define the meaning of gradient for the multigrid object.
//

void
MacOperator::setCoefficients (MultiFab*   area,
                              MultiFab&   rho,
                              int         rho_comp,
                              const Real* dx)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = boxArray(0);
    const DistributionMapping& dm = DistributionMap();
    BL_ASSERT(rho.boxArray() == ba);
    BL_ASSERT(rho.DistributionMap() == dm);
    //
    // First set scalar coeficients.
    //
    setScalars(0.0,1.0);
    //
    // Don't need to set a because alpha is set to zero.
    //
    const int n_grow = 0;

    AMREX_D_TERM(MultiFab bxcoef(area[0].boxArray(),dm,area[0].nComp(),n_grow);,
           MultiFab bycoef(area[1].boxArray(),dm,area[1].nComp(),n_grow);,
           MultiFab bzcoef(area[2].boxArray(),dm,area[2].nComp(),n_grow););
    AMREX_D_TERM(bxcoef.setVal(0);,
           bycoef.setVal(0);,
           bzcoef.setVal(0););

    for (MFIter rhomfi(rho); rhomfi.isValid(); ++rhomfi)
    {
        BL_ASSERT(ba[rhomfi.index()] == rhomfi.validbox());

        const Box& grd       = ba[rhomfi.index()];
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& bx        = bxcoef[rhomfi];
        FArrayBox& by        = bycoef[rhomfi];
        const FArrayBox& ax  = area[0][rhomfi];
        const FArrayBox& ay  = area[1][rhomfi];
        const FArrayBox& den = rho[rhomfi];

        DEF_LIMITS(bx,bx_dat,bxlo,bxhi);
        DEF_LIMITS(by,by_dat,bylo,byhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);

        const int* dlo      = den.loVect();
        const int* dhi      = den.hiVect();
        const Real* den_dat = den.dataPtr(rho_comp);

#if (BL_SPACEDIM == 2)
        FORT_MACCOEF(bx_dat,AMREX_ARLIM(bxlo),AMREX_ARLIM(bxhi),
                     by_dat,AMREX_ARLIM(bylo),AMREX_ARLIM(byhi),
                     ax_dat,AMREX_ARLIM(axlo),AMREX_ARLIM(axhi),
                     ay_dat,AMREX_ARLIM(aylo),AMREX_ARLIM(ayhi),
                     den_dat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),lo,hi,dx);
#endif
#if (BL_SPACEDIM == 3)
        FArrayBox& bz       = bzcoef[rhomfi];
        const FArrayBox& az = area[2][rhomfi];

        DEF_CLIMITS(az,az_dat,azlo,azhi);
        DEF_LIMITS(bz,bz_dat,bzlo,bzhi);

        FORT_MACCOEF(bx_dat,AMREX_ARLIM(bxlo),AMREX_ARLIM(bxhi),
                     by_dat,AMREX_ARLIM(bylo),AMREX_ARLIM(byhi),
                     bz_dat,AMREX_ARLIM(bzlo),AMREX_ARLIM(bzhi),
                     ax_dat,AMREX_ARLIM(axlo),AMREX_ARLIM(axhi),
                     ay_dat,AMREX_ARLIM(aylo),AMREX_ARLIM(ayhi),
                     az_dat,AMREX_ARLIM(azlo),AMREX_ARLIM(azhi),
                     den_dat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),lo,hi,dx);
#endif
    }
  
    AMREX_D_TERM(bCoefficients(bxcoef,0);,
           bCoefficients(bycoef,1);,
           bCoefficients(bzcoef,2););
}

//
// This function creates the initial rhs for use in the mac multgrid solve.
//

void
MacOperator::defRHS (MultiFab* area,
                     MultiFab& volume,
                     MultiFab& Rhs,
                     MultiFab* vel,
                     Real      scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    BL_ASSERT(Rhs.boxArray() == gbox[0]);

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        const FArrayBox& ax  = area[0][Rhsmfi];
        const FArrayBox& ay  = area[1][Rhsmfi];
        const FArrayBox& vol = volume[Rhsmfi];
        const FArrayBox& ux  = vel[0][Rhsmfi];
        const FArrayBox& uy  = vel[1][Rhsmfi];
        FArrayBox& rhs       = Rhs[Rhsmfi];

        DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);

#if (BL_SPACEDIM == 2)
        FORT_MACRHS(ux_dat,AMREX_ARLIM(uxlo),AMREX_ARLIM(uxhi),
                    uy_dat,AMREX_ARLIM(uylo),AMREX_ARLIM(uyhi),
                    ax_dat,AMREX_ARLIM(axlo),AMREX_ARLIM(axhi),
                    ay_dat,AMREX_ARLIM(aylo),AMREX_ARLIM(ayhi),
                    vol_dat,AMREX_ARLIM(vlo),AMREX_ARLIM(vhi), 
                    rhs_dat,AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),
                    lo,hi,&scale);
#endif
#if (BL_SPACEDIM == 3)
        const FArrayBox& az = area[2][Rhsmfi];
        DEF_CLIMITS(az,az_dat,azlo,azhi);

        const FArrayBox& uz = vel[2][Rhsmfi];
        DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);

        FORT_MACRHS(ux_dat,AMREX_ARLIM(uxlo),AMREX_ARLIM(uxhi),
                    uy_dat,AMREX_ARLIM(uylo),AMREX_ARLIM(uyhi),
                    uz_dat,AMREX_ARLIM(uzlo),AMREX_ARLIM(uzhi),
                    ax_dat,AMREX_ARLIM(axlo),AMREX_ARLIM(axhi),
                    ay_dat,AMREX_ARLIM(aylo),AMREX_ARLIM(ayhi),
                    az_dat,AMREX_ARLIM(azlo),AMREX_ARLIM(azhi),
                    vol_dat,AMREX_ARLIM(vlo),AMREX_ARLIM(vhi),
                    rhs_dat,AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),
                    lo,hi,&scale);
#endif
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Apply the mac pressure gradient to a velocity field.
// init, means that velocities are initialized here.
//

void
mac_vel_update (int              init,
                AMREX_D_DECL(FArrayBox& ux,
                       FArrayBox& uy,
                       FArrayBox& uz),
                const FArrayBox& phi,
                const FArrayBox* rhoptr,
                int              rho_comp,  
                const Box&       grd,
                int              level,
                int              n,
                const Real*      dx,
                Real             scale)
{
    const int* lo        = grd.loVect();
    const int* hi        = grd.hiVect();

    const FArrayBox& rho = *rhoptr;
    
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);

    const int* rlo      = rho.loVect();
    const int* rhi      = rho.hiVect();
    const Real* rho_dat = rho.dataPtr(rho_comp);
    
#if (BL_SPACEDIM == 2)
    FORT_MACUPDATE(&init,
                   ux_dat,AMREX_ARLIM(uxlo),AMREX_ARLIM(uxhi),
                   uy_dat,AMREX_ARLIM(uylo),AMREX_ARLIM(uyhi),
                   phi_dat,AMREX_ARLIM(p_lo),AMREX_ARLIM(p_hi),
                   rho_dat,AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),
                   lo,hi,dx,&scale);
#endif
#if (BL_SPACEDIM == 3)
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    
    FORT_MACUPDATE(&init,
                   ux_dat,AMREX_ARLIM(uxlo),AMREX_ARLIM(uxhi),
                   uy_dat,AMREX_ARLIM(uylo),AMREX_ARLIM(uyhi),
                   uz_dat,AMREX_ARLIM(uzlo),AMREX_ARLIM(uzhi),
                   phi_dat,AMREX_ARLIM(p_lo),AMREX_ARLIM(p_hi),
                   rho_dat,AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),
                   lo,hi,dx,&scale);
#endif
}

//
// Apply the mac pressure gradient to the divergent mac velocities.
// The resultant velocity field is nondivergent.
//

void
MacOperator::velUpdate (MultiFab*       Vel,
                        MultiFab&       Phi,
                        const MultiFab& Rho,
                        int             rho_comp,
                        const Real*     dx,
                        Real            scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    BL_ASSERT(Rho.boxArray() ==  gbox[0]);
    //
    // Set bndry data in ghost zones.
    //
    int apply_lev = 0;
    applyBC(Phi,0,1,apply_lev);

    for (MFIter Phimfi(Phi); Phimfi.isValid(); ++Phimfi)
    {
        const Box& grd = Phimfi.validbox();

        mac_vel_update(0, 
                       AMREX_D_DECL(Vel[0][Phimfi],Vel[1][Phimfi],Vel[2][Phimfi]),
                       Phi[Phimfi],
                       &(Rho[Phimfi]), rho_comp,  
                       grd, 0, Phimfi.index(), dx, scale );
    }
}

//
// Multiply by volume*rhs_scale since reflux step (which computed rhs)
// divided by volume.
//

void
MacOperator::syncRhs (const MultiFab& Volume,
                      MultiFab&       Rhs,
                      Real            rhs_scale,
                      const Real*     dx)
{
    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& rhs       = Rhs[Rhsmfi];
        const FArrayBox& vol = Volume[Rhsmfi];

        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);
        FORT_MACSYNCRHS(rhs_dat,AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),lo,hi,
                        vol_dat,AMREX_ARLIM(vlo),AMREX_ARLIM(vhi),&rhs_scale);
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Driver functions follow.
//

//
// A driver function for computing a level MAC solve.
//

void
mac_level_driver (AmrCore*        parent,
                  const MacBndry& mac_bndry,
		  const BCRec&    phys_bc,
                  const BoxArray& grids,
                  int             the_solver,
                  int             level,
                  int             Density,
                  const Real*     dx,
                  Real            dt,
                  Real            mac_tol,
                  Real            mac_abs_tol,
                  Real            rhs_scale,
                  MultiFab*       area,
                  MultiFab&       volume,
                  MultiFab&       S,
                  MultiFab&       Rhs,
                  MultiFab*       u_mac,
                  MultiFab*       mac_phi)
{
    MacOperator mac_op(parent,mac_bndry,dx);

    mac_op.setCoefficients(area,S,Density,dx);
    mac_op.defRHS(area,volume,Rhs,u_mac,rhs_scale);
    mac_op.maxOrder(max_order);

    if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
        amrex::Error("Can't use CGSolver with maxorder > 2");
    }
    //
    // Construct MultiGrid or CGSolver object and solve system.
    //
    if (the_solver == 1)
    {
        bool use_mg_precond = true;
        CGSolver mac_cg(mac_op,use_mg_precond);
        mac_cg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
    else if (the_solver == 2 )
    {
        amrex::Error("mac_level_driver::HypreABec not in this build");
    }
    else
    {
        MultiGrid mac_mg(mac_op);
        mac_mg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
    //
    // velUpdate will set bndry values for mac_phi.
    //
    mac_op.velUpdate(u_mac,*mac_phi,S,Density,dx,-dt/2.0);
}

//
// A driver function for computing a sync MAC solve.
//

void
mac_sync_driver (AmrCore*            parent,
                 const MacBndry& mac_bndry,
	         const BCRec&    phys_bc,
                 const BoxArray& grids,
                 int             the_solver,
                 int             level, 
                 const Real*     dx,
                 Real            dt,
                 Real            mac_sync_tol,
                 Real            mac_abs_tol,
                 Real            rhs_scale,
                 MultiFab*       area,
                 MultiFab&       volume,
                 MultiFab&       Rhs,
                 MultiFab*       rho_half,
                 MultiFab*       mac_sync_phi)
{
    MacOperator mac_op(parent,mac_bndry,dx);

    mac_op.maxOrder(max_order);
    mac_op.setCoefficients(area,*rho_half, 0, dx);
    mac_op.syncRhs(volume,Rhs,rhs_scale,dx);

    if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
        amrex::Error("Can't use CGSolver with maxorder > 2");
    }
    //
    // Now construct MultiGrid or CGSolver object to solve system.
    //
    if (the_solver == 1)
    {
        bool use_mg_precond = true;
        CGSolver mac_cg(mac_op,use_mg_precond);
        mac_cg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
    else if ( the_solver == 2 )
    {
        amrex::Error("mac_sync_driver: HypreABec not in this build");
    }
    else
    {
        MultiGrid mac_mg(mac_op);
        mac_mg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
    
    int mac_op_lev = 0;
    mac_op.applyBC(*mac_sync_phi,0,1,mac_op_lev);
}
