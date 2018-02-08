#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include "Constants.H"
#include "em_pic_F.H"

using namespace amrex;

IntVect Bx_nodal_flag(1,0,0);
IntVect By_nodal_flag(0,1,0);
IntVect Bz_nodal_flag(0,0,1);

IntVect Ex_nodal_flag(0,1,1);
IntVect Ey_nodal_flag(1,0,1);
IntVect Ez_nodal_flag(1,1,0);

IntVect jx_nodal_flag(0,1,1);
IntVect jy_nodal_flag(1,0,1);
IntVect jz_nodal_flag(1,1,0);

struct TestParams {
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
    bool verbose;
};

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
{
    int nx = nppc[0];
    int ny = nppc[1];
    int nz = nppc[2];

    int ix_part = i_part/(ny * nz);
    int iy_part = (i_part % (ny * nz)) % ny;
    int iz_part = (i_part % (ny * nz)) / ny;
    
    r[0] = (0.5+ix_part)/nx;
    r[1] = (0.5+iy_part)/ny;
    r[2] = (0.5+iz_part)/nz;    
}

void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
    Real ux_th = amrex::RandomNormal(0.0, u_std);
    Real uy_th = amrex::RandomNormal(0.0, u_std);
    Real uz_th = amrex::RandomNormal(0.0, u_std);
    
    u[0] = u_mean + ux_th;
    u[1] = u_mean + uy_th;
    u[2] = u_mean + uz_th;
}

void set_initial_conditions(Array<Real>& xp, Array<Real>& yp, Array<Real>& zp,
                            Array<Real>& uxp, Array<Real>& uyp, Array<Real>& uzp,
                            Array<Real>& w, const Geometry& geom, const BoxArray& ba,
                            const DistributionMapping& dm, const TestParams& parms)
{
    
    BL_PROFILE("set_initial_conditions");
    
    const int num_ppc = AMREX_D_TERM(parms.nppc[0], *parms.nppc[1], *parms.nppc[2]);
    const Real* dx = geom.CellSize();
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
    
    const Real thermal_momentum_std  = 0.01;
    const Real thermal_momentum_mean = 10;
    const Real density               = 1e25;
    
    MultiFab dummy_mf(ba, dm, 1, 0, MFInfo().SetAlloc(false));
    for(MFIter mfi(dummy_mf); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const Real* tile_lo = tile_realbox.lo();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) 
            {
                Real r[3];
                Real u[3];
                
                get_position_unit_cell(r, parms.nppc, i_part);
                get_gaussian_random_momentum(u, thermal_momentum_mean, thermal_momentum_std);
                
                Real x = tile_lo[0] + (iv[0] + r[0])*dx[0];
                Real y = tile_lo[1] + (iv[1] + r[1])*dx[1];
                Real z = tile_lo[2] + (iv[2] + r[2])*dx[2];
                
                xp.push_back(x);
                yp.push_back(y);
                zp.push_back(z);

                uxp.push_back(u[0]);
                uyp.push_back(u[1]);
                uzp.push_back(u[2]);

                w.push_back(density * scale_fac);                
            }
        }
    }
}

Real compute_dt(const Geometry& geom) 
{
    const static Real cfl = 1.0;
    const Real* dx = geom.CellSize();
    const Real deltat  = cfl * 1./( std::sqrt(D_TERM(  1./(dx[0]*dx[0]),
                                                       + 1./(dx[1]*dx[1]),
                                                       + 1./(dx[2]*dx[2]))) * PhysConst::c );
}

void evolve_electric_field(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez,
                           const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                           const MultiFab& jx, const MultiFab& jy, const MultiFab& jz,
                           const Geometry& geom, Real dt)
{
    BL_PROFILE("evolve_electric_field");

    static constexpr Real c2 = PhysConst::c*PhysConst::c;
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real foo = (PhysConst::c*PhysConst::c) * dt;

    const Real* dx = geom.CellSize();    
    const std::array<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi) 
    {    
        const Box& tbx  = mfi.tilebox(Ex_nodal_flag);
        const Box& tby  = mfi.tilebox(Ey_nodal_flag);
        const Box& tbz  = mfi.tilebox(Ez_nodal_flag);
       
        push_electric_field(tbx.loVect(), tbx.hiVect(), 
                            tby.loVect(), tby.hiVect(),
                            tbz.loVect(), tbz.hiVect(),
                            BL_TO_FORTRAN_3D(Ex[mfi]),
                            BL_TO_FORTRAN_3D(Ey[mfi]),
                            BL_TO_FORTRAN_3D(Ez[mfi]),
                            BL_TO_FORTRAN_3D(Bx[mfi]),
                            BL_TO_FORTRAN_3D(By[mfi]),
                            BL_TO_FORTRAN_3D(Bz[mfi]),
                            BL_TO_FORTRAN_3D(jx[mfi]),
                            BL_TO_FORTRAN_3D(jy[mfi]),
                            BL_TO_FORTRAN_3D(jz[mfi]),
                            &mu_c2_dt, &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2]);
    }
}

void evolve_magnetic_field(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                           MultiFab& Bx, MultiFab& By, MultiFab& Bz,
                           const Geometry& geom, Real dt)
{
    BL_PROFILE("evolve_magnetic_field");

    const Real* dx = geom.CellSize();
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi) 
    {    
        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);
       
        push_magnetic_field(tbx.loVect(), tbx.hiVect(), 
                            tby.loVect(), tby.hiVect(),
                            tbz.loVect(), tbz.hiVect(),
                            BL_TO_FORTRAN_3D(Ex[mfi]),
                            BL_TO_FORTRAN_3D(Ey[mfi]),
                            BL_TO_FORTRAN_3D(Ez[mfi]),
                            BL_TO_FORTRAN_3D(Bx[mfi]),
                            BL_TO_FORTRAN_3D(By[mfi]),
                            BL_TO_FORTRAN_3D(Bz[mfi]),
                            &dtsdx[0], &dtsdx[1], &dtsdx[2]);
    }
}

void push_and_depose_particles(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                               const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                               MultiFab& jx, MultiFab& jy, MultiFab& jz,
                               Array<Real>& xp,  Array<Real>& yp,  Array<Real>& zp,
                               Array<Real>& uxp, Array<Real>& uyp, Array<Real>& uzp,
                               Array<Real>& exp, Array<Real>& eyp, Array<Real>& ezp,
                               Array<Real>& bxp, Array<Real>& byp, Array<Real>& bzp,
                               Array<Real>& ginvp, Array<Real>& w, const Geometry& geom,
                               Real dt)
{
    BL_PROFILE_VAR_NS("field_gather", blp_fg);
    BL_PROFILE_VAR_NS("particle_push", blp_pp);
    BL_PROFILE_VAR_NS("current_deposition", blp_cd);

    int np = xp.size();
    const Real* dx  = geom.CellSize();
    const Real* plo = geom.ProbLo();

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi) 
    {

        BL_PROFILE_VAR_START(blp_fg);

        gather_magnetic_field(&np, xp.data(), yp.data(), zp.data(),
                              bxp.data(), byp.data(), bzp.data(),
                              BL_TO_FORTRAN_3D(Bx[mfi]),
                              BL_TO_FORTRAN_3D(By[mfi]),
                              BL_TO_FORTRAN_3D(Bz[mfi]),
                              plo, dx);
        
        gather_electric_field(&np, xp.data(), yp.data(), zp.data(),
                              exp.data(), eyp.data(), ezp.data(),
                              BL_TO_FORTRAN_3D(Ex[mfi]),
                              BL_TO_FORTRAN_3D(Ey[mfi]),
                              BL_TO_FORTRAN_3D(Ez[mfi]),
                              plo, dx);

        BL_PROFILE_VAR_STOP(blp_fg);

        BL_PROFILE_VAR_START(blp_pp);

        push_momentum_boris(&np, uxp.data(), uyp.data(), uzp.data(), ginvp.data(),
                            exp.data(), eyp.data(), ezp.data(),
                            bxp.data(), byp.data(), bzp.data(),
                            &PhysConst::q_e, &PhysConst::m_e, &dt);

        push_position_boris(&np, xp.data(), yp.data(), zp.data(),
                            uxp.data(), uyp.data(), uzp.data(), ginvp.data(), &dt);

        BL_PROFILE_VAR_STOP(blp_pp);

        BL_PROFILE_VAR_START(blp_cd);

        deposit_current(BL_TO_FORTRAN_3D(jx[mfi]),
                        BL_TO_FORTRAN_3D(jy[mfi]),
                        BL_TO_FORTRAN_3D(jz[mfi]),
                        &np, xp.data(), yp.data(), zp.data(),
                        uzp.data(), uyp.data(), uzp.data(), 
                        ginvp.data(), w.data(), &PhysConst::q_e,
                        plo, &dt, dx);

        BL_PROFILE_VAR_STOP(blp_cd);

    }
}

void push_particles_only(Array<Real>& xp,  Array<Real>& yp,  Array<Real>& zp,
                         Array<Real>& uxp, Array<Real>& uyp, Array<Real>& uzp,
                         Array<Real>& ginvp, Real dt)
{
    BL_PROFILE("particle_push");
    
    int np = xp.size();
    
    set_gamma(&np, uxp.data(), uyp.data(), uzp.data(), ginvp.data());
    
    push_position_boris(&np, xp.data(), yp.data(), zp.data(),
                        uxp.data(), uyp.data(), uzp.data(), ginvp.data(), &dt);
}

void enforce_periodic_bcs(Array<Real>& xp,  Array<Real>& yp,  Array<Real>& zp, 
                          const Geometry& geom)
{
    BL_PROFILE("enforce_periodic");
    
    int np = xp.size();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    enforce_periodic(&np, xp.data(), yp.data(), zp.data(), plo, phi);
}
                                                              
void test_em_pic(const TestParams& parms)
{

    BL_PROFILE("test_em_pic");

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0)); 
    IntVect domain_hi(AMREX_D_DECL(parms.ncell[0]-1,parms.ncell[1]-1,parms.ncell[2]-1)); 
    const Box domain(domain_lo, domain_hi);
    
    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 1; 
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    //  note - no domain decomposition right now
    //  ba.maxSize(parms.max_grid_size);
    if (parms.verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Number of boxes              : " << ba.size() << '\n' << '\n';
    }    
    DistributionMapping dm(ba);
    
    const int ng = 1;
    
    MultiFab Bx(amrex::convert(ba, Bx_nodal_flag), dm, 1, ng);
    MultiFab By(amrex::convert(ba, By_nodal_flag), dm, 1, ng);
    MultiFab Bz(amrex::convert(ba, Bz_nodal_flag), dm, 1, ng);
    
    MultiFab Ex(amrex::convert(ba, Bx_nodal_flag), dm, 1, ng);
    MultiFab Ey(amrex::convert(ba, By_nodal_flag), dm, 1, ng);
    MultiFab Ez(amrex::convert(ba, Bz_nodal_flag), dm, 1, ng);
    
    MultiFab jx(amrex::convert(ba, Bx_nodal_flag), dm, 1, ng);
    MultiFab jy(amrex::convert(ba, By_nodal_flag), dm, 1, ng);
    MultiFab jz(amrex::convert(ba, Bz_nodal_flag), dm, 1, ng);
    
    Bx.setVal(0.0); By.setVal(0.0); Bz.setVal(0.0);
    Ex.setVal(0.0); Ey.setVal(0.0); Ez.setVal(0.0);
    jx.setVal(0.0); jy.setVal(0.0); jz.setVal(0.0);
    
    Array<Real>  xp; Array<Real>  yp; Array<Real>  zp;
    Array<Real> uxp; Array<Real> uyp; Array<Real> uzp;
    Array<Real> exp; Array<Real> eyp; Array<Real> ezp;
    Array<Real> bxp; Array<Real> byp; Array<Real> bzp;
    
    Array<Real> ginvp;
    Array<Real> w;

    set_initial_conditions(xp, yp, zp, uxp, uyp, uzp, w, geom, ba, dm, parms);
    
    int np = xp.size();
    int nsteps = parms.nsteps;
    const Real dt = compute_dt(geom);
    bool synchronized = true;
    
    exp.resize(np, 0); eyp.resize(np, 0); ezp.resize(np, 0);
    bxp.resize(np, 0); byp.resize(np, 0); bzp.resize(np, 0);
    ginvp.resize(np, 0);
    
    for (int step = 0; step < nsteps; ++step) { 
        
        if (synchronized) {
            evolve_electric_field(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, 0.5*dt);
            push_particles_only(xp, yp, zp, uzp, uyp, uzp, ginvp, 0.5*dt);
            enforce_periodic_bcs(xp, yp, zp, geom);
            synchronized = false;
        }
        
        evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
        
        push_and_depose_particles(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz,
                                  xp, yp, zp, uxp, uyp, uzp, exp, eyp, ezp, bxp, byp, bzp,
                                  ginvp, w, geom, dt);
        
        enforce_periodic_bcs(xp, yp, zp, geom);
        
        evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
        
        if (step == nsteps - 1) {
            evolve_electric_field(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, dt);
            push_particles_only(xp, yp, zp, uzp, uyp, uzp, ginvp, -0.5*dt);
            synchronized = true;
        } else {
            evolve_electric_field(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, dt);
        }
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    ParmParse pp;
    
    TestParams parms;
    
    pp.get("ncell", parms.ncell);
    pp.get("nppc",  parms.nppc);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nsteps", parms.nsteps);

    parms.verbose = false;
    pp.query("verbose", parms.verbose);
    
    if (parms.verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl;
        std::cout << "Number of particles per cell : ";
        std::cout << parms.nppc  << std::endl;
        std::cout << "Size of domain               : ";
        std::cout << parms.ncell[0] << " " 
                  << parms.ncell[1] << " " 
                  << parms.ncell[2] << std::endl;
    }
    
    test_em_pic(parms);
    
    amrex::Finalize();
}
