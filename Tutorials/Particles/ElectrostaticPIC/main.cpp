#include <iostream>
#include <random>
#include <cassert>

#include <AMReX.H>
#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_MacBndry.H>
#include "AMReX_FillPatchUtil.H"
#include "AMReX_PlotFileUtil.H"

#include "ElectrostaticParticleContainer.H"

#include "electrostatic_pic_F.H"

using namespace amrex;

void WritePlotFile(const ScalarMeshData& rhs,
                   const ScalarMeshData& phi,
                   const VectorMeshData& E,
                   const ElectrostaticParticleContainer& pc,
                   const Array<Geometry>& geom,
                   int nstep)
{
    int num_output_comp = 5;
    int num_levels = rhs.size();
    IntVect cc_flag = IntVect::TheZeroVector();
    Array<std::unique_ptr<MultiFab> > output_cc(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        const BoxArray& nodal_ba = rhs[lev]->boxArray();
        output_cc[lev].reset(new MultiFab(amrex::convert(nodal_ba, cc_flag), 
                                          rhs[lev]->DistributionMap(), num_output_comp, 0));
        amrex::average_node_to_cellcenter(*output_cc[lev], 0, *rhs[lev],  0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 1, *phi[lev],  0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 2, *E[lev][0], 0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 3, *E[lev][1], 0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 4, *E[lev][2], 0, 1);
    }
    
    Array<std::string> varnames;
    varnames.push_back("rhs");
    varnames.push_back("phi");
    varnames.push_back("Ex");
    varnames.push_back("Ey");
    varnames.push_back("Ez");

    Array<std::string> particle_varnames;
    particle_varnames.push_back("weight");
    particle_varnames.push_back("vx");
    particle_varnames.push_back("vy");
    particle_varnames.push_back("vz");
    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");
    particle_varnames.push_back("Ez");
    
    Array<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);
    
    int output_levs = num_levels;
    
    Array<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputRR[lev] = IntVect(2, 2, 2);
    }

    const std::string& pltfile = amrex::Concatenate("plt", nstep, 5);    
    WriteMultiLevelPlotfile(pltfile, output_levs, GetArrOfConstPtrs(output_cc),
                            varnames, geom, 0.0, level_steps, outputRR);

    pc.Checkpoint(pltfile, "particle0", true, particle_varnames);

}

void computeE(      VectorMeshData& E,
              const ScalarMeshData& phi, 
              const Array<Geometry>& geom) {

    const int num_levels = E.size();
    const int finest_level = num_levels - 1;

    for (int lev = 0; lev < num_levels; ++lev) {

        const auto& gm = geom[lev];
        const Real* dx = gm.CellSize();

        for (MFIter mfi(*phi[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            compute_E_nodal(bx.loVect(), bx.hiVect(),
                            (*phi[lev] )[mfi].dataPtr(),
                            (*E[lev][0])[mfi].dataPtr(),
                            (*E[lev][1])[mfi].dataPtr(),
                            (*E[lev][2])[mfi].dataPtr(), dx);
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
        E[lev][2]->FillBoundary(gm.periodicity());
    }
}

void computePhi(ScalarMeshData& rhs, ScalarMeshData& phi,
                Array<BoxArray>& grids,
                Array<DistributionMapping>& dm,
                Array<Geometry>& geom) {

    int num_levels = rhs.size();

    for (int lev = 0; lev < num_levels; ++lev) {
        phi[lev]->setVal(0.0);
    }
    
    PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};   
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    NodeBilinear mapper;

    bool nodal = true;
    bool have_rhcc = false;
    int  nc = 0;
    int Ncomp = 1;
    int stencil = ND_CROSS_STENCIL;
    int verbose = 0;
    Array<int> mg_bc(2*BL_SPACEDIM, 1); // this means Dirichlet
    Real rel_tol = 1.0e-9;
    Real abs_tol = 1.0e-9;

    MGT_Solver solver(geom, mg_bc.dataPtr(), grids, dm, nodal,
                      stencil, have_rhcc, nc, Ncomp, verbose);
    
    solver.set_nodal_const_coefficients(1.0);
    
    solver.solve_nodal(amrex::GetArrOfPtrs(phi), amrex::GetArrOfPtrs(rhs), rel_tol, abs_tol);

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }
        
    for (int lev = 0; lev < num_levels-1; ++lev) {
        // info for coarse/fine interpolation
        PhysBCFunct cphysbc, fphysbc;
        int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
        int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
        Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
        NodeBilinear mapper;

        MultiFab tmp(phi[lev+1]->boxArray(), phi[lev+1]->DistributionMap(), 1, 1);
                    Array<MultiFab*> crse(1);
            crse[0] = phi[lev].get();
            Array<Real> ctime(1, 0.0);
            Array<Real> ftime(1, 0.0);            
            Array<MultiFab*> fine(1);
            fine[0] = phi[lev+1].get();
            amrex::FillPatchTwoLevels(tmp, 0.0, crse, ctime, fine, ftime, 
                                      0, 0, 1, 
                                      geom[lev], geom[lev+1], cphysbc, fphysbc,
                                      IntVect(2, 2, 2), &mapper, bcs);
            MultiFab::Copy(*phi[lev+1], tmp, 0, 0, 1, 1);
    }

    // Array<Geometry>            level_geom(1);
    // Array<BoxArray>            level_grids(1);
    // Array<DistributionMapping> level_dm(1);
    // Array<MultiFab*>           level_phi(1);
    // Array<MultiFab*>           level_rhs(1);
    
    // for (int lev = 0; lev < num_levels; ++lev) {
    //     level_phi[0]   = phi[lev].get();
    //     level_rhs[0]   = rhs[lev].get();
    //     level_geom[0]  = geom[lev];
    //     level_grids[0] = grids[lev];
    //     level_dm[0]    = dm[lev];
        
    //     MGT_Solver solver(level_geom, mg_bc.dataPtr(), level_grids, 
    //                       level_dm, nodal,
    //                       stencil, have_rhcc, nc, Ncomp, verbose);
        
    //     solver.set_nodal_const_coefficients(1.0);
        
    //     solver.solve_nodal(level_phi, level_rhs, rel_tol, abs_tol);

    //     if (lev < num_levels-1) {
    //         amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0, *phi[lev],
    //                                      0, 0, 1, geom[lev], geom[lev+1],
    //                                      cphysbc, fphysbc,
    //                                      IntVect(2, 2, 2), &mapper, bcs);
    //     }
    // }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int max_level, n_cell, max_grid_size, max_step, is_periodic[BL_SPACEDIM];
    Real dt;

    // inputs parameters
    {
        ParmParse pp;
        pp.get("max_level", max_level);
        pp.get("n_cell", n_cell);
        pp.get("max_grid_size", max_grid_size);       
        pp.get("max_step", max_step);
        pp.get("dt", dt);
    }
    
    assert(max_level < 2);

    int num_levels = max_level + 1;
    
    Array<int> rr(num_levels-1);
    for (int lev = 1; lev < num_levels; lev++)
        rr[lev-1] = 2;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,-20.0e-6);
        real_box.setHi(n, 20.0e-6);
    }
    
    // This says we are using Cartesian coordinates
    int coord = 0;
    
    // This sets the boundary conditions to be doubly or triply periodic
    for (int i = 0; i < BL_SPACEDIM; i++) {
        is_periodic[i] = 1;
    }

    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    Box domain(dom_lo, dom_hi);
   
    // make Geometry for each level
    Array<Geometry> geom(num_levels);
    geom[0].define(domain,&real_box,coord,is_periodic);
    for (int lev = 1; lev < num_levels; ++lev) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, coord, is_periodic);
    }
    
    // make grids for each level
    Array<BoxArray> grids(num_levels);
    grids[0].define(domain);
    if (num_levels > 1) {
        int n_fine = n_cell*rr[0];
        IntVect refined_lo(n_fine/4,n_fine/4,n_fine/4); 
        IntVect refined_hi(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[1].define(refined_patch);
    }
    
    for (int lev = 0; lev < num_levels; lev++) {
        grids[lev].maxSize(max_grid_size);
    }
    
    int Nghost = 1;
    int Ncomp  = 1;
    
    Array<DistributionMapping> dm(num_levels);
    Array<std::unique_ptr<MultiFab> > phi(num_levels);
    Array<std::unique_ptr<MultiFab> > rhs(num_levels);
    Array<std::array<std::unique_ptr<MultiFab>, 3> > eField(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();
        dm[lev].define(grids[lev]);

        rhs[lev].reset(new MultiFab(nba, dm[lev], Ncomp, Nghost));
        phi[lev].reset(new MultiFab(nba, dm[lev], Ncomp, Nghost));
        
        eField[lev][0].reset(new MultiFab(nba, dm[lev], Ncomp, Nghost));
        eField[lev][1].reset(new MultiFab(nba, dm[lev], Ncomp, Nghost));
        eField[lev][2].reset(new MultiFab(nba, dm[lev], Ncomp, Nghost));

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);
        
        eField[lev][0]->setVal(0.0);
        eField[lev][1]->setVal(0.0);
        eField[lev][2]->setVal(0.0);
    }
    
    ElectrostaticParticleContainer myPC(geom, dm, grids, rr);

    myPC.InitParticles();

    for (int step = 0; step <= max_step; ++step) {

        myPC.DepositCharge(rhs);

        computePhi(rhs, phi, grids, dm, geom);
        
        computeE(eField, phi, geom);

        //        WritePlotFile(rhs, phi, eField, myPC, geom, step);

        myPC.FieldGather(eField);

        myPC.writeParticles(step);
        
        myPC.Evolve(eField, rhs, dt);

        myPC.Redistribute();        
    }
    
    amrex::Finalize();
}
