#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>

#include <AMReX.H>

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_MacBndry.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "ElectrostaticParticleContainer.H"

#include "electrostatic_pic_F.H"

using namespace amrex;

class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time, int) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

void WritePlotFile(const ScalarMeshData& rhs,
                   const ScalarMeshData& phi,
                   const VectorMeshData& E,
                   const ElectrostaticParticleContainer& pc,
                   const Vector<Geometry>& geom,
                   int nstep)
{
    int num_output_comp = 2 + BL_SPACEDIM;
    int num_levels = rhs.size();
    IntVect cc_flag = IntVect::TheZeroVector();
    Vector<std::unique_ptr<MultiFab> > output_cc(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        const BoxArray& nodal_ba = rhs[lev]->boxArray();
        output_cc[lev].reset(new MultiFab(amrex::convert(nodal_ba, cc_flag), 
                                          rhs[lev]->DistributionMap(), num_output_comp, 0));
        amrex::average_node_to_cellcenter(*output_cc[lev], 0, *rhs[lev],  0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 1, *phi[lev],  0, 1);
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            amrex::average_node_to_cellcenter(*output_cc[lev], 2+i, *E[lev][i], 0, 1);
        }
    }
    
    Vector<std::string> varnames;
    varnames.push_back("rhs");
    varnames.push_back("phi");
    varnames.push_back("Ex");
    varnames.push_back("Ey");
#if BL_SPACEDIM == 3
    varnames.push_back("Ez");
#endif

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("weight");
    particle_varnames.push_back("vx");
    particle_varnames.push_back("vy");
#if BL_SPACEDIM == 3
    particle_varnames.push_back("vz");
#endif
    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");
#if BL_SPACEDIM == 3
    particle_varnames.push_back("Ez");
#endif
    
    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);
    
    int output_levs = num_levels;
    
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputRR[lev] = IntVect(D_DECL(2, 2, 2));
    }

    const std::string& pltfile = amrex::Concatenate("plt", nstep, 5);    
    WriteMultiLevelPlotfile(pltfile, output_levs, GetVecOfConstPtrs(output_cc),
                            varnames, geom, 0.0, level_steps, outputRR);

    pc.Checkpoint(pltfile, "particle0", true, particle_varnames);

}

void computeE(      VectorMeshData& E,
              const ScalarMeshData& phi, 
              const Vector<Geometry>& geom) {

    const int num_levels = E.size();

    for (int lev = 0; lev < num_levels; ++lev) {

        const auto& gm = geom[lev];
        const Real* dx = gm.CellSize();

        for (MFIter mfi(*phi[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            compute_E_nodal(bx.loVect(), bx.hiVect(),
                            (*phi[lev] )[mfi].dataPtr(),
                            (*E[lev][0])[mfi].dataPtr(),
                            (*E[lev][1])[mfi].dataPtr(),
#if BL_SPACEDIM == 3
                            (*E[lev][2])[mfi].dataPtr(),
#endif
                            dx);
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
#if BL_SPACEDIM == 3
        E[lev][2]->FillBoundary(gm.periodicity());                 
#endif
        //        VisMF::Write(*E[lev][0], amrex::MultiFabFileFullPrefix(lev, "tmp", "Level_", "Ex"));
    }
}
        
void zeroOutBoundary(MultiFab& input_data,
                     MultiFab& bndry_data,
                     const FabArray<BaseFab<int> >& mask) {
    bndry_data.setVal(0.0, 1);
    for (MFIter mfi(input_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        zero_out_bndry(bx.loVect(), bx.hiVect(),
                       input_data[mfi].dataPtr(),
                       bndry_data[mfi].dataPtr(),
                       mask[mfi].dataPtr());
    }
    bndry_data.FillBoundary();
}

void sumFineToCrseNodal(const MultiFab& fine, MultiFab& crse, 
                        const Geometry& cgeom, const IntVect& ratio) {
    
    const BoxArray& fine_BA = fine.boxArray();
    const DistributionMapping& fine_dm = fine.DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(ratio);
    
    MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
    coarsened_fine_data.setVal(0.0);
    
    for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Box& crse_box = coarsened_fine_data[mfi].box();
        const Box& fine_box = fine[mfi].box();
        sum_fine_to_crse_nodal(bx.loVect(), bx.hiVect(), ratio.getVect(),
                               coarsened_fine_data[mfi].dataPtr(), crse_box.loVect(), crse_box.hiVect(),
                               fine[mfi].dataPtr(), fine_box.loVect(), fine_box.hiVect());
    }
    
    crse.copy(coarsened_fine_data, cgeom.periodicity(), FabArrayBase::ADD);    
}

void fixRHSForSolve(Vector<std::unique_ptr<MultiFab> >& rhs,
                    const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
                    const Vector<Geometry>& geom, const IntVect& ratio) {
    int num_levels = rhs.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        MultiFab& fine_rhs = *rhs[lev];
        const FabArray<BaseFab<int> >& mask = *masks[lev];        
        const BoxArray& fine_ba = fine_rhs.boxArray();
        const DistributionMapping& fine_dm = fine_rhs.DistributionMap();
        MultiFab fine_bndry_data(fine_ba, fine_dm, 1, 1);
        zeroOutBoundary(fine_rhs, fine_bndry_data, mask);
    }
}

void computePhi(ScalarMeshData& rhs, ScalarMeshData& phi,
                Vector<BoxArray>& grids,
                Vector<DistributionMapping>& dm,
                Vector<Geometry>& geom,
                Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks) {

    int num_levels = rhs.size();

    Vector<std::unique_ptr<MultiFab> > tmp_rhs(num_levels);    
    for (int lev = 0; lev < num_levels; ++lev) {
        tmp_rhs[lev].reset(new MultiFab(rhs[lev]->boxArray(), dm[lev], 1, 0));
        MultiFab::Copy(*tmp_rhs[lev], *rhs[lev], 0, 0, 1, 0);
    }
    
    IntVect ratio(D_DECL(2, 2, 2));
    fixRHSForSolve(tmp_rhs, masks, geom, ratio);

    int verbose = 2;
    Real rel_tol = 1.0e-12;
    Real abs_tol = 1.0e-14;

    Vector<Geometry>            level_geom(1);
    Vector<BoxArray>            level_grids(1);
    Vector<DistributionMapping> level_dm(1);
    Vector<MultiFab*>           level_phi(1);
    Vector<const MultiFab*>     level_rhs(1);
    
    for (int lev = 0; lev < num_levels; ++lev) {
        level_phi[0]   = phi[lev].get();
        level_rhs[0]   = tmp_rhs[lev].get();
        level_geom[0]  = geom[lev];
        level_grids[0] = grids[lev];
        level_dm[0]    = dm[lev];

        MLNodeLaplacian linop(level_geom, level_grids, level_dm);

        linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet)});

        linop.setLevelBC(0, nullptr);

        MultiFab sigma(level_grids[0], level_dm[0], 1, 0);
        sigma.setVal(1.0);
        linop.setSigma(0, sigma);
        
        MLMG mlmg(linop);
        mlmg.setMaxIter(100);
        mlmg.setMaxFmgIter(0);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(0);

        mlmg.solve(level_phi, level_rhs, rel_tol, abs_tol);

        if (lev < num_levels-1) {

            NoOpPhysBC cphysbc, fphysbc;
#if BL_SPACEDIM == 3
            int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
            int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
#else
            int lo_bc[] = {INT_DIR, INT_DIR};
            int hi_bc[] = {INT_DIR, INT_DIR};
#endif
            Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
            NodeBilinear mapper;

            amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0, *phi[lev],
                                         0, 0, 1, geom[lev], geom[lev+1],
                                         cphysbc, 0, fphysbc, 0,
                                         IntVect(D_DECL(2, 2, 2)), &mapper, bcs, 0);
        }

        //        VisMF::Write(*phi[lev], amrex::MultiFabFileFullPrefix(lev, "tmp", "Level_", "phi"));
        //        VisMF::Write(*rhs[lev], amrex::MultiFabFileFullPrefix(lev, "tmp", "Level_", "rhs"));
    }
    
    for (int lev = 0; lev < num_levels; ++lev) {
        const Geometry& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }    
}

void
getLevelMasks(Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
              const Vector<BoxArray>& grids,
              const Vector<DistributionMapping>& dmap,
              const Vector<Geometry>& geom,
              const int ncells = 1) {
    int num_levels = grids.size();
    BL_ASSERT(num_levels == dmap.size());

    int covered = 0;
    int notcovered = 1;
    int physbnd = 1;
    int interior = 0;
    
    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();
        
        FabArray<BaseFab<int> > tmp_mask(nba, dmap[lev], 1, ncells);
        tmp_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                           covered, notcovered, physbnd, interior);
        masks[lev].reset(new FabArray<BaseFab<int> >(nba, dmap[lev], 1, 0));
        for (MFIter mfi(tmp_mask); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            build_mask(bx.loVect(), bx.hiVect(),
                       tmp_mask[mfi].dataPtr(), (*masks[lev])[mfi].dataPtr(), &ncells);
        }
    }
}

void main_main () {
    int max_level, n_cell, max_grid_size, particle_output_int, n_buffer, max_step, is_periodic[BL_SPACEDIM];
    Real dt;

    // inputs parameters
    {
        ParmParse pp;
        pp.get("max_level", max_level);
        pp.get("n_cell", n_cell);
        pp.get("n_buffer", n_buffer);
        pp.get("max_grid_size", max_grid_size);       
        pp.get("max_step", max_step);
        pp.get("particle_output_int", particle_output_int);
        pp.get("dt", dt);
    }
    
    assert(max_level < 2);

    int num_levels = max_level + 1;
    
    Vector<int> rr(num_levels-1);
    for (int lev = 1; lev < num_levels; lev++)
        rr[lev-1] = 2;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,-20.0e-6);
        real_box.setHi(n, 20.0e-6);
    }
    
    // This sets the boundary conditions to be doubly or triply periodic
    for (int i = 0; i < BL_SPACEDIM; i++) {
        is_periodic[i] = 0;
    }

    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    Box domain(dom_lo, dom_hi);
   
    // make Geometry for each level
    Vector<Geometry> geom(num_levels);
    geom[0].define(domain,&real_box,CoordSys::cartesian,is_periodic);
    for (int lev = 1; lev < num_levels; ++lev) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, CoordSys::cartesian, is_periodic);
    }
    
    // make grids for each level
    Vector<BoxArray> grids(num_levels);
    grids[0].define(domain);
    if (num_levels > 1) {
        int n_fine = n_cell*rr[0];
        IntVect refined_lo(D_DECL(3*n_fine/8,3*n_fine/8,3*n_fine/8)); 
        IntVect refined_hi(D_DECL(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[1].define(refined_patch);
    }
    
    for (int lev = 0; lev < num_levels; lev++) {
        grids[lev].maxSize(max_grid_size);
    }
    
    int Ncomp  = 1;
    
    Vector<DistributionMapping> dm(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    Vector<std::unique_ptr<MultiFab> > rhs(num_levels);
    Vector<std::array<std::unique_ptr<MultiFab>, BL_SPACEDIM> > eField(num_levels);

    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();
        dm[lev].define(grids[lev]);

        rhs[lev].reset(new MultiFab(nba, dm[lev], Ncomp, 1));
        phi[lev].reset(new MultiFab(nba, dm[lev], Ncomp, 2));
        
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            eField[lev][idim].reset(new MultiFab(nba, dm[lev], Ncomp, 1));
            eField[lev][idim]->setVal(0.0);
        }

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);        
    }

    Vector<std::unique_ptr<FabArray<BaseFab<int> > > > masks(num_levels); 
    getLevelMasks(masks, grids, dm, geom);

    Vector<std::unique_ptr<FabArray<BaseFab<int> > > > gather_masks(num_levels);
    getLevelMasks(gather_masks, grids, dm, geom, n_buffer + 1); // convert from num nodes to num cells
    
    ElectrostaticParticleContainer myPC(geom, dm, grids, rr);

    myPC.InitParticles();

    for (int step = 0; step <= max_step; ++step) {

        myPC.DepositCharge(rhs);

        computePhi(rhs, phi, grids, dm, geom, masks);
        
        computeE(eField, phi, geom);

        myPC.FieldGather(eField, gather_masks);

        if (step % particle_output_int == 0) myPC.writeParticles(step);

        //        WritePlotFile(rhs, phi, eField, myPC, geom, step);
        
        myPC.Evolve(eField, rhs, dt);

        myPC.Redistribute();        
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    main_main();

    amrex::Finalize();
}
