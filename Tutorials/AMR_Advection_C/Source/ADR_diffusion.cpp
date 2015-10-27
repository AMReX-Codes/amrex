#include <winstd.H>

#include "ADR.H"
#include "ADR_F.H"
#include "Diffusion.H"

using std::string;

void
ADR::add_diffusion_to_old_source (MultiFab& ext_src_old, MultiFab& OldDiffTerm, Real prev_time, int comp)
{
    getDiffusionTerm(prev_time,OldDiffTerm,comp);
    MultiFab::Add(ext_src_old,OldDiffTerm,comp,comp,1,0);
    geom.FillPeriodicBoundary(ext_src_old,comp,1);
}

void
ADR::time_center_diffusion(MultiFab& S_new, MultiFab& OldDiffTerm, Real cur_time, Real dt, int comp)
{
        // Correct the update so that it will be time-centered.
        MultiFab NewDiffTerm(grids,1,1);
        NewDiffTerm.setVal(0.);
        getDiffusionTerm(cur_time,NewDiffTerm,comp);

        NewDiffTerm.mult( 0.5*dt,   0,1,0);
        OldDiffTerm.mult(-0.5*dt,comp,1,0);

        // Subtract off half of the old source term, and add half of the new.
        MultiFab::Add(S_new,OldDiffTerm,comp,comp,1,0);
        MultiFab::Add(S_new,NewDiffTerm,   0,comp,1,0);
}

void
ADR::getDiffusionTerm (Real time, MultiFab& DiffTerm, int comp)
{
   MultiFab& S_old = get_old_data(State_Type);
   if (verbose > 1 && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calculating diffusion term of component " << comp << " at level " << level 
                << " at time " << time << std::endl;

   // Fill at this level.
   MultiFab species(grids,1,1,Fab_allocate);
   for (FillPatchIterator fpi    (*this,S_old,1,time,State_Type,   comp,1),
                          fpi_rho(*this,S_old,1,time,State_Type,Density,1);
                          fpi.isValid()&&fpi_rho.isValid();++fpi,++fpi_rho)
   {
       species[fpi].copy(fpi(),fpi().box());
       species[fpi].divide(fpi_rho(),0,0,1);
   }

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   for (int dir = 0; dir < BL_SPACEDIM ; dir++) {
      coeffs.set(dir,new MultiFab);
      BoxArray edge_boxes(grids);
      edge_boxes.surroundingNodes(dir);
      coeffs[dir].define(edge_boxes,1,0,Fab_allocate);
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   for (MFIter mfi(S_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = mfi.validbox();
       BL_FORT_PROC_CALL(FILL_DIFF_COEFF,fill_diff_coeff)
                (bx.loVect(), bx.hiVect(),
                 D_DECL(BL_TO_FORTRAN(coeffs[0][mfi]),
                        BL_TO_FORTRAN(coeffs[1][mfi]),
                        BL_TO_FORTRAN(coeffs[2][mfi])),
                 dx_fine);
   }

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   if (level == 0) {
      diffusion->applyop(species,DiffTerm,coeffs);
   } else if (level > 0) {
      // Fill at next coarser level, if it exists.
      const BoxArray& crse_grids = getLevel(level-1).boxArray();
      MultiFab Crse(crse_grids,1,1,Fab_allocate);
      for (FillPatchIterator fpi    (getLevel(level-1),Crse,1,time,State_Type,   comp,1),
                             fpi_rho(getLevel(level-1),Crse,1,time,State_Type,Density,1);
          fpi.isValid()&&fpi_rho.isValid(); ++fpi,++fpi_rho)
      {
        Crse[fpi].copy(fpi());
        Crse[fpi].divide(fpi_rho(),0,0,1);
      }
      diffusion->applyop(level,species,Crse,DiffTerm,coeffs);
   }

   // Multiply Lap(c) by Rho to make Rho Lap(c)
   for (FillPatchIterator fpi_rho(*this,S_old,1,time,State_Type,Density,1);
                                  fpi_rho.isValid();++fpi_rho)
   {
        DiffTerm[fpi_rho].mult(fpi_rho(),0,0,1);
   }
}
