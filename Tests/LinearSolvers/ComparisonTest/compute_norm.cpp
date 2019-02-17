#include <iomanip>

#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <COMP_NORM_F.H>

using namespace amrex;

void compute_norm(const Vector<MultiFab*>& soln,
		  const Vector<MultiFab*>& exac, 
		  const Vector<Geometry>& geom,
		  const Vector<BoxArray>& grids,
		  int nsoln, int iCpp, int iHyp)
{
  Vector<Real> twonorm(nsoln, 0.0);
  Vector<Real> maxnorm(nsoln, 0.0);
  Real volume=0.0;
  
  int nlevel = soln.size();
  int ref_ratio = 2.;

  Vector<Vector<Real> > levmaxnorm(nlevel);

  for (int ilev=0; ilev<nlevel; ilev++) {

    levmaxnorm[ilev].resize(nsoln, 0.0);

    BoxArray baf;

    if (ilev < nlevel-1) {
      baf = grids[ilev+1];
      baf.coarsen(ref_ratio);
    }

    for (MFIter mfi(*soln[ilev]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      
      const Box& bx = grids[ilev][i];

      FArrayBox mask(bx, 1);
      mask.setVal(1.0);
      if (ilev < nlevel-1) {
	const std::vector< std::pair<int,Box> > isects = baf.intersections(bx);

        for (int ii = 0; ii < isects.size(); ii++) {
          mask.setVal(0.0, isects[ii].second, 0);
        }
      }

      FArrayBox volbox(bx, 1);
      geom[ilev].GetVolume(volbox, grids[ilev], i, 0);

      BL_FORT_PROC_CALL(LST_COMP_NORM, lst_comp_norm)
	(bx.loVect(), bx.hiVect(),
	 BL_TO_FORTRAN((*soln[ilev])[mfi]),
	 BL_TO_FORTRAN((*exac[ilev])[mfi]),
	 BL_TO_FORTRAN(mask),
	 BL_TO_FORTRAN(volbox),
	 twonorm.dataPtr(),
	 levmaxnorm[ilev].dataPtr(),
	 &volume, &nsoln);
    }
  }

  for (int isol = 0; isol < nsoln; ++isol) {
      for (int ilev = 0; ilev < nlevel; ++ilev) {
          maxnorm[isol] = std::max(maxnorm[isol], levmaxnorm[ilev][isol]);
      }
  }

  ParallelDescriptor::ReduceRealSum(twonorm.dataPtr(), nsoln);
  ParallelDescriptor::ReduceRealSum(volume);
  ParallelDescriptor::ReduceRealMax(maxnorm.dataPtr(), nsoln);

  for (int ilev = 0; ilev < nlevel; ++ilev) {
      ParallelDescriptor::ReduceRealMax(levmaxnorm[ilev].data(), nsoln);
  }

  for (int i=0; i<nsoln; i++) {
      twonorm[i] = std::sqrt(twonorm[i] / volume); 
  }
 
  std::cout << std::setprecision(15);
  if (ParallelDescriptor::IOProcessor()) {
    if (iCpp >= 0) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "BoxLib_C: not implemented yet. " << std::endl;
      //      std::cout << "BoxLib_C: max-norm error = "<< maxnorm[iCpp] << std::endl;
      //      std::cout << "BoxLib_C:   2-norm error = "<< twonorm[iCpp] << std::endl;
    }
    if (iHyp >= 0) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Hypre: max-norm error = "<< maxnorm[iHyp] << std::endl;
      std::cout << "Hypre:   2-norm error = "<< twonorm[iHyp] << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  return;
}

