////test ported from chombo
#include <cmath>
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_MultiFab.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_BoxIterator.H"

using namespace amrex;


Real
getQuadValue(const IntVect& index, Real dx, int idir)
{
  Real xh = dx*(index[idir] + 1);
  Real xl = dx*(index[idir]    );
  //  Real vol = D_TERM(dx, *dx, *dx);
  Real aveXsq = (xh*xh*xh*xh - xl*xl*xl*xl)/(4.*dx);
  return aveXsq;
}

Real
getQuadValueMD(const IntVect& index, Real dx)
{
  Real xh = dx*(index[0] + 1);
  Real xl = dx*(index[0]    );
  Real yh = dx*(index[1] + 1);
  Real yl = dx*(index[1]    );
  //  Real vol = D_TERM(dx, *dx, *dx);
  Real vol = dx;
  Real aveXsq = (xh*xh*xh*xh - xl*xl*xl*xl)/(4.);
  Real aveYsq = (yh*yh*yh*yh - yl*yl*yl*yl)/(4.);
  Real ans = aveXsq + aveYsq;
#if CH_SPACEDIM==3
  Real zh = dx*(index[2] + 1);
  Real zl = dx*(index[2]    );
  Real aveZsq = (zh*zh*zh*zh - zl*zl*zl*zl)/(4.);
  ans += aveZsq;
#endif  
  return ans/vol;
}


int
CFInterpTest()
{
  int retflag = 0;

  int nxc = 32;
  int nxf = 64;
  int nref = 2;
  int nlo = 3*nxf/8;
  int nhi = 5*nxf/8;
  Box domc(IntVect::Zero, (nxc-1)*IntVect::Unit);
  Box domf(IntVect::Zero, (nxf-1)*IntVect::Unit);

  
  Box lev1Box(nlo*IntVect::Unit, (nhi-1)*IntVect::Unit);

  BoxArray bac(domc);
  BoxArray baf(lev1Box);
  DistributionMapping dmc(bac);
  DistributionMapping dmf(baf);

  
  amrex::Print() <<  "dblCoar:" << bac << endl;
  amrex::Print() <<  "dblFine:" << baf << endl;

  int ncomp = 1;
  Real dxf = 3.14159;
  Real dxc = nref*dxf;

  MultiFab phiSave(baf, dmf, 1, 0);

  MultiFab coarData(bac, dmc, 1, 0);
  MultiFab fineData(baf, dmf, 1, 1);
  Real tf = 1.0;
  Vector<Real> timevec(2, tf);
  int isrc = 0; int idst = 0; int inco = ncomp;
  Geometry geomC(domc);
  Geometry geomF(domf);
  BCRec bcs;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    bcs.setLo(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
    bcs.setHi(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
  }
    
  PhysBCFunct cphysbc(geomC,bcs,BndryFunctBase());
  PhysBCFunct fphysbc(geomF,bcs,BndryFunctBase());
  Interpolater* mapper = &quartic_interp;
  IntVect refrat = nref*IntVect::Unit;

  //testing each direction separately
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    //make coarse data = quadratic everywhere
    for (MFIter mfc(coarData); mfc.isValid(); ++mfc)
    {
      FArrayBox& coarFab=coarData[mfc];
      BoxIterator bit(coarFab.box());
      for (bit.reset(); bit.ok(); ++bit)
      {
        coarFab(bit(), 0) = getQuadValue(bit(), dxc, idir);
      }
    }

    //make fine data = quadratic
    //everywhere EXCEPT at ghost cells
    for (MFIter mff(fineData); mff.isValid(); ++mff)
    {
      FArrayBox&  fineFab = fineData[mff];
      fineFab.setVal(0.);
      Box fabIntBox= mff.validbox();
      BoxIterator bit(fabIntBox);
      for (bit.reset(); bit.ok(); ++bit)
      {
        fineFab(bit(), 0) = getQuadValue(bit(), dxf, idir);
      }
    }
    phiSave.copy(fineData);
    Vector<MultiFab*> coarmf(2, &coarData);
    Vector<MultiFab*> finemf(2, &phiSave);

    FillPatchTwoLevels(fineData, tf, 
                       coarmf, timevec,
                       finemf, timevec,
                       isrc, idst, inco,
                       geomC, geomF,
                       cphysbc, fphysbc, refrat,
                       mapper, bcs);


    for (MFIter mff(fineData); mff.isValid(); ++mff)
    {
      FArrayBox&  fineFab = fineData[mff];
      Box iterBox= mff.validbox();
      iterBox.grow(idir,1);
      iterBox &= domf;
      BoxIterator bit(iterBox);
      for (bit.reset(); bit.ok(); ++bit)
      {
        Real rightAns = getQuadValue(bit(), dxf, idir);
        Real fabVal = fineFab(bit(), 0);
        Real diff =  std::abs(fabVal - rightAns) / rightAns ;
        Real eps = 1.0e-4;
#ifdef CH_USE_FLOAT
        eps = sqrt(eps);
#endif
        if (diff > eps)
        {
          amrex::Print() << "1 QCFI failed  at"
                         << "   idir = " << idir
                         << ",  iv   = " << bit() << endl;

          amrex::Print() << "1 RightAns = " << rightAns
                         << ",  fabval = " << fabVal
                         << ",  diff   = " << diff
                         << endl;

          retflag = 1;
        }
        //if (retflag > 0) break;
      }
    }        
  }



  //now for multidim values
  //make coarse data = quadratic everywhere
  for (MFIter mfc(coarData); mfc.isValid(); ++mfc)
  {
    FArrayBox& coarFab=coarData[mfc];
    BoxIterator bit(coarFab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      coarFab(bit(), 0) = getQuadValueMD(bit(), dxc);
    }
  }

  //make fine data = quadratic
  //everywhere EXCEPT at ghost cells
  for (MFIter mff(fineData); mff.isValid(); ++mff)
  {
    FArrayBox&  fineFab = fineData[mff];
    fineFab.setVal(0.);
    Box fabIntBox= mff.validbox();
    BoxIterator bit(fabIntBox);
    for (bit.reset(); bit.ok(); ++bit)
    {
      fineFab(bit(), 0) = getQuadValueMD(bit(), dxf);
    }
  }

  phiSave.copy(fineData);
  Vector<MultiFab*> coarmf(2, &coarData);
  Vector<MultiFab*> finemf(2, &phiSave);

  FillPatchTwoLevels(fineData, tf, 
                     coarmf, timevec,
                     finemf, timevec,
                     isrc, idst, inco,
                     geomC, geomF,
                     cphysbc, fphysbc, refrat,
                     mapper, bcs);

  for (MFIter mff(fineData); mff.isValid(); ++mff)
  {
    FArrayBox&  fineFab = fineData[mff];
    Box iterBox= mff.validbox();
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      iterBox.grow(idir,1);
      iterBox &= domf;
      BoxIterator bit(iterBox);
      for (bit.reset(); bit.ok(); ++bit)
      {
        Real rightAns = getQuadValueMD(bit(), dxf);
        Real fabVal = fineFab(bit(), 0);
        Real diff =  std::abs(fabVal - rightAns) / rightAns ;
        Real eps = 1.0e-4;
#ifdef CH_USE_FLOAT
        eps = sqrt(eps);
#endif
        if (diff > eps)
        {
          amrex::Print() << "2 QCFI failed  at"
                         << "   idir = " << idir
                         << ",  iv   = " << bit() << endl;

          amrex::Print() << "2 RightAns = " << rightAns
                         << ",  fabval = " << fabVal
                         << ",  diff   = " << diff
                         << endl;

          retflag = 2;
        }
        //if (retflag > 0) break;
      }
    }
  }        
  return retflag;
}

int
main(int argc, char *argv[])
{
  amrex::Initialize(argc,argv);

  int eekflag = CFInterpTest();
  if (eekflag != 0)
  {
    amrex::Print() << " fourth order cfi test failed with error code: "
                   << eekflag << endl;
  }
  else
  {
    amrex::Print() <<  " fourth oder cfi test passed. " << endl;
  }
  amrex::Finalize();

  return 0;
}

