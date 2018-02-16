/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */
                                                       // 
#include <cmath>

#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_STLIF.H"
#include "AMReX_VisMF.H"
#include "AMReX_EBIndexSpace.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_EBLevelGrid.H"

using namespace amrex;

int makeGeometry(Box& a_domain,
                 Real& a_dx)
{
  int eekflag =  0;
  int maxbox;
  //parse input file

  ParmParse pp;
  pp.get("maxboxsize", maxbox);
  RealVect origin = RealVect::Zero;
  std::vector<int> n_cell;
  std::vector<Real> originvec;
  std::vector<Real> probhivec;
  ParmParse ppg("geometry");
  ppg.getarr("prob_lo", originvec, 0, SpaceDim);
  ppg.getarr("prob_hi", probhivec, 0, SpaceDim);

  pp.getarr("n_cell", n_cell, 0, SpaceDim);

  IntVect lo = IntVect::TheZeroVector();
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
  {
    if (n_cell[ivec] <= 0)
    {
      amrex::Print() << " bogus number of cells input = " << n_cell[ivec];
      return(-1);
    }
    hi[ivec] = n_cell[ivec] - 1;
    origin[ivec] = originvec[ivec];
  }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Real prob_hi = (probhivec[0]- originvec[0]);
  a_dx = prob_hi/n_cell[0];

  string stlfile;
  pp.get("STL_file", stlfile);

  STLIF implicit(stlfile);


  GeometryShop workshop(implicit);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
  ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);

  return eekflag;
}

/****/
void  fillEBMF(MultiFab& a_mf, const EBLevelGrid& a_eblg, const Real& a_dx)
{
  for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
  {
    Box valid = a_eblg.getDBL()[mfi];
    EBISBox ebisBox = a_eblg.getEBISL()[mfi];
    IntVectSet ivsbox(valid);
    Real pi = 4.*atan(1.0);
    a_mf[mfi].setVal(0.0);
    for(VoFIterator vofit(ivsbox, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const IntVect & iv = vofit().gridIndex();
      Real val = 1.0;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        Real x = (iv[idir]+0.5)*a_dx;
        val *= sin(pi*x);
      }
      Real kappa = ebisBox.volFrac(vofit());
      a_mf[mfi](iv, 0) = val;
      a_mf[mfi](iv, 1) = kappa;
    }
  }
}
/****/

int stlgeom()
{
  
  Box domain;
  Real dx;
  //make the initial geometry
  amrex::Print() << "making EBIS" << endl;
  makeGeometry(domain, dx);

  BoxArray ba(domain);
  DistributionMapping dm(ba);
  MultiFab mf(ba, dm, 2, 0);
  amrex::Print() << "filling a multifab with data" << endl;
  EBLevelGrid eblg(ba, dm, domain, 0);
  fillEBMF(mf, eblg, dx);

  string outfile("plt.mf.stl")  ;
  amrex::Print() << "outputting to " << outfile << endl;
  
  Vector<string> varnames(2);
  varnames[0] = string("contrived_sines");
  varnames[1] = string("vfrac");

  Geometry geom(domain);
  Real time = 0; int level_step = 0;  //fake stuff to get fulfill weird API

  WriteSingleLevelPlotfile(outfile, mf, varnames, geom, time, level_step);
  
  return 0;
}


int
main(int argc,char **argv)
{
  amrex::Initialize(argc,argv);

  int eekflag = stlgeom();

  if (eekflag != 0)
    {
      cout << "non zero eek detected = " << eekflag << endl;
      cout << "stl test failed" << endl;
    }
  else
    {
      cout << "stl test passed" << endl;
    }

  return 0;
}


