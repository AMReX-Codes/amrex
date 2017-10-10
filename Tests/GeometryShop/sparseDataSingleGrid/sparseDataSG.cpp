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

#include <cmath>

#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_EBGraph.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_BaseIFFAB.H"
namespace amrex
{
  int rightAns(const IntVect& a_iv)
  {
    int retval = 1;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval += a_iv[idir];
    }
    return retval;
  }
////////////
  int checkIVFAB(const BaseIVFAB<int> & a_dat1,
                 const BaseIVFAB<int> & a_dat2,
                 const EBGraph        & a_ebgraph)
  {
    const Box& domain = a_ebgraph.getDomain();
    IntVectSet ivs = a_ebgraph.getIrregCells(domain);
    for (VoFIterator vofit(ivs, a_ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      int correct = rightAns(vof.gridIndex());
      int val1 = a_dat1(vof, 0);
      int val2 = a_dat2(vof, 0);
      if(val1 != correct)
      {
        return -1;
      }
      if(val2 != correct)
      {
        return -2;
      }
    }
    return 0;
  }
  void fillIVFAB(BaseIVFAB<int>       & a_dat1,
                 const EBGraph        & a_ebgraph)
  {
    const Box& domain = a_ebgraph.getDomain();
    IntVectSet ivs = a_ebgraph.getIrregCells(domain);
    for (VoFIterator vofit(ivs, a_ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      int correct = rightAns(vof.gridIndex());
      a_dat1(vof, 0) = correct;;
    }
    return ;
  }
  void fillIFFAB(BaseIFFAB<int>       & a_dat1,
                 const EBGraph        & a_ebgraph)
  {
    const Box& domain = a_ebgraph.getDomain();
    IntVectSet ivs = a_ebgraph.getIrregCells(domain);
    int direction = a_dat1.direction();
    for(FaceIterator faceit(ivs, a_ebgraph, direction, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      int correct = rightAns(face.gridIndex(Side::Hi));
      a_dat1(face, 0) = correct;
    }
    return ;
  }
  int checkIFFAB(const BaseIFFAB<int> & a_dat1,
                 const BaseIFFAB<int> & a_dat2,
                 const EBGraph        & a_ebgraph)
  {
    const Box& domain = a_ebgraph.getDomain();
    IntVectSet ivs = a_ebgraph.getIrregCells(domain);
    int direction = a_dat1.direction();
    for(FaceIterator faceit(ivs, a_ebgraph, direction, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      int val1 = a_dat1(face, 0);
      int val2 = a_dat2(face, 0);
      int correct = rightAns(face.gridIndex(Side::Hi));
      if(val1 != correct)
      {
        return -1;
      }
      if(val2 != correct)
      {
        return -2;
      }
    }
    return 0;
  }
////////////
  int checkStuff()
  {
    Real radius = 0.5;
    Real domlen = 1;
    std::vector<Real> centervec(SpaceDim);
    std::vector<int>  ncellsvec(SpaceDim);

    ParmParse pp;
    pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
    pp.get(   "sphere_radius", radius);
    pp.getarr("sphere_center", centervec, 0, SpaceDim);
    pp.get("domain_length", domlen);                     
    RealVect center;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centervec[idir];
    }
    bool insideRegular = false;
    SphereIF sphere(radius, center, insideRegular);
    int verbosity = 0;

    pp.get("verbosity", verbosity);
    GeometryShop gshop(sphere, verbosity);
    BaseFab<int> regIrregCovered;
    Vector<IrregNode> nodes;

    IntVect ivlo = IntVect::TheZeroVector();
    IntVect ivhi;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      ivhi[idir] = ncellsvec[idir] - 1;
    }


    Box domain(ivlo, ivhi);
    Box validRegion = domain;
    Box ghostRegion = domain; //whole domain so ghosting does not matter
    RealVect origin = RealVect::Zero;
    Real dx = domlen/ncellsvec[0];
    gshop.fillGraph(regIrregCovered, nodes, validRegion, ghostRegion, domain, origin, dx);


    EBGraph ebgraph(domain, 1);
    ebgraph.buildGraph(regIrregCovered, nodes, domain, domain);

    IntVectSet ivsirreg = ebgraph.getIrregCells(domain);
    BaseIVFAB<int> ivf1(ivsirreg, ebgraph, 1);
    BaseIVFAB<int> ivf2(ivsirreg, ebgraph, 1);
    fillIVFAB(ivf1, ebgraph);
    ivf2.copy(ivf1, domain, 0, domain, 0, 1);
    checkIVFAB(ivf1, ivf2, ebgraph);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<int> iffab1(ivsirreg, ebgraph, idir, 1);
      BaseIFFAB<int> iffab2(ivsirreg, ebgraph, idir, 1);
      fillIFFAB(iffab1, ebgraph);
      iffab2.copy(iffab1, domain, 0, domain, 0, 1);
      checkIFFAB(iffab1, iffab2, ebgraph);
    }
    return 0;
  }
}
int
main(int argc,char **argv)
{
  amrex::Initialize(argc,argv);
  {
    int eekflag = amrex::checkStuff();

    if (eekflag != 0)
    {
      cout << "non zero eek detected = " << eekflag << endl;
      cout << "sphere test failed" << endl;
    }
    else
    {
      cout << "sphere test passed" << endl;
    }
  }
  amrex::Finalize();

  return 0;
}


