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
#include "AMReX_FlatPlateGeom.H"

using namespace amrex;
using std::cout;
using std::endl;


int checkFlatPlate()
{
  Real domlen = 1;
  std::vector<int>  ncellsvec(SpaceDim);
  std::vector<Real>  platelovec(SpaceDim);
  std::vector<Real>  platehivec(SpaceDim);

  ParmParse pp;
  int normalDir;
  Real plateLoc;
    
  pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
  pp.getarr(  "plate_lo"     , platelovec, 0, SpaceDim);
  pp.getarr(  "plate_hi"     , platehivec, 0, SpaceDim);
  pp.get(  "plate_location"   , plateLoc);
  pp.get("domain_length", domlen);                     // 
  pp.get("plate_normal" , normalDir);                  // 

  IntVect ivlo = IntVect::TheZeroVector();
  IntVect ivhi;
  RealVect plateLo, plateHi;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      ivhi[idir] = ncellsvec[idir] - 1;
      plateLo[idir] = platelovec[idir];
      plateHi[idir] = platehivec[idir];
    }
  Box domain(ivlo, ivhi);
  Real dx = domlen/ncellsvec[0];

  FlatPlateGeom geometryGenerator(normalDir, plateLoc, plateLo, plateHi);
  BaseFab<int> regIrregCovered(domain,1);
  Vector<IrregNode> nodes;
  RealVect origin = RealVect::Zero;
  NodeMap nodemap;
  geometryGenerator.fillGraph(regIrregCovered, nodes, nodemap, domain, domain, domain, origin, dx);
  for(int inode = 0; inode < nodes.size(); inode++)
  {
    std::cout << "(" << inode << "): ";
    const IrregNode& node = nodes[inode];
    std::cout << "vof = (" << node.m_cell << "," << node.m_cellIndex << "), volfrac = " << node.m_volFrac << endl;
  }

  return 0;
}
int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  const char* in_file = "flatplate.inputs";
  //parse input file
  ParmParse pp;
  pp.Initialize(0, NULL, in_file);

  // check volume and surface area of approximate sphere
  int eekflag = checkFlatPlate();

  if (eekflag != 0)
    {
      cout << "non zero eek detected = " << eekflag << endl;
      cout << "plate test failed" << endl;
    }
  else
    {
      cout << "plate test passed" << endl;
    }

  return 0;
}


