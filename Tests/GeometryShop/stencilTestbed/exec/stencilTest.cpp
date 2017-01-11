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
#include "AMReX_EBCellFAB.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_TestbedUtil.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

using namespace amrex;
using std::cout;
using std::endl;

////////////
//gets regular grid stencil for kappa*div(grad phi) with neumann eb and domamin bcs.
void getRegularStencil(VoFStencil           & a_stencil,
                       const IntVect        & a_iv,
                       const Box            & a_domain,
                       const Real           & a_dx)
{
  a_stencil.clear();
  Real dxinvsq = 1.0/a_dx/a_dx;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      IntVect ivlo = a_iv - BASISV(idir); 
      IntVect ivhi = a_iv + BASISV(idir); 
      if(a_domain.contains(ivlo))
        {
          a_stencil.add(VolIndex(ivlo, 0), dxinvsq);
        }
      if(a_domain.contains(ivhi))
        {
          a_stencil.add(VolIndex(ivhi, 0), dxinvsq);
        }
    }
  a_stencil.add(VolIndex(a_iv, 0), -2*SpaceDim*dxinvsq);
}



void getFluxStencil(VoFStencil           & a_stencil,
                    const FaceIndex      & a_face,
                    const IrregNode      & a_node,
                    const int            & a_index,
                    const BaseFab<int>   & a_regIrregCovered,
                    const Box            & a_domain,
                    const Real           & a_dx)
{
  
  FaceStencil interpSten = TestbedUtil::getInterpStencil(a_face, 
                                                         a_node, a_index, a_regIrregCovered,
                                                         a_domain, a_dx);

  a_stencil.clear();
  for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      if(!face.isBoundary())
        {
          
          faceCentSten.add(face.getVoF(Side::Hi),  1.0/a_dx, 0);
          faceCentSten.add(face.getVoF(Side::Lo), -1.0/a_dx, 0);
        }
      else
        {
          //no flux at boundaries for this test
        }
      faceCentSten *= weight;
      a_stencil += faceCentSten;
    }

}
////////////
void getIrregularStencil(VoFStencil           & a_stencil,
                         const IrregNode      & a_node,
                         const BaseFab<int>   & a_regIrregCovered,
                         const Box            & a_domain,
                         const Real           & a_dx)
{
  //gets stencil for kappa*div(grad phi) with neumann eb and domamin bcs.
  a_stencil.clear();
  const IntVect& iv = a_node.m_cell;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int index = IrregNode::index(idir, sit());
          const std::vector<int>& arcs = a_node.m_arc[index];
          //this checks for both boundary faces and covered cells
          if((arcs.size() > 0) && (arcs[0] >= 0))
            {
              int isign = sign(sit());
              //assuming no multi-valued stuff here
              IntVect ivneigh = iv + isign*BASISV(idir);
              FaceIndex face(VolIndex(iv, 0), VolIndex(ivneigh, 0));
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, face, a_node, index, a_regIrregCovered, a_domain, a_dx);
              Real areaFrac = a_node.m_areaFrac[index][0];
              fluxStencil *= Real(isign)*areaFrac/a_dx;
              a_stencil += fluxStencil;
            }
        }
    }
}
////////////
void 
defineGeometry(BaseFab<int>           & a_regIrregCovered,
               std::vector<IrregNode> & a_nodes,
               const Real             & a_radius,
               const RealVect         & a_center,
               const Box              & a_domain,
               const Real             & a_dx)
{
  //inside regular tells whether domain is inside or outside the sphere
  bool insideRegular = true;
  SphereIF sphere(a_radius, a_center, insideRegular);
  int verbosity = 0;
  ParmParse pp;
  pp.get("verbosity", verbosity);

  GeometryShop gshop(sphere, verbosity);
  BaseFab<int> regIrregCovered;
  std::vector<IrregNode> nodes;
  Box validRegion = a_domain;
  Box ghostRegion = a_domain; //whole domain so ghosting does not matter
  RealVect origin = RealVect::Zero;
  gshop.fillGraph(a_regIrregCovered, a_nodes, validRegion, ghostRegion, a_domain, origin, a_dx);
}

void getStencils(BaseFab<VoFStencil>          & a_stencil,
                 const BaseFab<int>           & a_regIrregCovered,
                 const std::vector<IrregNode> & a_nodes,
                 const Box                    & a_domain,
                 const Real                   & a_dx)
{
  a_stencil.resize(a_domain, 1);
  for(BoxIterator boxit(a_domain); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      VoFStencil pointSten;
      if(a_regIrregCovered(iv, 0) == 1)
        {
          //regular cells are set to 1
          getRegularStencil(pointSten, iv, a_domain, a_dx);
        }
      else if (a_regIrregCovered(iv, 0) == -1)
        {
          //coveredCells do not get a stencil
        }
      else if(a_regIrregCovered(iv, 0) == 0)
        {
          //remaining cells are irregular--doing those separately.
        }
      else
        {
          Abort("bogus value in regirregcovered");
        }
      a_stencil(iv, 0) = pointSten;
    }

  for(int ivec = 0; ivec < a_nodes.size(); ivec++)
    {
      const IntVect& iv = a_nodes[ivec].m_cell;
      if(a_regIrregCovered(iv, 0) != 0)
        {
          Abort("regirregcovered and nodes inconsistent");
        }
      VoFStencil pointSten;
      getIrregularStencil(pointSten, a_nodes[ivec], a_regIrregCovered, a_domain, a_dx);

      a_stencil(iv, 0) = pointSten;
    }
}

int testStuff()
{
  int eekflag =  0;
  Real radius = 0.5;
  Real domlen = 1;
  std::vector<Real> centervec(SpaceDim);
  std::vector<int>  ncellsvec(SpaceDim);

  ParmParse pp;
  pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
  pp.get(   "sphere_radius", radius);
  pp.getarr("sphere_center", centervec, 0, SpaceDim);
  pp.get("domain_length", domlen);                     // 

  IntVect ivlo = IntVect::TheZeroVector();
  IntVect ivhi;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      ivhi[idir] = ncellsvec[idir] - 1;
    }

  Box domain(ivlo, ivhi);

  Real dx = domlen/ncellsvec[0];

  BaseFab<int> regIrregCovered;
  std::vector<IrregNode>  nodes;
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centervec[idir];
    }
  defineGeometry(regIrregCovered, nodes, radius, center, domain, dx); 

  
  EBCellFAB src(domain, 1);
  src.setVal(0.0);
  EBCellFAB dst(domain, 1);
  BaseFab<VoFStencil> stencils;
  {
    BL_PROFILE("getting_stencils");
    getStencils(stencils, regIrregCovered, nodes, domain, dx);
  }

  {
    BL_PROFILE("pointwise_apply_everywhere");
    TestbedUtil::applyStencilPointwise(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  }
  {
    BL_PROFILE("fortran_plus_irreg_pointwise_apply");
    TestbedUtil::applyStencilFortranPlusPointwise(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  }
  {
    BL_PROFILE("aggsten_everywhere");
    TestbedUtil::applyStencilAllAggSten(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  }
  {
    BL_PROFILE("fortran_plus_aggsten_at_irreg");
    TestbedUtil::applyStencilFortranPlusAggSten(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  }
  return eekflag;
}


int
main(int argc,char **argv)
{
  amrex::Initialize(argc,argv);
  {
    BL_PROFILE_VAR("main()", pmain);
                                      
    int eekflag = testStuff();

    if (eekflag != 0)
      {
        cout << "non zero eek detected = " << eekflag << endl;
        cout << "sphere test failed" << endl;
      }
    else
      {
        cout << "stencil test passed" << endl;
      }

    BL_PROFILE_VAR_STOP(pmain);
  }
  amrex::Finalize();
  return 0;
}


