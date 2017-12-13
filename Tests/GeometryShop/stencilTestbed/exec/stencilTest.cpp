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

#include "lapl_nd_F.H"
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ArrayLim.H>
#include <map>
#include <fstream>
#include "stencilTest_F.H"


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
  // Assumes data outside domain lives at cell center as is correct
  a_stencil.clear();
  Real dxinvsq = 1.0/a_dx/a_dx;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    IntVect ivlo = a_iv - BASISV(idir); 
    IntVect ivhi = a_iv + BASISV(idir); 
    a_stencil.add(VolIndex(ivlo, 0), dxinvsq);
    a_stencil.add(VolIndex(ivhi, 0), dxinvsq);
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
defineGeometry(BaseFab<int>      & a_regIrregCovered,
               Vector<IrregNode> & a_nodes,
               const Real        & a_radius,
               const RealVect    & a_center,
               const Box         & a_domain,
               const Real        & a_dx)
{
  //inside regular tells whether domain is inside or outside the sphere
  //bool insideRegular = true;
  bool insideRegular = false;
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

struct FaceData
{
  FaceData(const RealVect& a_centroid,
           Real            a_aperature) : m_centroid(a_centroid), m_aperature(a_aperature) {}
  FaceData() {}
  RealVect m_centroid;
  Real m_aperature;
};

void
applyStencilAllFortran(EBCellFAB                       & a_dst,
		       const EBCellFAB                 & a_src,
		       const BaseFab<VoFStencil>       & a_stencil,
		       const BaseFab<int>              & a_regIrregCovered,
		       const std::vector<IrregNode>    & a_nodes,
		       const Box                       & a_domain,
		       const Real                      & a_dx)
{
  BoxArray ba(a_domain);
  DistributionMapping dm(ba);
  MultiFab srcMF(ba,dm,1,1);
  IntVect tilesize(AMREX_D_DECL(10240,8,32));
  const BaseFab<Real> & regSrc = a_src.getSingleValuedFAB();
  BaseFab<Real>       & regDst = a_dst.getSingleValuedFAB();

  int num_tiles = srcMF.getTileArray(tilesize)->tileArray.size();
  std::vector<std::vector<std::map<IntVect,FaceData> > > faceData(SpaceDim);
  std::vector<bool> has_irreg(num_tiles);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    faceData[idir].resize(num_tiles);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
  {
    const Box& tbx = mfi.tilebox();

    // Find all partial faces in this tile
    has_irreg[mfi.tileIndex()] = false;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int n=0; n<a_nodes.size(); ++n)
      {
        const IrregNode& node = a_nodes[n];
        const IntVect& iv = node.m_cell;
        if (tbx.contains(iv))
        {
          for (SideIterator sit; sit.ok(); ++sit)
          {
            int index = IrregNode::index(idir, sit());
            const std::vector<int>& arcs = node.m_arc[index];

            //this excludes boundary and covered faces, as well as those with aperature = 1
            if((arcs.size() > 0) && (arcs[0] >= 0) && (node.m_areaFrac[index][0] < 1))
            {
              int isign = sign(sit());
              IntVect ivface = isign < 0 ? iv : iv + BASISV(idir);
              faceData[idir][mfi.tileIndex()][ivface] = 
                                       FaceData(node.m_faceCentroid[index][0],node.m_areaFrac[index][0]);
              has_irreg[mfi.tileIndex()] = true;
            }
          }
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
  {
    const Box& tbx = mfi.tilebox();
    int i = mfi.tileIndex();
    if (has_irreg[i])
    {
      FArrayBox fd[BL_SPACEDIM];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        fd[idir].resize(amrex::surroundingNodes(tbx,idir),SpaceDim+1); // comps: aperature, centroid[0], centroid[1]
        fd[idir].setVal(0);
        fd[idir].setVal(1,0);

        for (std::map<IntVect,FaceData>::const_iterator it=faceData[idir][i].begin();
             it!=faceData[idir][i].end(); ++it)
        {
          fd[idir](it->first,0) = it->second.m_aperature;
          for (int idir1 = 0; idir1 < SpaceDim; ++idir1)
          {
            fd[idir](it->first,idir1+1) = it->second.m_centroid[idir1];
          }
        }
      }

      lapleb_MSD(BL_TO_FORTRAN_N(regDst,0), 
                 BL_TO_FORTRAN_N(regSrc,0),
                 AMREX_D_DECL(BL_TO_FORTRAN_N(fd[0],0),
                        BL_TO_FORTRAN_N(fd[1],0),
                        BL_TO_FORTRAN_N(fd[2],0)),
                 tbx.loVect(), tbx.hiVect(), &a_dx);
    }
    else
    {
      lapl_MSD(BL_TO_FORTRAN_N(regDst,0), 
               BL_TO_FORTRAN_N(regSrc,0),
               tbx.loVect(), tbx.hiVect(), &a_dx);
    }
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

  Real dx;
  dx = domlen/ncellsvec[0];
  Vector<Real> probLo(SpaceDim,0);

  BaseFab<int> regIrregCovered;
  Vector<IrregNode>  nodes;
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    center[idir] = centervec[idir];
  }

  amrex::Print() << "Define geometry\n";
  BL_PROFILE_VAR("define_geometry",dg);
  defineGeometry(regIrregCovered, nodes, radius, center, domain, dx); 
  BL_PROFILE_VAR_STOP(dg);
  
  EBCellFAB src(grow(domain, 1), 1); //for a ghost cell
  BL_PROFILE_VAR("init_data",init);
  BaseFab<Real>& regSrc = src.getSingleValuedFAB();
  init_phi(BL_TO_FORTRAN_N(regSrc,0),
           src.box().loVect(), src.box().hiVect(),
           &(probLo[0]), &dx);
  BL_PROFILE_VAR_STOP(init);

  EBCellFAB dst(domain, 1);
  BaseFab<VoFStencil> stencils;

  amrex::Print() << "Getting stencils\n";
  BL_PROFILE_VAR("getting_stencils",gs);
  getStencils(stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(gs);

  amrex::Print() << "Pointwise apply everywhere\n";
  BL_PROFILE_VAR("pointwise_apply_everywhere",pae);
  TestbedUtil::applyStencilPointwise(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(pae);

  amrex::Print() << "Fortran + irreg pointwise\n";
  BL_PROFILE_VAR("fortran_plus_irreg_pointwise_apply",fpip);
  TestbedUtil::applyStencilFortranPlusPointwise(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(fpip);

  amrex::Print() << "Aggsten everywhere\n";
  BL_PROFILE_VAR("aggsten_everywhere",ae);
  TestbedUtil::applyStencilAllAggSten(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(ae);

  amrex::Print() << "Fortran + irreg aggsten\n";
  BL_PROFILE_VAR("fortran_plus_aggsten_at_irreg",fpia);
  TestbedUtil::applyStencilFortranPlusAggSten(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(fpia);

  amrex::Print() << "Fortran everywhere\n";
  BL_PROFILE_VAR("fortran_everywhere",fe);
  applyStencilAllFortran(dst, src, stencils, regIrregCovered, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(fe);

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


