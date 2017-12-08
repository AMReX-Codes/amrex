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
#include "dirichletTest_F.H"


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

struct EBBndryData
{
  EBBndryData(const RealVect& a_normal,
              const RealVect& a_bndry_centroid,
              Real            a_value)
    : m_normal(a_normal), m_bndry_centroid(a_bndry_centroid), m_value(a_value) {}
  EBBndryData() {}
  RealVect m_normal, m_bndry_centroid;
  Real m_value;
};

void
applyStencilAllFortran(EBCellFAB                       & a_dst,
		       const EBCellFAB                 & a_src,
		       const BaseFab<VoFStencil>       & a_stencil,
		       const BaseFab<int>              & a_regIrregCovered,
		       const std::vector<IrregNode>    & a_nodes,
		       const Box                       & a_domain,
		       Real                              a_dx)
{
  BoxArray ba(a_domain);
  DistributionMapping dm(ba);
  MultiFab srcMF(ba,dm,1,1);
  IntVect tilesize(AMREX_D_DECL(10240,8,32));

  int num_tiles = srcMF.getTileArray(tilesize)->tileArray.size();
  std::vector<std::vector<int> > tileNodes(num_tiles);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
  {
    const Box& tbx = mfi.tilebox();

    for (int n=0; n<a_nodes.size(); ++n)
    {
      const IrregNode& node = a_nodes[n];
      const IntVect& iv = node.m_cell;
      if (tbx.contains(iv))
      {
        tileNodes[mfi.tileIndex()].push_back(n);
      }
    }
  }

  std::vector<std::vector<Stencil> > stencilData(num_tiles);
  std::vector<std::vector<Real> >    stencilBData(num_tiles);
  std::vector<std::vector<IntVect> > stencilBase(num_tiles);

  std::cout << "IntVect is standard layout? " << std::is_standard_layout<IntVect>::value << " size = " << sizeof(IntVect) << '\n';
  std::cout << "Stencil is standard layout? " << std::is_standard_layout<Stencil>::value << " size = " << sizeof(Stencil) <<  '\n';
  std::cout << "RealVect is standard layout? " << std::is_standard_layout<RealVect>::value << " size = " << sizeof(RealVect) <<  '\n';

  std::vector<RealVect> bndry_normals;
  std::vector<RealVect> bndry_centroids;
  std::vector<IntVect> ivs;
  std::vector<Real> bndry_values;
  std::vector<Real> gtmp;

  // Build stencils
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
  {
    int tid = mfi.tileIndex();
    int num_tileNodes = tileNodes[tid].size();
    if (num_tileNodes > 0)
    {
      stencilData[tid].resize(num_tileNodes);
      stencilBData[tid].resize(num_tileNodes);
      stencilBase[tid].resize(num_tileNodes);
      bndry_normals.resize(num_tileNodes);
      bndry_centroids.resize(num_tileNodes);
      ivs.resize(num_tileNodes);
      bndry_values.resize(num_tileNodes);
      for (int n=0; n<num_tileNodes; ++n)
      {
        const IrregNode& node = a_nodes[tileNodes[tid][n]];
	RealVect& normal = bndry_normals[n];
        normal = RealVect(AMREX_D_DECL(0,0,0));
        for (int idir=0; idir<SpaceDim; ++idir)
        {
          for (SideIterator sit; sit.ok(); ++sit)
          {
            int index = IrregNode::index(idir, sit());
            const std::vector<int>& arcs = node.m_arc[index];
            if((arcs.size() > 0) && (arcs[0] >= 0))
            {
              normal[idir] += sign(sit()) * node.m_areaFrac[index][0];
            }
          }
        }
	bndry_centroids[n] = node.m_bndryCentroid;
	ivs[n] = node.m_cell;

        bndry_values[n] = 1;
      }

      get_bndry_grad_stencil(&(stencilData[tid][0]), &(stencilBData[tid][0]), &(stencilBase[tid][0]),
			     &(bndry_normals[0]), &(bndry_centroids[0]), &(ivs[0]), &num_tileNodes, &a_dx);

      gtmp.resize(num_tileNodes);
      const BaseFab<Real> &  rega_src = a_src.getSingleValuedFAB();
      apply_bndry_grad_stencil(&(gtmp[0]), &(bndry_values[0]),
                               BL_TO_FORTRAN_N(rega_src,0),
                               &(stencilData[tid][0]), &(stencilBData[tid][0]), &(stencilBase[tid][0]),
                               &num_tileNodes);
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

  std::vector<Real> dx(SpaceDim);

  for (int idir=0; idir<SpaceDim; ++idir)
  {
    dx[idir] = domlen/ncellsvec[idir];
  }
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
  defineGeometry(regIrregCovered, nodes, radius, center, domain, dx[0]); 
  BL_PROFILE_VAR_STOP(dg);
#if 0  
  EBCellFAB src(grow(domain, 1), 1); //for a ghost cell
  BL_PROFILE_VAR("init_data",init);
  BaseFab<Real> &  regsrc = src.getSingleValuedFAB();
  init_phi(BL_TO_FORTRAN_N(regsrc,0),
           src.box().loVect(), src.box().hiVect(),
           &(probLo[0]), &(dx[0]));
  BL_PROFILE_VAR_STOP(init);

  EBCellFAB dst(domain, 1);
  BaseFab<VoFStencil> stencils;

  amrex::Print() << "Getting stencils\n";
  BL_PROFILE_VAR("getting_stencils",gs);
  getStencils(stencils, regIrregCovered, nodes, domain, dx[0]);
  BL_PROFILE_VAR_STOP(gs);

  amrex::Print() << "Fortran everywhere\n";
  BL_PROFILE_VAR("fortran_everywhere",fe);
  applyStencilAllFortran(dst, src, stencils, regIrregCovered, nodes, domain, dx[0]);
  BL_PROFILE_VAR_STOP(fe);
#endif
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


