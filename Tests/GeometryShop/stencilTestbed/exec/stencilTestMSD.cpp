#include <cmath>

#include "AMReX_GeometryShop.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_EBCellFAB.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ArrayLim.H>
#include "stencilTestMSD_F.H"

#include <fstream>

using namespace amrex;
using std::cout;
using std::endl;

void 
defineGeometry(BaseFab<int>           & a_regIrregCovered,
               std::vector<IrregNode> & a_nodes,
               const Real             & a_radius,
               const RealVect         & a_center,
               const Box              & a_domain,
               const Real             & a_dx)
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

struct FaceData
{
  FaceData(const RealVect& a_centroid,
           Real            a_aperature) : m_centroid(a_centroid), m_aperature(a_aperature) {}
  FaceData() {}
  RealVect m_centroid;
  Real m_aperature;

  std::ostream& Write(std::ostream& ofs) const
    {
      ofs.write((char*)(&m_aperature), sizeof(Real));
      ofs.write((char*)(m_centroid.dataPtr()), SpaceDim * sizeof(Real));
      return ofs;
    }
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

  std::ostream& Write(std::ostream& ofs) const
    {
      ofs.write((char*)(&m_value), sizeof(Real));
      ofs.write((char*)(m_normal.dataPtr()), SpaceDim * sizeof(Real));
      ofs.write((char*)(m_bndry_centroid.dataPtr()), SpaceDim * sizeof(Real));
      return ofs;
    }
};

void
applyStencilAllFortran(EBCellFAB                       & a_dst,
		       const EBCellFAB                 & a_src,
		       const std::vector<IrregNode>    & a_nodes,
		       const Box                       & a_domain,
		       const Real                      & a_dx)
{

  BL_PROFILE_VAR_NS("fortran only prep",fp);
  BL_PROFILE_VAR_NS("fortran only sparse-to-dense",fp2);
  BL_PROFILE_VAR_NS("fortran only eb",feb);
  BL_PROFILE_VAR_NS("fortran only reg",fr);

  BL_PROFILE_VAR_START(fp);
  BoxArray ba(a_domain);
  DistributionMapping dm(ba);
  MultiFab srcMF(ba,dm,1,1);
  IntVect tilesize(AMREX_D_DECL(10240,8,32));

  BL_PROFILE_VAR_START(fp);
  int num_tiles = srcMF.getTileArray(tilesize)->tileArray.size();
  std::vector<bool> has_irreg(num_tiles);
  std::vector<bool> has_eb(num_tiles);

  std::vector<std::vector<std::map<IntVect,FaceData> > > faceData(SpaceDim);
  std::vector<std::map<IntVect,EBBndryData> > bndryData;

  Real eb_dirichlet_value = 1;


  bool compute_geom = true;
  ParmParse pp;
  pp.query("compute_geom",compute_geom);
  if (compute_geom)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      faceData[idir].resize(num_tiles);
    }
    bndryData.resize(num_tiles);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
    {
      const Box& tbx = mfi.tilebox();
      int tid = mfi.tileIndex();

      const Box& ovlp = tbx & srcMF.boxArray()[mfi.index()];

      // Find all partial faces in this tile
      has_irreg[tid] = false;
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
                faceData[idir][tid][ivface] = FaceData(node.m_faceCentroid[index][0],node.m_areaFrac[index][0]);
                has_irreg[tid] = true;
              }
            }
          }
        }
      }

      has_eb[tid] = false;
      for (int n=0; n<a_nodes.size(); ++n)
      {
        const IrregNode& node = a_nodes[n];
        const IntVect& iv = node.m_cell;
        if (tbx.contains(iv))
        {
          has_eb[tid] = true;

          RealVect normal(AMREX_D_DECL(0,0,0));
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

          bndryData[tid][iv] = EBBndryData(normal,node.m_bndryCentroid,eb_dirichlet_value);
        }
      }

    }

    // Write faceData to TD files
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      std::string filename=amrex::Concatenate("TD_",idir,1);
      std::ofstream ofs(filename.c_str(),std::ofstream::binary);

      int num_tiles = faceData[idir].size();
      ofs.write((char*)(&num_tiles), sizeof(int));
      for (int i=0; i<faceData[idir].size(); ++i)
      {
        int num_ivs = faceData[idir][i].size();
        ofs.write((char*)(&num_ivs), sizeof(int));
        for (std::map<IntVect,FaceData>::const_iterator it=faceData[idir][i].begin(); it!= faceData[idir][i].end(); ++it)
        {
          ofs.write((char*)(it->first.getVect()), SpaceDim * sizeof(int));
          const FaceData& face_data = it->second;
          ofs.write((char*)(&face_data.m_aperature), sizeof(Real));
          ofs.write((char*)(face_data.m_centroid.dataPtr()), SpaceDim * sizeof(Real));
        }
      }
    }

    // Write eb bndry data to TD file
    std::string filename("TD_C");
    std::ofstream ofs(filename.c_str(),std::ofstream::binary);

    int num_tiles = bndryData.size();
    ofs.write((char*)(&num_tiles), sizeof(int));
    for (int i=0; i<bndryData.size(); ++i)
    {
      int num_ivs = bndryData[i].size();
      ofs.write((char*)(&num_ivs), sizeof(int));
      for (std::map<IntVect,EBBndryData>::const_iterator it=bndryData[i].begin(); it!= bndryData[i].end(); ++it)
      {
        ofs.write((char*)(it->first.getVect()), SpaceDim * sizeof(int));
        const EBBndryData& bndry_data = it->second;
        ofs.write((char*)(&bndry_data.m_value), sizeof(Real));
        ofs.write((char*)(bndry_data.m_normal.dataPtr()), SpaceDim * sizeof(Real));
        ofs.write((char*)(bndry_data.m_bndry_centroid.dataPtr()), SpaceDim * sizeof(Real));
      }
    }
  }
  else
  {
    // Read/reconstruct faceData from TD files
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      std::string filename=amrex::Concatenate("TD_",idir,1);
      std::ifstream ifs(filename.c_str(),std::ofstream::binary);
      if ( ifs.is_open() )
      {
        int num_tiles = -1;
        ifs.read((char*)(&num_tiles), sizeof(int));

        if (num_tiles > 0)
        {
          faceData[idir].resize(num_tiles);
          if (idir == 0)
          {
            has_irreg.resize(num_tiles, false);
          }
        }

        for (int i=0; i<faceData[idir].size(); ++i)
        {
          int num_ivs = -1;
          ifs.read((char*)(&num_ivs), sizeof(int));

          if (num_ivs > 0) has_irreg[i] = true;

          for (int niv=0; niv<num_ivs; ++niv)
          {
            int idxs[SpaceDim];
            ifs.read((char*)(idxs), SpaceDim * sizeof(int));
            IntVect iv(AMREX_D_DECL(idxs[0],idxs[1],idxs[2]));
            
            Real ap;
            ifs.read((char*)(&ap), sizeof(Real));

            Real c[SpaceDim];
            ifs.read((char*)(c), SpaceDim * sizeof(Real));
            RealVect cent(AMREX_D_DECL(c[0],c[1],c[2]));

            faceData[idir][i][iv] = FaceData(cent, ap);
          }
        }
      }
      else
      {
        amrex::Abort("Re-run with compute_geom=t first");
      }
    }

    std::string filename("TD_C");
    std::ifstream ifs(filename.c_str(),std::ofstream::binary);
    if ( ifs.is_open() )
    {
      int num_tiles = -1;
      ifs.read((char*)(&num_tiles), sizeof(int));

      if (num_tiles > 0)
      {
        bndryData.resize(num_tiles);
        has_eb.resize(num_tiles, false);

        for (int i=0; i<bndryData.size(); ++i)
        {
          int num_ivs = -1;
          ifs.read((char*)(&num_ivs), sizeof(int));

          if (num_ivs > 0) has_eb[i] = true;

          for (int niv=0; niv<num_ivs; ++niv)
          {
            int idxs[SpaceDim];
            ifs.read((char*)(idxs), SpaceDim * sizeof(int));
            IntVect iv(AMREX_D_DECL(idxs[0],idxs[1],idxs[2]));
            
            Real val;
            ifs.read((char*)(&val), sizeof(Real));

            Real c[SpaceDim];
            ifs.read((char*)(c), SpaceDim * sizeof(Real));
            RealVect normal(AMREX_D_DECL(c[0],c[1],c[2]));

            ifs.read((char*)(c), SpaceDim * sizeof(Real));
            RealVect centroid(AMREX_D_DECL(c[0],c[1],c[2]));

            bndryData[i][iv] = EBBndryData(normal, centroid, val);
          }
        }
      }
    }
    else
    {
      amrex::Abort("Re-run with compute_geom=t first");
    }
  }
  BL_PROFILE_VAR_STOP(fp);

  int num_iters = 1;
  pp.query("num_iters",num_iters);
  for (int iter=0; iter<num_iters; ++iter)
  { 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(srcMF,tilesize); mfi.isValid(); ++mfi)
    {
      const Box& tbx = mfi.tilebox();
      int i = mfi.tileIndex();
      if (has_irreg[i] || has_eb[i])
      {
        // Build regular version of sparse irregular data
        BL_PROFILE_VAR_START(fp2);
        FArrayBox fd[BL_SPACEDIM], cd;
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          fd[idir].resize(amrex::surroundingNodes(tbx,idir),SpaceDim+1);
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
        cd.resize(tbx,2*SpaceDim+1); // value, normal, centroid
        cd.setVal(0); // normal = 0 => no eb
        for (std::map<IntVect,EBBndryData>::const_iterator it=bndryData[i].begin();
             it!=bndryData[i].end(); ++it)
        {
          cd(it->first,0) = it->second.m_value;
          for (int idir1 = 0; idir1 < SpaceDim; ++idir1)
          {
            cd(it->first,idir1+1) = it->second.m_normal[idir1];
            cd(it->first,idir1+1+SpaceDim) = it->second.m_bndry_centroid[idir1];
          }
        }
        

        BL_PROFILE_VAR_STOP(fp2);

        // Apply irregular operator
        BL_PROFILE_VAR_START(feb);
        lapleb_MSD(BL_TO_FORTRAN_N(a_dst,0), 
                   BL_TO_FORTRAN_N(a_src,0),
                   BL_TO_FORTRAN_N(cd,0),
                   AMREX_D_DECL(BL_TO_FORTRAN_N(fd[0],0),
                          BL_TO_FORTRAN_N(fd[1],0),
                          BL_TO_FORTRAN_N(fd[2],0)),
                   tbx.loVect(), tbx.hiVect(), &a_dx);
        BL_PROFILE_VAR_STOP(feb);
      }
      else
      {
        // Apply regular operator
        BL_PROFILE_VAR_START(fr);
        lapl_MSD(BL_TO_FORTRAN_N(a_dst,0), 
                 BL_TO_FORTRAN_N(a_src,0),
                 tbx.loVect(), tbx.hiVect(), &a_dx);
        BL_PROFILE_VAR_STOP(fr);
      }
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

  Real dx = domlen/ncellsvec[0];
  Vector<Real> probLo(SpaceDim,0);

  BaseFab<int> regIrregCovered;
  std::vector<IrregNode>  nodes;
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    center[idir] = centervec[idir];
  }

  EBCellFAB src(grow(domain, 1), 1); //for a ghost cell
  BL_PROFILE_VAR("init_data",init);
  init_phi(BL_TO_FORTRAN_N(src,0),
           src.box().loVect(), src.box().hiVect(),
           &(probLo[0]), &dx);
  BL_PROFILE_VAR_STOP(init);

  EBCellFAB dst(domain, 1);

  bool compute_geom = true;
  pp.query("compute_geom",compute_geom);
  if (compute_geom) {
    amrex::Print() << "Define geometry\n";
    BL_PROFILE_VAR("define_geometry",dg);
    defineGeometry(regIrregCovered, nodes, radius, center, domain, dx); 
    BL_PROFILE_VAR_STOP(dg);
  }

  amrex::Print() << "Fortran everywhere\n";
  BL_PROFILE_VAR("fortran_everywhere",fe);
  applyStencilAllFortran(dst, src, nodes, domain, dx);
  BL_PROFILE_VAR_STOP(fe);

  if (pp.countval("outfile"))
  {
    std::string outfile;
    pp.get("outfile",outfile);
    std::ofstream ofs(outfile.c_str());
    dst.writeOn(ofs);
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


