
#include <cmath>

#include "AMReX_GeometryShop.H"
#include "AMReX_WrappedGShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_BoxArray.H"
#include "AMReX_FabArray.H"
#include "AMReX_EBGraph.H"
namespace amrex
{
////////////
  int checkGraph(const EBGraph& a_ebgraph)
  {
    const Box& domain = a_ebgraph.getDomain();
    Box region = a_ebgraph.getRegion();
    region &= domain;
    ///account for regular and covered cells
    for (BoxIterator bit(region); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      //check to see that there is never a regular cell next
      //to a covered cell
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          IntVect ivshift = BASISV(idir);
          ivshift[idir] *= sign(sit());
          IntVect ivneigh = iv + ivshift;
          if (region.contains(ivneigh))
          {
            if((a_ebgraph.isRegular(iv) && a_ebgraph.isCovered(ivneigh) )  ||
               (a_ebgraph.isCovered(iv) && a_ebgraph.isRegular(ivneigh)))
            {
              return -2;
            }

            if((a_ebgraph.isRegular(iv)      && a_ebgraph.isMultiValued(ivneigh)) ||
               (a_ebgraph.isRegular(ivneigh) && a_ebgraph.isMultiValued(iv)     ))
            {
              return -3;
            }
          }
          std::vector<VolIndex> vofs = a_ebgraph.getVoFs(iv);
          for(int ivof = 0; ivof < vofs.size(); ivof++)
          {
            const VolIndex& vof = vofs[ivof];
            if(vof.gridIndex() != iv)
            {
              return -6;
            }
            if(vof.cellIndex() < 0)
            {
              return -7;
            }
            //check to see if the vof matches the faces
            std::vector<FaceIndex> faces = a_ebgraph.getFaces(vof, idir, sit());
            for(int iface = 0; iface < faces.size(); iface++)
            {
              const FaceIndex& face = faces[iface];
              VolIndex startVoF = face.getVoF(flip(sit()));
              if(startVoF != vof)
              {
                return -8;
              }
              VolIndex neighVoF = face.getVoF(sit());
              if(neighVoF.gridIndex() != ivneigh)
              {
                return -4;
              }
              if(domain.contains(ivneigh) && (neighVoF.cellIndex() < 0))
              {
                return -5;
              }
              if(!domain.contains(ivneigh) && (neighVoF.cellIndex() >= 0))
              {
                return -9;
              }
            }
          }
        }
      }
    }
    return 0;
  }
////////////
  int checkGraph(int a_ishop)
  {
    Real radius = 0.5;
    Real domlen = 1;
    std::vector<Real> centervec(SpaceDim);
    std::vector<int>  ncellsvec(SpaceDim);
    int maxgrid;
    ParmParse pp;
    pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
    pp.get(   "sphere_radius", radius);
    pp.get(   "max_grid_size", maxgrid);
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
    
    IntVect ivlo = IntVect::TheZeroVector();
    IntVect ivhi;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      ivhi[idir] = ncellsvec[idir] - 1;
    }

    Box domain(ivlo, ivhi);
    RealVect origin = RealVect::Zero;
    Real dx = domlen/ncellsvec[0];

    BoxArray ba(domain);
    ba.maxSize(maxgrid);
    DistributionMapping dm(ba);
    
    FabArray<EBGraph> allgraphs(ba, dm, 1, 0);
    pp.get("verbosity", verbosity);
    GeometryShop oldschool(sphere, verbosity);
    WrappedGShop wrappedgs(sphere, verbosity);

    GeometryService* gshop;
    if(a_ishop == 0)
    {
      amrex::Print() << "using GeometryShop for geometry generation" << endl;
      gshop = &oldschool;
    }
    else
    {
      amrex::Print() << "using WrappedGShop for geometry generation" << endl;
      gshop = &wrappedgs;
    }
    
    for(MFIter mfi(allgraphs); mfi.isValid(); ++mfi)
    {
      
      Box validRegion = mfi.validbox();
      Box ghostRegion = validRegion;
      ghostRegion.grow(1);
      ghostRegion &= domain;
      BaseFab<int> regIrregCovered;
      Vector<IrregNode> nodes;
      NodeMap nodemap;
      gshop->fillGraph(regIrregCovered, nodes, nodemap, validRegion, ghostRegion, domain, origin, dx);
      EBGraph& ebgraph = allgraphs[mfi];
      ebgraph.buildGraph(regIrregCovered, nodes, validRegion, domain);
    }
    //check graphs at each box
    for(MFIter mfi(allgraphs); mfi.isValid(); ++mfi)
    {
      EBGraph& ebgraph = allgraphs[mfi];
      int eekflag = checkGraph(ebgraph);
      if(eekflag != 0)
      {
        return 10 + eekflag;
      }
    }
    
    //now lets make a second one and see if copy works
    FabArray<EBGraph> secondgraph(ba, dm, 1, 4);
    secondgraph.copy(allgraphs, 0, 0, 1);
    secondgraph.FillBoundary();
    for(MFIter mfi(allgraphs); mfi.isValid(); ++mfi)
    {
      EBGraph& ebgraph = secondgraph[mfi];
      int eekflag = checkGraph(ebgraph);
      if(eekflag != 0)
      {
        return 20 + eekflag;
      }
    }
    
    //now let's check the answer
    return 0;
  }

}
int
main(int argc,char **argv)
{

  amrex::Initialize(argc,argv);
  {
    for(int ishop = 0; ishop <2; ishop++)
    {
      int eekflag = amrex::checkGraph(ishop);

      if (eekflag != 0)
      {
        cout << "non zero eek detected = " << eekflag << endl;
        cout << "sphere test failed" << endl;
        return eekflag;
      }
    }


    cout << "sphere test passed" << endl;
  }
  amrex::Finalize();

  return 0;
}


