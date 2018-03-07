#include "AMReX_DirichletConductivityDomainBC.H"
#include "AMReX_EBArith.H"

namespace amrex
{
/*****/
  DirichletConductivityDomainBC::
  DirichletConductivityDomainBC()
  {
  }
/*****/
  DirichletConductivityDomainBC::
  ~DirichletConductivityDomainBC()
  {
  }
/*****/
  void 
  DirichletConductivityDomainBC::
  fillPhiGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               bool           a_homogeneous)
  {
    Box grownBox = a_valid;
    grownBox.grow(1);

    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        Box choppedBox = grownBox;
        choppedBox.grow(idir,-1);
        Box toRegion = EBArith::adjCellBox(choppedBox, idir, sit(), 1);

        if(!m_eblg.getDomain().contains(toRegion))
        {
          for (BoxIterator bit(toRegion); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();
            //fake vof just to get the location
            VolIndex vof(iv, 0);
            RealVect loc = EBArith::getVoFLocation(vof, m_dx, RealVect::Zero);
            int isign = sign(sit());
            IntVect ivneigh = iv - isign*BASISV(idir);
            Real value = bcvaluefunc(loc, idir, sit());
            if(a_homogeneous) value = 0.;

            a_phi(iv, 0) = -a_phi(ivneigh, 0)  + 2.*value;

          }
        }
      } 

    }//end loop over directions
  }

/*****/
  void
  DirichletConductivityDomainBC::
  getFaceFlux(Real&                 a_faceFlux,
              const VolIndex&       a_vof,
              const MFIter    &     a_mfi,
              const EBCellFAB&      a_phi,
              const int&            a_idir,
              const Side::LoHiSide& a_side,
              const bool&           a_useHomogeneous)
  {
    const EBISBox& ebisBox = a_phi.getEBISBox();
    const Box& domainBox = ebisBox.getDomain();
    a_faceFlux = 0.0;
    vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    if (faces.size() > 0)
    {
      if (faces.size()==1)
      {
        IntVectSet cfivs;
        FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                         cfivs,
                                                         ebisBox,
                                                         domainBox);
        for (int isten=0; isten < faceSten.size(); isten++)
        {
          const Real& weight = faceSten.weight(isten);
          const FaceIndex& face = faceSten.face(isten);

          Real thisFaceFlux;
          getFaceGradPhi(thisFaceFlux,face, a_mfi, a_phi, a_idir, a_side, a_useHomogeneous);

          a_faceFlux += thisFaceFlux*weight;
        }
        a_faceFlux *= ebisBox.areaFrac(faces[0]);
      }
      else
      {
        amrex::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
      }
    }

    Real bcoave = 0;
    Real areaTot = 0;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      Real bcoFace  = (*m_bcoef)[a_mfi][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
      areaTot += areaFrac;
    }
    if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
    a_faceFlux *= bcoave;
  }

/*****/
  void
  DirichletConductivityDomainBC::
  getFaceGradPhi(Real&                 a_faceFlux,
                 const FaceIndex&      a_face,
                 const MFIter    &     a_mfi,
                 const EBCellFAB&      a_phi,
                 const int&            a_idir,
                 const Side::LoHiSide& a_side,
                 const bool&           a_useHomogeneous)
  {
    int iside = -sign(a_side);
    const Real ihdx = 2.0 / m_dx;

    Real value = -1.e99;
    if (a_useHomogeneous)
    {
      value = 0.0;
    }
    else if (m_isFunction)
    {
      RealVect point = EBArith::getFaceLocation(a_face,m_dx,m_probLo);
      value = bcvaluefunc(point, a_idir, a_side);
    }
    else
    {
      if (m_onlyHomogeneous)
      {
        value = 0;
      }
      else
      {
        value = m_value;
      }
    }

    const VolIndex& vof = a_face.getVoF(flip(a_side));

    a_faceFlux = iside * ihdx * (a_phi(vof,0) - value);
  }
}
