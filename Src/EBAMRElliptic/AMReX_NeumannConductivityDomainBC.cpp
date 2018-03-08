#include "AMReX_NeumannConductivityDomainBC.H"
#include "AMReX_EBArith.H"
namespace amrex
{
  NeumannConductivityDomainBC::
  NeumannConductivityDomainBC()
  {
  }

  NeumannConductivityDomainBC::
  ~NeumannConductivityDomainBC()
  {
  }

  void 
  NeumannConductivityDomainBC::
  fillPhiGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               bool           a_homogeneous)
  {
    Box grownBox = a_valid;
    grownBox.grow(1);

    for (int idir=0; idir< BL_SPACEDIM; ++idir)
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
            if(a_homogeneous) value = 0;
            a_phi(iv, 0) = a_phi(ivneigh, 0)  + m_dx*value;
          }
        }
      } 
    }//end loop over directions
  }
  /*****/
  void
  NeumannConductivityDomainBC::
  getFaceFlux(Real&                 a_faceFlux,
              const VolIndex&       a_vof,
              const MFIter    &     a_mfi,
              const EBCellFAB&      a_phi,
              const int&            a_idir,
              const Side::LoHiSide& a_side,
              const bool&           a_useHomogeneous)
  {
    const EBISBox& ebisBox = a_phi.getEBISBox();

    Real totalMassFlux = 0.0;
    vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    for (int i = 0; i < faces.size(); i++)
    {
      const RealVect centroid = ebisBox.centroid(faces[i]);
      Real thisFaceFlux;
      getFaceGradPhi(thisFaceFlux,faces[i],a_mfi, a_phi, a_idir, a_side, a_useHomogeneous);
      totalMassFlux += thisFaceFlux;
    }
    if(faces.size() > 1)
    {
      totalMassFlux /= faces.size();
    }
    a_faceFlux = totalMassFlux;
   
    Real bcoave = 0;
    Real areaTot = 0.0;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      areaTot += areaFrac;
      Real bcoFace  = (*m_bcoef)[a_mfi][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
    }
    if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
    a_faceFlux *= bcoave;
  }

  /*****/
  void
  NeumannConductivityDomainBC::
  getFaceGradPhi(Real&                 a_faceFlux,
                 const FaceIndex&      a_face,
                 const MFIter    &     a_mfi,
                 const EBCellFAB&      a_phi,
                 const int&            a_idir,
                 const Side::LoHiSide& a_side,
                 const bool&           a_useHomogeneous)
  {
    const int iside = -sign(a_side);

    Real flux = -1.e99;
    if (a_useHomogeneous)
    {
      flux = 0.0;
    }
    else if (m_isFunction)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point  = EBArith::getFaceLocation(a_face,m_dx,m_probLo);
      point[a_face.direction()] = 0.0;//make this not depend on whatever ebisbox is returning for centroid in the face direction.
      flux = bcvaluefunc(point, a_idir, a_side);
    }
    else
    {
      flux = m_value;
    }
    a_faceFlux = iside*flux;
  }

}
