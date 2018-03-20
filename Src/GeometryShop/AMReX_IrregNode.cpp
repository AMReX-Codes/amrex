#include "AMReX_IrregNode.H"
#include <cmath>

namespace amrex
{

  Real IrregNode::bndryArea() const
  {
      Real irregArea=0.0;
      for (int idir=0; idir < SpaceDim; idir++)
        {
          int loindex = index(idir, Side::Lo);
          int hiindex = index(idir, Side::Hi);
          Real loArea = 0;   
          Real hiArea = 0;
          const Vector<Real>& loAreas = m_areaFrac[loindex];
          const Vector<Real>& hiAreas = m_areaFrac[hiindex];
          for(int ilo = 0; ilo < loAreas.size(); ilo++)
            {
              loArea += loAreas[ilo];
            }
          for(int ihi = 0; ihi < hiAreas.size(); ihi++)
            {
              hiArea += hiAreas[ihi];
            }
          irregArea += (hiArea-loArea)*(hiArea-loArea);
        }

      return std::sqrt(irregArea);
  }

  void IrregNode::makeRegular(const IntVect& iv, const Real& a_dx, const Box& a_domain)
  {
    m_cell = iv;
    m_volFrac = 1.0;
    m_cellIndex = 0;
    m_volCentroid   = RealVect::Zero;
    m_bndryCentroid = RealVect::Zero;
    //low sides
    for (int i=0; i<SpaceDim; i++)
      {
        IntVect otherIV = iv - BASISV(i);
        if(a_domain.contains(otherIV))
        {
          m_arc[i].resize(1,0);
        }
        else
        {
          m_arc[i].resize(1,-1);
        }
        m_areaFrac[i].resize(1,1.0);
        RealVect faceCenter = RealVect::Zero;
        faceCenter[i] = -0.5;
        m_faceCentroid[i].resize(1,faceCenter);
      }
    //hi sides
    for (int i=0; i<SpaceDim; i++)
      {
        IntVect otherIV = iv + BASISV(i);
        if(a_domain.contains(otherIV))
        {
          m_arc[i+SpaceDim].resize(1,0);
        }
        else
        {
          m_arc[i+SpaceDim].resize(1,-1);
        }
        m_areaFrac[i+SpaceDim].resize(1,1.0);
        RealVect faceCenter = RealVect::Zero;
        faceCenter[i] = 0.5;
        m_faceCentroid[i+SpaceDim].resize(1,faceCenter);
      }


    setMomentsToRegular(a_dx);
  }

/*******************************/
  void 
  IrregNode::
  setNormalizedStuff(const Real& a_dx)
  {
    Real fullCellVolume = D_TERM(a_dx, *a_dx, *a_dx);
    Real fullFaceArea   = D_TERM(1.0,  *a_dx, *a_dx);
    Real volScaleFactor  = 1./fullCellVolume;
    Real areaScaleFactor = 1./fullFaceArea;

    Real volume = m_volumeMoments[IvSpaceDim::Zero];
    m_volFrac = volume*volScaleFactor;
    m_volCentroid   = RealVect::Zero;
    m_bndryCentroid = RealVect::Zero;

    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if(volume > 0.);
      {
        ///volume centroid
        m_volCentroid[idir] = m_volumeMoments[BASISV_TM<int,SpaceDim>(idir)];

        //divide by the volume
        m_volCentroid[idir] /= volume;

        //convert to relative coordinate 
        m_volCentroid[idir] /= a_dx;

      }
      ///boundary centroid
      Real area = m_EBMoments[IvSpaceDim::Zero];
      if(area > 0.) //can be zero if not really cut.
      {
        m_bndryCentroid[idir] = m_EBMoments[BASISV_TM<int,SpaceDim>(idir)];

        m_bndryCentroid[idir] /= area;

        //convert to relative coordinate 
        m_bndryCentroid[idir] /= a_dx;
      }
    }

    for(int ifacedir = 0; ifacedir < SpaceDim; ifacedir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        //areafrac
        int iilist = this->index(ifacedir, sit());
        Real area = m_faceMoments[iilist][IvSDMinOne::Zero];
        m_areaFrac[iilist].resize(1);
        m_areaFrac[iilist][0] = area*areaScaleFactor;

        m_faceCentroid[iilist].resize(1);
        m_faceCentroid[iilist][0] = RealVect::Zero;
        if(area > 0.) //can be zero if there is not really a face
        {
          //face centroids
          int iindex = 0;
          for (int idir = 0; idir < SpaceDim; ++idir)
          {
            if(idir != ifacedir)
            {
              IvSDMinOne mono = BASISV_TM<int,SpaceDim-1>(iindex);
              m_faceCentroid[iilist][0][idir] = m_faceMoments[iilist][mono];
              //normalize  by area
              m_faceCentroid[iilist][0][idir] /= area;
              //normalize by dx
              m_faceCentroid[iilist][0][idir] /= a_dx;
              iindex++;
            }
          }
        }
      }
    }
  }
}
