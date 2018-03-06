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
  }

}
