#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#include "IrregNode.H"
#include "EB_TYPEDEFS.H"
#include "NamespaceHeader.H"

/*******************************/
/*******************************/
IrregNode::IrregNode()
{
}
/*******************************/
/*******************************/
IrregNode::~IrregNode()
{
}


std::ostream& operator<< (std::ostream&  a_os,
                          const IrregNode& a_iv)
{
  a_os<<a_iv.m_cell<<" index:"<<a_iv.m_cellIndex<<" volFrac:"<<a_iv.m_volFrac
      <<" centroid:"<<a_iv.m_volCentroid<<" m_bndryCentroid:"<<a_iv.m_bndryCentroid<<"\n"
      <<"x-arcsLo:"<<a_iv.m_arc[0]<<" x-areaFracsLo:"<<a_iv.m_areaFrac[0]<<"\n"
      <<"y-arcsLo:"<<a_iv.m_arc[1]<<" y-areaFracsLo:"<<a_iv.m_areaFrac[1]<<"\n"
      <<"x-arcsHi:"<<a_iv.m_arc[2]<<" x-areaFracsHi:"<<a_iv.m_areaFrac[2]<<"\n"
      <<"y-arcsHi:"<<a_iv.m_arc[3]<<" y-areaFracsHi:"<<a_iv.m_areaFrac[3]<<"\n";
  return a_os;

  return a_os;
}
/*******************************/
/*******************************/
int IrregNode::
index(int a_idir, Side::LoHiSide a_sd)
{
  CH_assert(a_idir >= 0 && a_idir < SpaceDim);
  int retval;
  if (a_sd == Side::Lo)
    {
      retval = a_idir;
    }
  else
    {
      retval = a_idir + SpaceDim;
    }
  return retval;
}
/*******************************/
/*******************************/

void IrregNode::makeRegular(const IntVect& iv, const Real& a_dx)
{
  setMomentsToRegular(a_dx);
  m_cell = iv;
  m_volFrac = 1.0;
  m_cellIndex = 0;
  m_volCentroid = IntVect::Zero;
  m_bndryCentroid = IntVect::Zero;
  //low sides
  for (int i=0; i<SpaceDim; i++)
    {
      m_arc[i].resize(1,0);
      m_areaFrac[i].resize(1,1.0);
      RealVect faceCenter = IntVect::Zero;
      faceCenter[i] = -0.5;
      m_faceCentroid[i].resize(1,faceCenter);
    }
  //hi sides
  for (int i=0; i<SpaceDim; i++)
    {
      m_arc[i+SpaceDim].resize(1,0);
      m_areaFrac[i+SpaceDim].resize(1,1.0);
      RealVect faceCenter = IntVect::Zero;
      faceCenter[i] = 0.5;
      m_faceCentroid[i+SpaceDim].resize(1,faceCenter);
    }
}

void IrregNode::faceReserve(int location, int size)
{
  if (m_arc[location].size() < size)
    {
      m_arc[location].resize(size);
      m_areaFrac[location].resize(size);
      m_faceCentroid[location].resize(size);
    }
}

/*******************************/
void 
IrregNode::
setMomentsToZero()
{
  m_volumeMoments.setToZero();
  m_EBMoments.setToZero();
  for(int idir  = 0; idir < SpaceDim; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
  for(int iface  = 0; iface < 2*SpaceDim; iface++)
    {
      m_faceMoments[iface].setToZero();
    }
}
/*******************************/
void 
IrregNode::
setMomentsToRegular(const Real& a_dx)
{
  //no eb for regular
  m_EBMoments.setToZero();
  for(int idir  = 0; idir < SpaceDim; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
  m_volumeMoments.setRegular(a_dx);
  for(int iface  = 0; iface < 2*SpaceDim; iface++)
    {
      m_faceMoments[iface].setRegular(a_dx);
    }
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
/*******************************/

#include "NamespaceFooter.H"
