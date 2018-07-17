#include "AMReX_FaceIndex.H"
#include "AMReX.H"
#include <cassert>
#include <iostream>

namespace amrex
{
  /*****************************************/
  FaceIndex::FaceIndex()
  {
    m_isDefined = false;
    m_isBoundary = false;

  }
  /*****************************************/
  const bool&
  FaceIndex::isDefined() const
  {
    return m_isDefined;
  }
  /******************************/
  const int&
  FaceIndex::direction() const
  {
    return m_direction;
  }

  /*****************************************/
  const bool&
  FaceIndex::isBoundary() const
  {
    return m_isBoundary;
  }
  /*****************************************/
  FaceIndex::FaceIndex(const FaceIndex& a_facein)
  {
    define(a_facein);
  }
  /*****************************************/
  void
  FaceIndex::define(const FaceIndex& a_facein)
  {
    //must be allowed to propogate undefinednitude
    //because of vector<vol>
    m_isDefined = a_facein.m_isDefined;
    m_isBoundary = a_facein.m_isBoundary;
    m_loiv = a_facein.m_loiv;
    m_hiiv = a_facein.m_hiiv;
    m_loIndex = a_facein.m_loIndex;
    m_hiIndex = a_facein.m_hiIndex;
    m_direction = a_facein.m_direction;
  }

  /*****************************************/
  FaceIndex::FaceIndex(const VolIndex& a_vof1,
                       const VolIndex& a_vof2,
                       const int& a_direction)
  {
    define(a_vof1, a_vof2, a_direction);
  }

  /*****************************************/
  FaceIndex::FaceIndex(const VolIndex& a_vof1,
                       const VolIndex& a_vof2)
  {
    define(a_vof1, a_vof2);
  }

  /*****************************************/
  void
  FaceIndex::define(const VolIndex& a_vof1,
                    const VolIndex& a_vof2,
                    const int& a_direction)
  {
    assert((a_direction >= 0) &&
              (a_direction < SpaceDim));
    //one vof or the other has to be in the domain.
    assert((a_vof1.cellIndex() >= 0) || (a_vof2.cellIndex() >= 0));

    //if either vof is outside the domain--shown by a negative
    //cell index, the face is a bounary face
    m_isBoundary =
      ((a_vof1.cellIndex() < 0) || (a_vof2.cellIndex() < 0));
    m_isDefined  = true;

    const IntVect& iv1 = a_vof1.gridIndex();
    const IntVect& iv2 = a_vof2.gridIndex();

    int isign = 1;
    if (iv1[a_direction] < iv2[a_direction])
      {
        m_loiv = iv1;
        m_hiiv = iv2;
        m_loIndex = a_vof1.cellIndex();
        m_hiIndex = a_vof2.cellIndex();
        isign = -1;
      }
    else
      {
        m_loiv = iv2;
        m_hiiv = iv1;
        m_loIndex = a_vof2.cellIndex();
        m_hiIndex = a_vof1.cellIndex();
        isign =  1;
      }

    m_direction = a_direction;
    //this check needs to go away in the land
    //of periodic boundary conditions
    if ((isign*(iv1-iv2)) != BASISV(a_direction))
      {
        std::cerr << "FaceIndex constructor--bad arguments" << std::endl;
        std::cerr << "iv1, iv2, direction =" << iv1 << iv2 << a_direction << std::endl;
        Abort("FaceIndex constructor--bad arguments");
      }
  }
  /*****************************************/
  void
  FaceIndex::define(const VolIndex& a_vof1,
                    const VolIndex& a_vof2)
  {
    const IntVect& iv1 = a_vof1.gridIndex();
    const IntVect& iv2 = a_vof2.gridIndex();

    IntVect div = (iv1 - iv2);
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        div[idir] = std::abs(div[idir]);
      }
    int dir = -1;

    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (div == BASISV(idir))
          {
            dir = idir;
            break;
          }
      }

    if (dir < 0)
      Error("a_vof1 and a_vof2 are not neighbors!");

    define(a_vof1, a_vof2, dir);
  }

  /*****************************************/
  int FaceIndex::faceSign(const VolIndex& a_vof) const
  {
    if (a_vof.gridIndex() == m_loiv)
      return 1;
    else if (a_vof.gridIndex() == m_hiiv)
      return -1;
    else
      return 0;
  }

  /*****************************************/
  const int&
  FaceIndex::cellIndex(const Side::LoHiSide& a_sd) const
  {
    assert(isDefined());
    const int* retval;
    if (a_sd == Side::Lo)
      retval=  &m_loIndex;
    else
      retval=  &m_hiIndex;

    return *retval;
  }
  /*****************************************/
  const IntVect&
  FaceIndex::gridIndex(const Side::LoHiSide& a_sd) const
  {
    assert(isDefined());
    const IntVect* retval;
    if (a_sd == Side::Lo)
      retval=  &m_loiv;
    else
      retval=  &m_hiiv;

    return *retval;
  }
  /*****************************************/
  VolIndex
  FaceIndex::getVoF(const Side::LoHiSide& a_sd) const
  {
    assert(isDefined());
    VolIndex retval;
    if (a_sd == Side::Lo)
      retval=  VolIndex(m_loiv, m_loIndex);
    else
      retval=  VolIndex(m_hiiv, m_hiIndex);

    return retval;
  }
  /*****************************************/
  FaceIndex::~FaceIndex()
  {
  }

  /*****************************************/
  FaceIndex&
  FaceIndex::operator= (const FaceIndex& a_facein)
  {
    if (&a_facein != this)
      {
        //assert(a_facein.isDefined());
        //must be allowed to propogate undefinednitude
        //because of vector<vol>
        m_isDefined = a_facein.m_isDefined;
        m_isBoundary = a_facein.m_isBoundary;
        m_loiv = a_facein.m_loiv;
        m_hiiv = a_facein.m_hiiv;
        m_loIndex = a_facein.m_loIndex;
        m_hiIndex = a_facein.m_hiIndex;
        m_direction = a_facein.m_direction;
      }
    return *this;
  }

  /******************************/
  bool FaceIndex::operator!=(const FaceIndex& a_facein) const
  {
    return (!(*this == a_facein));
  }

  template < >
  int linearSize<FaceIndex>(const FaceIndex& findex)
  {
    int numints = 2*(SpaceDim + 1);
    return numints*sizeof(int);
  }
//FaceIndex specialization of linearIn
  template < >
  void linearIn<FaceIndex>(FaceIndex& a_outputT, const void* const inBuf)
  {
    int* ibuf = (int*)inBuf;
    IntVect ivlo, ivhi;
    int celllo, cellhi;
    celllo = *ibuf; ibuf++;
    cellhi = *ibuf; ibuf++;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      ivlo[idir] = *ibuf; ibuf++;
      ivhi[idir] = *ibuf; ibuf++;
    }
    a_outputT = FaceIndex(VolIndex(ivlo, celllo), VolIndex(ivhi, cellhi));
  }

//FaceIndex specialization of linearOut
  template < >
  void linearOut<FaceIndex>(void* const a_outBuf, const FaceIndex& a_inputT)
  {
    
    int* ibuf = (int*)a_outBuf;
    IntVect ivlo = a_inputT.gridIndex(Side::Lo);
    IntVect ivhi = a_inputT.gridIndex(Side::Hi);
    int celllo   = a_inputT.cellIndex(Side::Lo);
    int cellhi   = a_inputT.cellIndex(Side::Hi);
    *ibuf = celllo; ibuf++;
    *ibuf = cellhi; ibuf++;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      *ibuf = ivlo[idir]; ibuf++;
      *ibuf = ivhi[idir]; ibuf++;
    }
  }

  ostream& operator<<( ostream& out, const FaceIndex& fi )
  {
    out << "[(" << fi.m_loiv << "):" << fi.m_loIndex << ";"
        << " (" << fi.m_hiiv << "):" << fi.m_hiIndex << "]";
    return out;
  }

}
