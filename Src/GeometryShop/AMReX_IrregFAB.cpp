#include <cmath>
#include <cstdlib>

#include "AMReX_IrregFAB.H"

namespace amrex
{
///
  IrregFAB::
  IrregFAB():BaseIVFAB<Real>()
  {
  }

///
  IrregFAB::
  IrregFAB(const Box&        a_region,
           const EBGraph&    a_ebgraph,
           const int&        a_nvarin):BaseIVFAB<Real>()
  {
    define(a_region, a_ebgraph, a_nvarin);
  }

///
  IrregFAB::
  ~IrregFAB()
  {
  }

///
  void
  IrregFAB::
  define(const Box&        a_region,
         const EBGraph&    a_ebgraph,
         const int&        a_nvarin)
  {
    m_region = a_region;
    IntVectSet ivs = a_ebgraph.getIrregCells(a_region);
    BaseIVFAB<Real>::define(ivs, a_ebgraph, a_nvarin);
  }
/**********************/
/**********************/
  IrregFAB&
  IrregFAB::
  operator+=(const IrregFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);

    plus(a_src, 0, 0, m_nComp);

    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  plus(const IrregFAB& a_src,
       int a_srccomp,
       int a_destcomp,
       int a_numcomp)
  {
    Box locRegion = a_src.m_region & m_region;
    plus(a_src, locRegion, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  applyOp(const IrregFAB& a_src,
          const Box& a_region,
          int a_srccomp,
          int a_destcomp,
          int a_numcomp,
          IrregFAB::arithOp& a_op)
  {
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.m_nComp);
    BL_ASSERT(a_destcomp + a_numcomp <= m_nComp);
    const Box& locRegion = a_region;

    if (!locRegion.isEmpty())
    {
      IntVectSet ivsMulti = a_src.getIVS();
      ivsMulti &= getIVS();
      ivsMulti &= locRegion;
      for(int ivof = 0; ivof < m_vofs.size(); ivof++)
      {
        const VolIndex& vof = m_vofs[ivof];
        const IntVect& iv = vof.gridIndex();
        if(ivsMulti.contains(iv))
        {
          for (int icomp = 0; icomp < a_numcomp; ++icomp)
          {
            a_op.func((*this)(vof, a_destcomp+icomp),
                      a_src(  vof, a_srccomp+icomp));
          }
        }
      }
    }
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  applyOp(const Real& a_src,
          int a_srccomp,
          int a_destcomp,
          int a_numcomp,
          IrregFAB::arithOp& a_op)
  {
    BL_ASSERT(a_destcomp + a_numcomp <= m_nComp);

    Real* l = dataPtr(a_destcomp);
    int nvof = numVoFs();

    for (int i=0; i<a_numcomp*nvof; i++)
    {
      a_op.func(l[i], a_src);
    }

    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  plus(const IrregFAB& a_src,
       const Box& a_region,
       int a_srccomp,
       int a_destcomp,
       int a_numcomp)
  {
    IrregFAB::additionOp op;
    this->applyOp(a_src, a_region, a_srccomp, a_destcomp, a_numcomp, op);
    return *this;
  }
/**********************/
  IrregFAB&
  IrregFAB::
  operator-=(const IrregFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);

    minus(a_src, 0, 0, m_nComp);

    return *this;
  }
/**********************/
  IrregFAB&
  IrregFAB::
  minus(const IrregFAB& a_src,
        int a_srccomp,
        int a_destcomp,
        int a_numcomp)
  {
    minus(a_src, m_region, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }
/**********************/
  IrregFAB&
  IrregFAB::
  minus(const IrregFAB& a_src,
        const Box& a_region,
        int a_srccomp,
        int a_destcomp,
        int a_numcomp)
  {
    IrregFAB::subtractionOp op;
    applyOp(a_src, a_region, a_srccomp, a_destcomp, a_numcomp, op);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  operator*=(const IrregFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);

    mult(a_src, 0, 0, m_nComp);

    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  mult(const IrregFAB& a_src,
       const Box& a_region,
       int a_srccomp,
       int a_destcomp,
       int a_numcomp)
  {
    IrregFAB::multiplicationOp op;
    applyOp(a_src, a_region, a_srccomp, a_destcomp, a_numcomp, op);
    return *this;
  }
/**********************/
  IrregFAB&
  IrregFAB::
  mult(const IrregFAB& a_src,
       int a_srccomp,
       int a_destcomp,
       int a_numcomp)
  {
    mult(a_src, m_region, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::
  operator/=(const IrregFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);

    divide(a_src, 0, 0, m_nComp);

    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::divide(const IrregFAB& a_src,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp)
  {
    Box region = m_region;
    region &= a_src.m_region;
    divide(a_src, region, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }
/**********************/
  IrregFAB&
  IrregFAB::
  divide(const IrregFAB& a_src,
         const Box& a_region,
         int a_srccomp,
         int a_destcomp,
         int a_numcomp)
  {
    IrregFAB::divisionOp op;
    applyOp(a_src, a_region, a_srccomp, a_destcomp, a_numcomp, op);
    return *this;
  }

/**********************/
/**********************/
  IrregFAB&
  IrregFAB::operator+=(const Real& a_src)
  {

    plus(a_src, 0, 0, m_nComp);
    return *this;
  }

/**********************/
/**********************/
  IrregFAB&
  IrregFAB::operator-=(const Real& a_src)
  {
    minus(a_src, 0, 0, m_nComp);
    return *this;
  }

/**********************/
/**********************/
  IrregFAB&
  IrregFAB::operator*=(const Real& a_src)
  {
    mult(a_src, 0, 0, m_nComp);
    return *this;
  }

/**********************/

  IrregFAB&
  IrregFAB::operator/=(const Real& a_src)
  {
    divide(a_src, 0, 0, m_nComp);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::plus(const Real& a_src, int isrc,int idst, int inco)
  {
    IrregFAB::additionOp op;
    applyOp(a_src, isrc, idst, inco, op);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::minus(const Real& a_src, int isrc,int idst, int inco)
  {
    IrregFAB::subtractionOp op;
    applyOp(a_src, isrc, idst, inco, op);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::mult(const Real& a_src, int isrc,int idst, int inco)
  {
    IrregFAB::multiplicationOp op;
    applyOp(a_src, isrc, idst, inco, op);
    return *this;
  }

/**********************/
  IrregFAB&
  IrregFAB::divide(const Real& a_src, int isrc,int idst, int inco)
  {
    IrregFAB::divisionOp op;
    applyOp(a_src, isrc, idst, inco, op);
    return *this;
  }

/**********************/
}

