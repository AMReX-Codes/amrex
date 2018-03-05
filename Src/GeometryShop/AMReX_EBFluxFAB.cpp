
#include "AMReX_EBFluxFAB.H"

namespace amrex
{

  void
  EBFluxFAB::
  clone(const EBFluxFAB& a_input)
  {
    define(a_input.m_ebisBox, a_input.m_region, a_input.m_nComp);
    this->setVal(0.);
    (*this) += a_input;
  }

  std::size_t 
  EBFluxFAB::
  nBytes (const Box& bx, int start_comp, int ncomps) const
  {
    size_t retval = 0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_fluxes[idir].nBytes(bx, start_comp, ncomps);
    }
    retval += m_irrFlux.nBytes(bx, start_comp, ncomps);
    return retval;
  }

  std::size_t 
  EBFluxFAB::
  copyToMem (const Box& srcbox,
             int        srccomp,
             int        numcomp,
             void*      dst) const
  {
    size_t retval = 0;
    unsigned char* buffer = (unsigned char*) dst;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_fluxes[idir].copyToMem(srcbox, srccomp, numcomp, buffer);
      buffer += retval;
    }
    retval += m_irrFlux.copyToMem(srcbox, srccomp, numcomp, buffer);
    return retval;

  }
// ---------------------------------------------------------
  std::size_t 
  EBFluxFAB::
  copyFromMem (const Box&  dstbox,
               int         dstcomp,
               int         numcomp,
               const void* src)
  {
    size_t retval = 0;
    unsigned char* buffer = (unsigned char*) src;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_fluxes[idir].copyFromMem(dstbox, dstcomp, numcomp, buffer);
      buffer += retval;
    }
    retval += m_irrFlux.copyFromMem(dstbox, dstcomp, numcomp, buffer);
    return retval;
  }

// ---------------------------------------------------------
  int
  EBFluxFAB::
  nComp() const
  {
    return m_nComp;
  }

// ---------------------------------------------------------
  const Box&
  EBFluxFAB::getRegion() const
  {
    return m_region;
  }

// ---------------------------------------------------------
  EBFaceFAB&
  EBFluxFAB::operator[] (const int dir)
  {
    BL_ASSERT(m_nComp >0);
    BL_ASSERT(dir < SpaceDim);

    return m_fluxes[dir];
  }

// ---------------------------------------------------------
  const EBFaceFAB&
  EBFluxFAB::operator[] (const int dir)  const
  {
    BL_ASSERT(m_nComp >0);
    BL_ASSERT(dir < SpaceDim);

    return m_fluxes[dir];
  }

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
  EBFluxFAB::EBFluxFAB()
  {
  }

// ---------------------------------------------------------
  EBFluxFAB::EBFluxFAB(const EBISBox& a_ebisBox,
                       const Box& a_region, int a_nComp)
  {
    define(a_ebisBox, a_region, a_nComp);
  }
  void
  EBFluxFAB::define(const EBISBox& a_ebisBox,
                    const Box& a_region, int a_nComp)
  {
    m_isDefined = true;
    m_ebisBox = a_ebisBox;
    m_region = a_region & a_ebisBox.getRegion();
    m_region &= a_ebisBox.getDomain();
    m_nComp = a_nComp;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_fluxes[idir].define(a_ebisBox, m_region, idir, a_nComp);
    }
    m_irrFlux.define(m_region, a_ebisBox.getEBGraph(), a_nComp);
  }
// ---------------------------------------------------------
  bool
  EBFluxFAB::isDefined() const
  {
    return m_isDefined;
  }
// ---------------------------------------------------------
  EBFluxFAB::~EBFluxFAB()
  {
  }

// ---------------------------------------------------------
  void
  EBFluxFAB::setVal(const Real& val)
  {
    BL_ASSERT(m_nComp > 0);

    for (int dir = 0; dir < SpaceDim; dir++)
    {
      m_fluxes[dir].setVal(val);
    }
    m_irrFlux.setVal(val);
  }

// ---------------------------------------------------------
  EBFluxFAB& 
  EBFluxFAB::
  copy(const EBFluxFAB&  a_src,
       const Box&        a_srcbox,
       int               a_srccomp,
       const Box&        a_dstbox,
       int               a_dstcomp,
       int               a_numcomp)
  {
    for (int idir=0; idir<SpaceDim; idir++)
    {
      m_fluxes[idir].copy(a_src.m_fluxes[idir], a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
    }
    m_irrFlux.copy(a_src.m_irrFlux, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
    return *this;
  }
}
