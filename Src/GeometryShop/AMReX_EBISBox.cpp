
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */


#include "AMReX_EBISBox.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_SPMD.H"
#include "AMReX_PolyGeom.H"

namespace amrex
{
/*******************************/
  const Box&
  EBISBox::getRegion() const
  {
    return m_graph.getRegion();
  }
/*******************************/
  EBISBox::EBISBox()
  {
  }
/*******************************/
  EBISBox::~EBISBox()
  {
  }
/*******************************/
  IntVectSet
  EBISBox::getMultiCells(const Box& a_subbox) const
  {
    return m_graph.getMultiCells(a_subbox);
  }
/*******************************/
  IntVectSet
  EBISBox::getIrregIVS(const Box& a_subbox) const
  {
    return m_graph.getIrregCells(a_subbox);
  }

///
/**
   Returns the irregular cells that have non-zero boundary area
*/

  IntVectSet EBISBox::boundaryIVS(const Box& a_subbox) const
  {
    IntVectSet ivs = getIrregIVS(a_subbox);
    IntVectSet rtn = ivs;

    IVSIterator it(ivs);
    for (;it.ok(); ++it)
    {
      if (m_data.bndryArea(VolIndex(it(),0)) == 0) rtn -= it();
    }
    return rtn;
  }

/*******************************/
  const Box&
  EBISBox::getDomain() const
  {
    return m_graph.getDomain();
  }
/*******************************/

/*******************************/
  Vector<VolIndex>
  EBISBox::getVoFs(const VolIndex& a_vof,
                   const int& a_dir,
                   const Side::LoHiSide& a_sd,
                   const int& a_steps) const
  {
    return m_graph.getVoFs(a_vof, a_dir, a_sd, a_steps);
  }
/*******************************/

/*******************************/
  int
  EBISBox::numFaces(const VolIndex& a_vof,
                    const int& a_idir,
                    const Side::LoHiSide& a_sd) const
  {
    Vector<FaceIndex> faces = getFaces(a_vof, a_idir, a_sd);
    int retval = faces.size();
    return retval;
  }
/*******************************/
  Real
  EBISBox::volFrac(const VolIndex& a_vof) const
  {
    Real retval;
    if (isRegular(a_vof.gridIndex()))
    {
      retval = 1.0;
    }
    else if (isCovered(a_vof.gridIndex()))
    {
      retval = 0.0;
    }
    else
    {
      retval = m_data.volFrac(a_vof);
    }

    return retval;
  }
/*******************************/
  Real
  EBISBox::areaFracScaling(const VolIndex& a_vof) const
  {
    Real alphaMax = 0;
    for (int idir=0; idir<SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        Vector<FaceIndex> faces = getFaces(a_vof, idir, sit());
        for (int iface=0; iface<faces.size(); iface++)
        {
          alphaMax = std::max(areaFrac(faces[iface]), alphaMax);
        }
      }
    }
    if (alphaMax > 0)
    {
      return (1./alphaMax);
    }
    else
    {
      return 1.;
    }
  }
/*******************************/
  Real
  EBISBox::sumArea(const VolIndex& a_vof,
                   const int& a_idir,
                   const Side::LoHiSide& a_sd) const
  {
    Real retval=0;
    if (isRegular(a_vof.gridIndex()))
    {
      retval = 1.0;
    }
    else if (isCovered(a_vof.gridIndex()))
    {
      retval = 0.0;
    }
    else
    {
      retval = 0.0;
      Vector<FaceIndex> faces = getFaces(a_vof, a_idir, a_sd);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        retval += areaFrac(faces[iface]);
      }
    }
    return retval;
  }
/*******************************/
  Vector<VolIndex>
  EBISBox::refine(const VolIndex& a_coarVoF) const
  {
    return(m_graph.refine(a_coarVoF));
  }
/*******************************/
  VolIndex
  EBISBox::coarsen(const VolIndex& a_fineVoF) const
  {
    return(m_graph.coarsen(a_fineVoF));
  }
/*******************************/
  bool
  EBISBox::isConnected(const VolIndex& a_vof1,
                       const VolIndex& a_vof2) const
  {
    return m_graph.isConnected(a_vof1, a_vof2);
  }
/*******************************/
  Real
  EBISBox::areaFrac(const FaceIndex& a_face) const
  {
    Real retval;

    Box region = m_graph.getRegion();
    const IntVect& loiv = a_face.gridIndex(Side::Lo);
    const IntVect& hiiv = a_face.gridIndex(Side::Hi);
    if (region.contains(loiv) && isRegular(loiv) )
    {
      retval = 1.0;
    }
    else if (region.contains(hiiv) && isRegular(hiiv))
    {
      retval = 1.0;
    }
    else if (region.contains(loiv) && isCovered(loiv))
    {
      retval =  0.0;
    }
    else if (region.contains(hiiv) && isCovered(hiiv))
    {
      retval =  0.0;
    }
    else
    {
      retval = m_data.areaFrac(a_face);
    }
    return retval;
  }
/*******************************/
  RealVect
  EBISBox::normal(const VolIndex& a_vof) const
  {
    RealVect retval;
    const IntVect& iv = a_vof.gridIndex();
    if (isRegular(iv)) 
    {
      retval = BASISREALV(0);
    }
    else if (isCovered(iv))
    {
      retval = BASISREALV(0);
    }
    else
    {
      retval = m_data.normal(a_vof);
      Real tol = 1.0e-20;
      //this can happen if the cell is really regular
      bool allZeros = AMREX_D_TERM((std::abs(retval[0]) < tol), && (std::abs(retval[1]) < tol), && (std::abs(retval[2]) < tol));

      if(allZeros)
      {
        retval = BASISREALV(0);
      }
    } //end else (vof is irregular)

    return retval;
  }

  /*******************************/
  RealVect
  EBISBox::centroid(const FaceIndex& a_face) const
  {
    RealVect retval;

    const IntVect& loiv = a_face.gridIndex(Side::Lo);
    const IntVect& hiiv = a_face.gridIndex(Side::Hi);
    Box region = m_graph.getRegion();
    if (region.contains(loiv) && isRegular(loiv))
    {
      retval = RealVect::Zero;
    }
    else if (region.contains(loiv) && isCovered(loiv))
    {
      retval = RealVect::Unit;
    }
    else if (region.contains(hiiv) && isRegular(hiiv))
    {
      retval = RealVect::Zero;
    }
    else if (region.contains(hiiv) && isCovered(hiiv))
    {
      retval = RealVect::Unit;
    }
    else
    {
      retval = m_data.centroid(a_face);
    }
    return retval;

  }
/*******************************/
  RealVect
  EBISBox::centroid(const VolIndex& a_vof) const
  {
    RealVect retval;
    const IntVect& iv = a_vof.gridIndex();
    if (isRegular(iv) || (isCovered(iv)))
    {
      retval = RealVect::Zero;
    }
    else
    {
      retval = m_data.centroid(a_vof);
    }
    return retval;
  }
/*******************************/
  RealVect
  EBISBox::bndryCentroid(const VolIndex& a_vof) const
  {
    RealVect retval;
    const IntVect& iv = a_vof.gridIndex();
    if (isRegular(iv))
    {
      retval = RealVect::Unit;
      retval *= -1.0;
    }
    else if (isCovered(iv))
    {
      retval = RealVect::Unit;
      retval *= -1.0;
    }
    else
    {
      retval = m_data.bndryCentroid(a_vof);
    }
    return retval;
  }

  /*******************************/
  Real
  EBISBox::bndryArea(const VolIndex& a_vof) const
  {
    //  Real retval = PolyGeom::bndryArea(a_vof, *this);

    Real retval;
    const IntVect& iv = a_vof.gridIndex();
    if (isRegular(iv))
    {
      retval = 0.0;
    }
    else if (isCovered(iv))
    {
      retval = -1.0;
    }
    else
    {
      retval = m_data.bndryArea(a_vof);
    }
    return retval;
  }

/*******************************/
  void
  EBISBox::define(const EBGraph&  a_graph,
                  const EBData&   a_data)
  {
    m_graph = a_graph;
    m_data  = a_data;
  }
/*******************************/
  void EBISBox::setDomain(const Box& a_domain)
  {
    m_graph.setDomain(a_domain);
  }
/*******************************/
  void EBISBox::setToAllRegular()
  {
    m_graph.setToAllRegular();
  }
/*******************************/
  void EBISBox::setToAllCovered()
  {
    m_graph.setToAllCovered();
  }
/*******************************/
  Vector<FaceIndex>
  EBISBox::refine(const FaceIndex& a_coarFace,const EBISBox& a_fineEBISBox) const
  {
    return m_graph.refine(a_coarFace, a_fineEBISBox.m_graph);
  }
/*******************************/
  FaceIndex
  EBISBox::coarsen(const FaceIndex& a_fineFace) const
  {
    return m_graph.coarsen(a_fineFace);
  }
/*******************************/
  EBISBox& EBISBox::operator=(const EBISBox& a_ebiin)
  {
    if (&a_ebiin != this)
    {
      m_graph = a_ebiin.m_graph;
      m_data  = a_ebiin.m_data;
    }
    return *this;
  }
/*******************************/
  EBISBox::EBISBox(const EBISBox& a_ebiin)
  {
    m_graph = a_ebiin.m_graph;
    m_data  = a_ebiin.m_data;
  }
/*******************************/
  const EBGraph& EBISBox::
  getEBGraph() const
  {
    return m_graph;
  }
/*******************************/
  const EBData& EBISBox::
  getEBData() const
  {
    return m_data;
  }
/*******************************/
  bool EBISBox::
  operator==(const EBISBox& a_ebiin)
  {
    return (m_graph == a_ebiin.m_graph) &&  (m_data == a_ebiin.m_data);
  }
/*******************************/
}
