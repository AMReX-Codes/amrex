
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

#include "AMReX_STLAsciiReader.H"
#include "AMReX_STLExplorer.H"
#include "AMReX_STLMesh.H"
#include "AMReX_STLBox.H"
#include "AMReX_CellEdge.H"
#include "AMReX_STLIF.H"


namespace amrex
{
  STLIF::STLIF(const string& a_filename)

  {
    BL_PROFILE("STLIF::STLIF_file");

    m_filename = a_filename;

    makeExplorer();
  }

  STLIF::STLIF(const STLIF& a_inputIF)
  {
    BL_PROFILE("STLIF::STLIF_copy");

    m_filename = a_inputIF.m_filename;

    makeExplorer();
  }

  STLIF::~STLIF()
  {
  }

  Real STLIF::value(const RealVect& a_point) const
  {
    Real retval = 0.0;

    amrex::Error("STLIF::value should never be called");

    return retval;
  }

  BaseIF* STLIF::newImplicitFunction() const
  {
    BL_PROFILE("STLIF::newImplicitFunction");

    STLIF* dataFilePtr = new STLIF(m_filename);

    return static_cast<BaseIF*>(dataFilePtr);
  }

  std::shared_ptr<STLExplorer> STLIF::getExplorer() const
  {
    if (!m_explorer)
    {
      amrex::Error("STLIF::getExplorer - STLExplorer not defined yet");
    }

    return m_explorer;
  }

  void STLIF::makeExplorer()
  {
    BL_PROFILE("STLIF::makeExplorer");

    shared_ptr<STLMesh> mesh;

    STLAsciiReader reader(m_filename);
    mesh = reader.GetMesh();


    m_explorer = shared_ptr(new STLExplorer(mesh));
  }
}

