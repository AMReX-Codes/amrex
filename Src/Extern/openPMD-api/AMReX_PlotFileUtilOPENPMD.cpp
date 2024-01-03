#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <AMReX_PlotFileUtilOPENPMD.H>

#include <fstream>
#include <iomanip>

namespace amrex::openpmd_api
{
    /* global handler, activate with InitHandler() & deactivate with CloseHandler() */
    std::unique_ptr< AMReX_openPMDHandler > m_OpenPMDHandler = nullptr;

    std::unique_ptr<AMReX_openPMDHandler> InitUserHandler(const std::string& prefix)
    {
      std::string filePath;
      if (prefix.empty())
        {
          ParmParse pp;
          pp.query("openpmd_directory", filePath);
        }
      else {
        filePath = prefix;
      }

      return std::make_unique<AMReX_openPMDHandler>(filePath);
    }

    void CloseUserHandler(std::unique_ptr<AMReX_openPMDHandler>& userHandler)
    {
      if (userHandler == nullptr) { return; }

      userHandler.reset(nullptr);
    }

    void InitHandler(const std::string& prefix)
    {
      std::string filePath;
      if (prefix.empty())
        {
          ParmParse pp;
          pp.query("openpmd_directory", filePath);
        }
      else {
        filePath = prefix;
      }

      if (m_OpenPMDHandler == nullptr) {
        m_OpenPMDHandler = std::make_unique<AMReX_openPMDHandler>(filePath);
      } else if (m_OpenPMDHandler->m_Writer != nullptr) {
        if (m_OpenPMDHandler->m_Writer->m_openPMDPrefix !=  filePath) {
          m_OpenPMDHandler = std::make_unique<AMReX_openPMDHandler>(filePath);
        }
      }
      // already using the directory, no action needed
    }

    void UseCustomWriter(AMReX_openPMDWriter* w)
    {
      BL_ASSERT ( m_OpenPMDHandler != nullptr );
      BL_ASSERT ( w != nullptr );

      // so the openpmd filepath assigned from input file is still in use
      w->m_openPMDPrefix = m_OpenPMDHandler->m_Writer->m_openPMDPrefix;
      m_OpenPMDHandler->m_Writer.reset(w);
    }

    void CloseHandler()
    {
      m_OpenPMDHandler.reset();
    }

    void SetStep(int ts)
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr)) {
        return;
      }

      m_OpenPMDHandler->m_Writer->SetStep(ts);
    }

    void CloseStep(int ts)
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr)) {
        return;
      }

      m_OpenPMDHandler->m_Writer->CloseStep(ts);
    }

    void WriteSingleLevel (//const std::string &plotfilename,
                           const MultiFab &mf,
                           const Vector<std::string> &varnames,
                           const Geometry &geom,
                           Real t,
                           //int level_step,
                           const std::string &versionName,
                           const std::string &levelPrefix,
                           const std::string &mfPrefix,
                           const Vector<std::string>& extra_dirs)
    {
      Vector<const MultiFab*> v_mf(1,&mf);
      Vector<Geometry> v_geom(1,geom);
      //Vector<int> v_level_steps(1,level_step);
      Vector<IntVect> ref_ratio;

      WriteMultiLevel(v_mf, varnames, v_geom, t, /*v_level_steps,*/ ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
    }

    void WriteMultiLevel (
                          //int nlevels, // will write all levels in mf & geom
                          const Vector<const MultiFab*> &mf,
                          const Vector<std::string> &varnames,
                          const Vector<Geometry> &geom,
                          Real time,
                          //const Vector<int> &level_steps,
                          const Vector<IntVect> & /*ref_ratio*/,
                          const std::string & /*versionName*/,
                          const std::string & /*levelPrefix*/,
                          const std::string & /*mfPrefix*/,
                          const Vector<std::string>& /*extra_dirs*/)
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr)) { return; }

      BL_ASSERT ( geom.size() == mf.size() );
      BL_ASSERT ( mf[0]->nComp() <= varnames.size() );

      m_OpenPMDHandler->m_Writer->WriteMesh(varnames,
                                            mf, //amrex::GetVecOfConstPtrs(mf),
                                            geom,
                                            //level_steps[0],
                                            time);
    }
} // namespace amrex::openpmd_api
