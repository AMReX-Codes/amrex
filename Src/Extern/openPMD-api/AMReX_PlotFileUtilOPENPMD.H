#ifndef AMREX_PLOTFILE_OPENPMD_API_H
#define AMREX_PLOTFILE_OPENPMD_API_H

#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileDataImpl.H>

#ifdef AMREX_PARTICLES
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#endif

#ifdef AMREX_USE_OPENPMD_API
#include <openPMD/openPMD.hpp>
#endif

#include <string>
#include <memory>

namespace amrex::openpmd_api
{
    namespace helper
    {
      inline void
      setOpenPMDUnit ( openPMD::Mesh mesh, std::string const& field_name )
      {
        if (field_name[0] == 'E'){
          mesh.setUnitDimension({
              {openPMD::UnitDimension::L,  1},
              {openPMD::UnitDimension::M,  1},
              {openPMD::UnitDimension::T, -3},
              {openPMD::UnitDimension::I, -1},
            });
        } else if (field_name[0] == 'B'){
          mesh.setUnitDimension({
              {openPMD::UnitDimension::M,  1},
              {openPMD::UnitDimension::I, -1},
              {openPMD::UnitDimension::T, -2}
              });
        }else if (field_name[0] == 'j'){ // current
          mesh.setUnitDimension({
              {openPMD::UnitDimension::L, -2},
              {openPMD::UnitDimension::I,  1},
            });
        } else if (field_name.substr(0,3) == "rho"){ // charge density
          mesh.setUnitDimension({
              {openPMD::UnitDimension::L, -3},
              {openPMD::UnitDimension::I,  1},
              {openPMD::UnitDimension::T,  1},
            });
        }
      } // name space helper

      // from warpx Utils/RelativeCellPosition.cpp
      inline std::vector< double >
      getRelativeCellPosition(amrex::MultiFab const& mf)
      {
        amrex::IndexType const idx_type = mf.ixType();

        std::vector< double > relative_position(AMREX_SPACEDIM, 0.0);
        // amrex::CellIndex::CELL means: 0.5 from lower corner for that index/direction
        // amrex::CellIndex::NODE means: at corner for that index/direction
        // WarpX::do_nodal means: all indices/directions on CellIndex::NODE
        for (int d = 0; d < AMREX_SPACEDIM; d++)
          {
            if (idx_type.cellCentered(d)) {
              relative_position.at(d) = 0.5;
            }
          }
        return relative_position;
      }

      inline std::vector<std::uint64_t>
      getReversedVec( const IntVect& v )
      {
        // Convert the IntVect v to and std::vector u
        std::vector<std::uint64_t> u = {
          AMREX_D_DECL(
                       static_cast<std::uint64_t>(v[0]),
                       static_cast<std::uint64_t>(v[1]),
                       static_cast<std::uint64_t>(v[2])
                       )
        };
        // Reverse the order of elements, if v corresponds to the indices of a
        // Fortran-order array (like an AMReX FArrayBox)
        // but u is intended to be used with a C-order API (like openPMD)
        std::reverse( u.begin(), u.end() );

        return u;
      }

      inline std::vector<double>
      getReversedVec( const Real* v )
      {
        // Convert Real* v to and std::vector u
        std::vector<double> u = {
          AMREX_D_DECL(
                       static_cast<double>(v[0]),
                       static_cast<double>(v[1]),
                       static_cast<double>(v[2])
                       )
        };
        // Reverse the order of elements, if v corresponds to the indices of a
        // Fortran-order array (like an AMReX FArrayBox)
        // but u is intended to be used with a C-order API (like openPMD)

        std::reverse( u.begin(), u.end() );

        return u;
      }
    }

    struct AMReX_VarNameParser; // forward declaration
    class  AMReX_openPMDWriter;
    class  AMReX_openPMDHandler;

    /*
    */

    void UseCustomWriter(AMReX_openPMDWriter* w);

    std::unique_ptr<AMReX_openPMDHandler> InitUserHandler(const std::string& prefix);
    void CloseUserHandler(std::unique_ptr<AMReX_openPMDHandler>& userHandler);

    void InitHandler(const std::string& prefix);
    //void SetCustomizer(AMReX_Optional* f);

    void SetStep(int ts);
    void CloseStep(int ts);
    void CloseHandler();

    void WriteSingleLevel (//const std::string &plotfilename,
                           const MultiFab &mf,
                           const Vector<std::string> &varnames,
                           const Geometry &geom,
                           Real t,
                           //int level_step,
                           const std::string &ignored_versionName = "HyperCLaw-V1.1",
                           const std::string &ignored_levelPrefix = "Level_",
                           const std::string &ignored_mfPrefix = "Cell",
                           const Vector<std::string>& ignored_extra_dirs = Vector<std::string>());

     void WriteMultiLevel (
                           //int nlevels, will all levels in mf & geom
                           const Vector<const MultiFab*> &mf,
                           const Vector<std::string> &varnames,
                           const Vector<Geometry> &geom,
                           Real time,
                           //const Vector<int> &level_steps,
                           const Vector<IntVect> &ref_ratio,
                           const std::string &ignored_versionName = "HyperCLaw-V1.1",
                           const std::string &ignored_levelPrefix = "Level_",
                           const std::string &ignored_mfPrefix = "Cell",
                           const Vector<std::string>& ignored_extra_dirs = Vector<std::string>());


    ////////////////////////////
    struct AMReX_VarNameParser
    {
      AMReX_VarNameParser(std::string const& varname);
      virtual ~AMReX_VarNameParser () = default;
      AMReX_VarNameParser (AMReX_VarNameParser const&) = default;
      AMReX_VarNameParser (AMReX_VarNameParser &&) = delete;
      AMReX_VarNameParser& operator= (AMReX_VarNameParser const&) = default;
      AMReX_VarNameParser& operator= (AMReX_VarNameParser &&) = delete;

      void GetFieldNameModeInt(const std::string& varname);
      void GetMeshCompNames (int meshLevel);

      [[nodiscard]] std::vector< std::string > getFieldComponentLabels () const
      {
        using vs = std::vector< std::string >;
        if (m_ThetaMode)
          {
            // if we write individual modes
            vs fieldComponents{"r", "t", "z"};
            return fieldComponents;
          }
        else
          {
            // if we just write reconstructed fields at theta=0 or are Cartesian
            // note: 1D3V and 2D3V simulations still have 3 components for the fields
            vs fieldComponents{"x", "y", "z"};
            return fieldComponents;
          }
      }

      [[nodiscard]] virtual std::vector< std::string >
      getFieldAxisLabels () const
      {
        using vs = std::vector< std::string >;

        // temporary resolution
#if AMREX_SPACEDIM == 1
        vs const axisLabels{"z"};
#elif AMREX_SPACEDIM == 2
        vs const axisLabels{"x", "z"};
#elif AMREX_SPACEDIM == 3
        vs const axisLabels{"x", "y", "z"};  // x varies fastest in memory
#else
#  error Unable to label more than 3d
        // no labels to be addressed
        vs const axisLabels{}
#endif
        // revert to C order (fastest varying index last)
        return {axisLabels.rbegin(), axisLabels.rend()};
      } // getFieldAxisLabels


      /////  members
      std::string m_FieldName;
      std::string m_CompName;
      bool m_ThetaMode = false;
      int m_ModeIndex=-1;

    };


    ////////////////////////////
    ////////////////////////////
    class AMReX_openPMDWriter {
    public:
      AMReX_openPMDWriter(std::string prefix,
                          openPMD::IterationEncoding ie,
                          std::string filetype,
                          std::string openpmdOptions);

      AMReX_openPMDWriter() = default;
      AMReX_openPMDWriter (AMReX_openPMDWriter const&) = delete;
      AMReX_openPMDWriter (AMReX_openPMDWriter &&) = delete;
      AMReX_openPMDWriter& operator= (AMReX_openPMDWriter const&) = delete;
      AMReX_openPMDWriter& operator= (AMReX_openPMDWriter &&) = delete;

      virtual ~AMReX_openPMDWriter();

      [[nodiscard]] virtual openPMD::Iteration GetIteration (int const iteration) const
      {
        return m_Series->writeIterations()[iteration];
      }

      virtual void SetStep(int ts);
      virtual void CloseStep(int ts);
      virtual void Init(openPMD::Access access);

      virtual void WriteMesh(const std::vector<std::string>& varnames,
                             const amrex::Vector<const amrex::MultiFab*>& mf,
                             const amrex::Vector<amrex::Geometry>& geom,
                             //const Vector<int> &iteration,
                             //const int iteration, /* note:  all levels are outputting the same step */
                             double time ) const;

      template<class PC>
      void DumpParticles(PC const& pc,
                         const std::string& name,
                         //const int iteration,
                         const amrex::Vector<int>& write_real_comp, // size = NStructReal + NArrayReal
                         const amrex::Vector<int>& write_int_comp,  // size = NStructInt  + NArrayInt
                         const amrex::Vector<std::string>& real_comp_names,
                         const amrex::Vector<std::string>&  int_comp_names
                         ) const;

      template<class PC, typename FUNC_pc, typename FUNC_pit>
      void DumpParticles(PC const& pc,
                         const std::string& name,
                         //const int iteration,
                         const amrex::Vector<int>& write_real_comp, // size = NStructReal + NArrayReal
                         const amrex::Vector<int>& write_int_comp,  // size = NStructInt  + NArrayInt
                         const amrex::Vector<std::string>& real_comp_names,
                         const amrex::Vector<std::string>&  int_comp_names,
                         FUNC_pc&& pc_func,
                         FUNC_pit&& posId_func
                         ) const;

      [[nodiscard]] virtual bool AllocatePtlProperties(openPMD::ParticleSpecies& currSpecies,
                                         const amrex::Vector<int>& write_real_comp,
                                         const amrex::Vector<std::string>& real_comp_names,
                                         const amrex::Vector<int>& write_int_comp,
                                         const amrex::Vector<std::string>& int_comp_names,
                                         unsigned long long np) const;

      virtual void SetParticleSpecieAttributes(openPMD::ParticleSpecies& /*currSpecies*/) const
      {}

      void SetupPos(openPMD::ParticleSpecies& species, const unsigned long long& np) const;
      [[nodiscard]] virtual unsigned long long GetGrandOffset() const; // to adjust offset/total with data already flushed,  example use case: BTD

      // called once because it is constant valued record
      virtual void SetupConstant(openPMD::ParticleSpecies& currSpecies,
                                 const unsigned long long& np) const;

      [[nodiscard]] virtual std::vector< std::string >
      getParticlePositionComponentLabels() const
      {
#if (AMREX_SPACEDIM == 3)
        std::vector< std::string > positionComponents{"x", "y", "z"};
#elif (AMREX_SPACEDIM == 2)
        std::vector< std::string > positionComponents{"x", "y"};
#elif (AMREX_SPACEDIM == 1)
        std::vector< std::string > positionComponents{"x"};
#else
#   error Unknown dimensionality.
#endif
        return positionComponents;
      }

      template<class PC>
      void
      SetupRealProperties(PC const& pc,
                          openPMD::ParticleSpecies& currSpecies,
                          const amrex::Vector<int>& write_real_comp,
                          const amrex::Vector<std::string>& real_comp_names,
                          const amrex::Vector<int>& write_int_comp,
                          const amrex::Vector<std::string>& int_comp_names,
                          unsigned long long np) const;
      /*
      */
      template<class PIt>
      void
      SaveRealProperty (PIt& pti,
                        openPMD::ParticleSpecies& currSpecies,
                        unsigned long long offset,
                        const amrex::Vector<int>& write_real_comp,
                        const amrex::Vector<std::string>& real_comp_names,
                        const amrex::Vector<int>& write_int_comp,
                        const amrex::Vector<std::string>& int_comp_names) const;


      std::unique_ptr<openPMD::Series> m_Series = nullptr;
      std::string m_openPMDPrefix = std::string();

      int m_openPMDMinDigits = 6;
      std::string m_openPMDDatasetOptions = "{}";

      openPMD::IterationEncoding  m_openPMDEncoding = openPMD::IterationEncoding::fileBased;

      std::string m_openPMDFileType = "bp";
      std::string m_openPMDSeriesOptions = std::string();

    protected:
      void CompSetup(int lev,
                     openPMD::Container< openPMD::Mesh >& meshes,
                     amrex::Geometry& full_geom,
                     const std::vector<std::string>& varnames,
                     const amrex::MultiFab* mf) const;

      static void CompStorage(int lev,
                       openPMD::Container< openPMD::Mesh >& meshes,
                       amrex::Geometry& full_geom,
                       const std::vector<std::string>& varnames,
                       const amrex::MultiFab* mf);


      virtual void SetupFields (openPMD::Container< openPMD::Mesh >& meshes,
                        amrex::Geometry& full_geom) const;

      virtual void SetupMeshComp (openPMD::Mesh& mesh,
                          const amrex::Geometry& full_geom,
                          amrex::MultiFab const& mf,
                          const AMReX_VarNameParser& varName) const;

      void GetFileName(std::string& filepath) const;

      template<class PIt>
      void StoreAoS_Real(PIt& pti,
                         openPMD::ParticleSpecies& currSpecies,
                         const amrex::Vector<int>& write_real_comp,
                         const amrex::Vector<std::string>& real_comp_names,
                         unsigned long long offset) const;

      template<class PIt>
      void
      SavePosId(PIt& pti,
                openPMD::ParticleSpecies& currSpecies,
                unsigned long long offset) const;

      template<class PIt>
      void
      SaveSpecieId(PIt& pti,
                  openPMD::ParticleSpecies& currSpecies,
                  unsigned long long offset) const;

      template<class PIt>
      void
      StoreAoS_Int(PIt& pti,
                   openPMD::ParticleSpecies& currSpecies,
                   const amrex::Vector<int>& write_int_comp,
                   const amrex::Vector<std::string>& int_comp_names,
                   unsigned long long offset) const;


      template<class PIt>
      void StoreSoAReal(PIt& pti,
                        openPMD::ParticleSpecies& currSpecies,
                        const amrex::Vector<int>& write_real_comp,
                        const amrex::Vector<std::string>& real_comp_names,
                        unsigned long long offset)  const;


      template<class PIt>
      void StoreSoAInt(PIt& pti,
                       openPMD::ParticleSpecies& currSpecies,
                       const amrex::Vector<int>& write_int_comp,
                       const amrex::Vector<std::string>& int_comp_names,
                       unsigned long long offset) const;

      int m_CurrentStep = -1;

    };



    ////////////////////////////
    class AMReX_openPMDHandler {
    public:
      AMReX_openPMDHandler(std::string const& prefix = std::string());

      std::unique_ptr<AMReX_openPMDWriter> m_Writer;

    private:
      void CreateWriter(const std::string& prefix = std::string());

    }; // class AMReX_openPMDHandler


} // namespace amrex::openpmd_api


#include <AMReX_PlotFileUtilOPENPMD_PTLImpl.H>

#endif // AMREX_PLOTFILE_OPENPMD_API_H
