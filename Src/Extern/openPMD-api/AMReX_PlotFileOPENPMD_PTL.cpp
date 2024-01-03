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
#include <openPMD/openPMD.hpp>

#include <regex>
#include <fstream>
#include <iomanip>


namespace amrex::openpmd_api {

    bool AMReX_openPMDWriter::AllocatePtlProperties(openPMD::ParticleSpecies& currSpecies,
                                                    const amrex::Vector<int>& write_real_comp,
                                                    const amrex::Vector<std::string>& real_comp_names,
                                                    const amrex::Vector<int>& write_int_comp,
                                                    const amrex::Vector<std::string>& int_comp_names,
                                                    const unsigned long long np) const
    {
      SetupPos(currSpecies, np);

      // Allocate _all_ datasets of dtype.
      // handle scalar and non-scalar records by name
      auto const lf_compRecordInit = [&currSpecies](const amrex::Vector<int>& write_comp,
                                                    const amrex::Vector<std::string>& comp_names,
                                                    openPMD::Dataset& dtype)
      {
        auto const min_counter = std::min(write_comp.size(), comp_names.size());
        for (int i = 0; i < min_counter; ++i)
          {
            if (write_comp[i]) {
              helper::getComponentRecord(currSpecies,  comp_names[i]).resetDataset(dtype);
            }
          }
      };
      auto dtype_real = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np}, m_openPMDDatasetOptions);
      lf_compRecordInit(write_real_comp, real_comp_names, dtype_real);

      auto dtype_int  = openPMD::Dataset(openPMD::determineDatatype<int>(), {np}, m_openPMDDatasetOptions);
      lf_compRecordInit(write_int_comp, int_comp_names, dtype_int);

      return true;
    }

    void AMReX_openPMDWriter::SetupPos(openPMD::ParticleSpecies& currSpecies,
                                       const unsigned long long& np) const
    {
      auto realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np}, m_openPMDDatasetOptions);
      auto idType = openPMD::Dataset(openPMD::determineDatatype< uint64_t >(), {np}, m_openPMDDatasetOptions);

      auto const positionComponents = /*helper::*/getParticlePositionComponentLabels();
      for( auto const& comp : positionComponents )
        {
          currSpecies["position"][comp].resetDataset( realType );
        }

      auto const * const scalar = openPMD::RecordComponent::SCALAR;
      currSpecies["id"][scalar].resetDataset( idType );
    }


    unsigned long long AMReX_openPMDWriter::GetGrandOffset() const
    {
      return 0;
    }

    // from warpx  SetConstParticleRecordsEDPIC ()
    //
    // constant values, just call before flushing
    //
    void AMReX_openPMDWriter::SetupConstant(openPMD::ParticleSpecies& currSpecies,
                                            const unsigned long long& np) const
    {
      auto realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np}, m_openPMDDatasetOptions);

      auto const positionComponents = getParticlePositionComponentLabels();
      for( auto const& comp : positionComponents ) {
        currSpecies["positionOffset"][comp].resetDataset( realType );
      }

      // make constant
      using namespace amrex::literals;
      for( auto const& comp : positionComponents ) {
        currSpecies["positionOffset"][comp].makeConstant( 0._prt );
      }

      currSpecies["positionOffset"].setUnitDimension( helper::getUnitDimension("positionOffset") );
      currSpecies["position"].setUnitDimension( helper::getUnitDimension("position") );

      SetParticleSpecieAttributes(currSpecies);
    }

}
