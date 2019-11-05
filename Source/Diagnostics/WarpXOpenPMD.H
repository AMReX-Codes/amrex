#ifndef WARPX_OPEN_PMD_H_
#define WARPX_OPEN_PMD_H_

#include <MultiParticleContainer.H> // has AMReX Vector etc used below
#include <WarpX.H>

#include <openPMD/openPMD.hpp>

class WarpXOpenPMDParticle
{
public: 
  WarpXOpenPMDParticle(const std::string& dir);
  ~WarpXOpenPMDParticle();

  void SaveContainerPlots(const std::unique_ptr<MultiParticleContainer>&);
  void SavePlotFile(const std::unique_ptr<WarpXParticleContainer>& pc,
		    const std::string& name,
		    int iteration,
		    const amrex::Vector<int>& write_real_comp,
		    const amrex::Vector<int>& write_int_comp,
		    const amrex::Vector<std::string>& real_comp_names,
		    const amrex::Vector<std::string>&  int_comp_names) const;
private:
  // get the offset in the overall particle id collection
  void GetParticleOffsetOfProcessor(const long& numParticles, 
				    unsigned long& offset, 
				    unsigned long long& sum)  const ;

  void SetupPos(openPMD::ParticleSpecies& currSpecies, 
		const unsigned long long& np) const ;

  void SetupRealProperties(openPMD::ParticleSpecies& currSpecies, 
			   const amrex::Vector<int>& write_real_comp,
			   const amrex::Vector<std::string>& real_comp_names,
			   unsigned long long np) const;

  void SaveRealProperty(WarpXParIter& pti, //int, int,
			openPMD::ParticleSpecies& currSpecies, 
			unsigned long offset,
			unsigned long long numParticles,
			const amrex::Vector<int>& write_real_comp,
			const amrex::Vector<std::string>& real_comp_names) const;				      


  //std::string m_Dir;
  openPMD::Series* m_Series;

  int m_MPIRank = 0;
  int m_MPISize = 1;

  int m_NumSoARealAttributes = PIdx::nattribs; 
  int m_NumAoSRealAttributes = 0; // WarpX definition

};


#endif // WARPX_OPEN_PMD_H
