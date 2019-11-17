#ifndef WARPX_OPEN_PMD_H_
#define WARPX_OPEN_PMD_H_

#include <MultiParticleContainer.H> // has AMReX Vector etc used below
//#include <WarpX.H>

#include <openPMD/openPMD.hpp>

//
// helper class
//
class Timer
{
public:
  Timer(const char* tag) {m_Tag = tag; m_Start = amrex::second();}
  ~Timer() {m_End = amrex::second(); amrex::Print()<<m_Tag<<" took:"<<m_End - m_Start<<" seconds\n";}
private:
  amrex::Real m_Start;
  amrex::Real m_End;  
  std::string m_Tag;
};


//
//
class WarpXParticleCounter
{
public:
  WarpXParticleCounter(const std::unique_ptr<WarpXParticleContainer>& pc);  
  unsigned long GetTotalNumParticles() {return m_Total;}

  std::vector<unsigned long long> m_ParticleOffsetAtRank;;
  std::vector<unsigned long long> m_ParticleSizeAtRank;;
private:
  // get the offset in the overall particle id collection
  void GetParticleOffsetOfProcessor(const long& numParticles,
				    unsigned long long& offset,
				    unsigned long long& sum)  const ;


  int m_MPIRank = 0;
  int m_MPISize = 1;  
  
  unsigned long long m_Total = 0;

  std::vector<unsigned long long> m_ParticleCounterByLevel;
};


//
//
//
class WarpXOpenPMDPlot
{
public:
  WarpXOpenPMDPlot(bool, const std::string& filetype);
  //WarpXOpenPMDPlot(const std::string& dir, const std::string& fileType);
  ~WarpXOpenPMDPlot();

  void SetStep(int ts);
  void SetOpenPMDType(const std::string& type);

  void WriteOpenPMDParticles(const std::unique_ptr<MultiParticleContainer>&);

  void WriteOpenPMDFields(//const std::string& filename,
              const std::vector<std::string>& varnames,
              const amrex::MultiFab& mf,
              const amrex::Geometry& geom,
              const int iteration, const double time ) const;


private:
  void Init(const std::string& filename,
	    openPMD::AccessType accessType);
  
  /*
  // get the offset in the overall particle id collection
  void GetParticleOffsetOfProcessor(const long& numParticles,
				    unsigned long long& offset,
				    unsigned long long& sum)  const ;
  */
  void SetupPos(openPMD::ParticleSpecies& currSpecies,
        const unsigned long long& np) const ;

  void SetupRealProperties(openPMD::ParticleSpecies& currSpecies,
			   const amrex::Vector<int>& write_real_comp,
			   const amrex::Vector<std::string>& real_comp_names,
			   unsigned long long np) const;

  void SaveRealProperty(WarpXParIter& pti, //int, int,
            openPMD::ParticleSpecies& currSpecies,
            unsigned long long offset,
            unsigned long long numParticles,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<std::string>& real_comp_names) const;


  void SavePlotFile(const std::unique_ptr<WarpXParticleContainer>& pc,
            const std::string& name,
            int iteration,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<int>& write_int_comp,
            const amrex::Vector<std::string>& real_comp_names,
            const amrex::Vector<std::string>&  int_comp_names) const;


  void GetFileName(std::string& filename, int ts);
  //std::string m_Dir;
  openPMD::Series* m_Series;

  int m_MPIRank = 0;
  int m_MPISize = 1;

  int m_NumSoARealAttributes = PIdx::nattribs;
  int m_NumAoSRealAttributes = 0; // WarpX definition

  bool m_OneFilePerTS = true;
  std::string m_OpenPMDFileType = "bp";
  int m_CurrentStep  = -1;
};


#endif // WARPX_OPEN_PMD_H
