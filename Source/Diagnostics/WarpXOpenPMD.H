#ifndef WARPX_OPEN_PMD_H_
#define WARPX_OPEN_PMD_H_

#include <MultiParticleContainer.H> // has AMReX Vector etc used below

#include <openPMD/openPMD.hpp>

//
// helper class
//
class Timer
{
public:
  Timer(const char* tag) {m_Tag = tag; m_Start = amrex::second();}
  ~Timer() {
      m_End = amrex::second();
      amrex::ParallelDescriptor::ReduceRealMax(m_End, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print()<<m_Tag<<" took:"<<m_End - m_Start<<" seconds\n";
  }
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
  /** get the offset in the overall particle id collection
  *
  * @param[out] numParticles particles on this processor  / amrex fab
  * @param[out] offset particle offset over all, mpi-global amrex fabs
  * @param[out] sum number of all particles from all amrex fabs
  */
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
/** Writer logic for openPMD particles and fields */
class WarpXOpenPMDPlot
{
public:
  // not  using const string, to  allow std::move to be effective
  WarpXOpenPMDPlot(bool, std::string& filetype);

  //WarpXOpenPMDPlot(const std::string& dir, const std::string& fileType);
  ~WarpXOpenPMDPlot();

  void SetStep(int ts);

  void WriteOpenPMDParticles(const std::unique_ptr<MultiParticleContainer>&);

  void WriteOpenPMDFields(
              const std::vector<std::string>& varnames,
              const amrex::MultiFab& mf,
              const amrex::Geometry& geom,
              const int iteration, const double time ) const;


private:
  void Init(//const std::string& filename,
        openPMD::AccessType accessType);

  /** This function sets up the entries for storing the particle positions in an openPMD  species
  *
  * @param[in]  currSpecies The openPMD species
  * @param[in]  np  Number of particles
  */
  void SetupPos(openPMD::ParticleSpecies& currSpecies,
        const unsigned long long& np) const ;

  /** This function sets up the entries for particle properties
   *
   * @param[in] currSpecies The openPMD species
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   * @param[in] np  Number of particles
   */
  void SetupRealProperties(openPMD::ParticleSpecies& currSpecies,
               const amrex::Vector<int>& write_real_comp,
               const amrex::Vector<std::string>& real_comp_names,
               unsigned long long np) const;

  /** This function saves the values of the entries for particle properties
   *
   * @param[in] pti WarpX particle iterator
   * @param[in] currSpecies The openPMD species to save to
   * @param[in] offset offset to start saving  the particle iterator contents
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   */
  void SaveRealProperty(WarpXParIter& pti, //int, int,
            openPMD::ParticleSpecies& currSpecies,
            unsigned long long offset,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<std::string>& real_comp_names) const;

  /** This function saves the plot file
   *
   * @param[in] pc WarpX particle container
   * @param[in] iteration timestep
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   * @param[in] write_int_comp The int attribute ids, from WarpX
   * @param[in] int_comp_names The int attribute names, from WarpX
   */
  void SavePlotFile(const std::unique_ptr<WarpXParticleContainer>& pc,
            const std::string& name,
            int iteration,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<int>& write_int_comp,
            const amrex::Vector<std::string>& real_comp_names,
            const amrex::Vector<std::string>&  int_comp_names) const;

  // no need for ts in the name, openPMD  handles it
  void GetFileName(std::string& filename);

  //std::string m_Dir;
  std::unique_ptr<openPMD::Series> m_Series;

  int m_MPIRank = 0;
  int m_MPISize = 1;

  int m_NumSoARealAttributes = PIdx::nattribs; //! WarpX' addition particle attributes in SoA
  int m_NumAoSRealAttributes = 0; //! WarpX definition: no additional attributes in particle AoS

  bool m_OneFilePerTS = true;  //! write in openPMD fileBased manner for individual time steps
  std::string m_OpenPMDFileType = "bp"; //! MPI-parallel openPMD backend: bp or h5
  int m_CurrentStep  = -1;
};


#endif // WARPX_OPEN_PMD_H
