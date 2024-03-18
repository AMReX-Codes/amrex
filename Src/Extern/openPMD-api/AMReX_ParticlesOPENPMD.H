#ifndef AMREX_PTL_OPENPMD_API_H
#define AMREX_PTL_OPENPMD_API_H

/** @file This file is included from AMReX_ParticleCounter.H
*   and belongs to the class ParticleContainer_impl.
*   The namespace amrex was already declared.
*
*   Therefore all the definitions in this file belongs to namespace amrex
*/

struct AMReX_PtlCounter
{
  int m_MPIRank = 0;
  int m_MPISize = 1;

  unsigned long long m_Total = 0;

  std::vector<unsigned long long> m_ParticleCounterByLevel;

  [[nodiscard]] unsigned long long GetTotalNumParticles () const { return m_Total;}

  std::vector<unsigned long long> m_ParticleOffsetAtRank;
  std::vector<unsigned long long> m_ParticleSizeAtRank;
};


void CountParticles()
{

  m_PtlCounter.m_MPISize = amrex::ParallelDescriptor::NProcs();
  m_PtlCounter.m_MPIRank = amrex::ParallelDescriptor::MyProc();

  m_PtlCounter.m_ParticleCounterByLevel.resize(this->finestLevel()+1);
  m_PtlCounter.m_ParticleOffsetAtRank.resize(this->finestLevel()+1);
  m_PtlCounter.m_ParticleSizeAtRank.resize(this->finestLevel()+1);

  auto lf_GetParticleOffsetOfProcessor = [&](const long long& numParticles,
                                             unsigned long long& offset,
                                             unsigned long long& sum)  -> void
  {
    std::vector<long long> result(m_PtlCounter.m_MPISize,  0);
    amrex::ParallelGather::Gather (numParticles, result.data(), -1, amrex::ParallelDescriptor::Communicator());

    sum = 0;
    offset = 0;
    for (int i=0;  i<int(result.size());  i++) {
      sum +=  result[i];
      if ( i < m_PtlCounter.m_MPIRank) {
        offset +=  result[i];
      }
    }
  };

  for (auto currentLevel = 0; currentLevel <= this->finestLevel(); currentLevel++)
    {
      // numParticles in this processor
      auto numParticles = static_cast<long long>(this->NumberOfParticlesAtLevel(currentLevel, true, true));

      unsigned long long offset=0; // offset of this level
      unsigned long long sum=0; // numParticles in this level (sum from all processors)
      lf_GetParticleOffsetOfProcessor(numParticles, offset,  sum);

      m_PtlCounter.m_ParticleCounterByLevel[currentLevel] = sum;
      m_PtlCounter.m_ParticleOffsetAtRank[currentLevel] = offset;
      m_PtlCounter.m_ParticleSizeAtRank[currentLevel] = numParticles;

      // adjust offset, it should be numbered after particles from previous levels
      for (auto lv=0; lv<currentLevel; lv++)
        {
          m_PtlCounter.m_ParticleOffsetAtRank[currentLevel] += m_PtlCounter.m_ParticleCounterByLevel[lv];
        }

      m_PtlCounter.m_Total += sum;
    }
}

AMReX_PtlCounter m_PtlCounter;
#endif // AMREX_PTL_OPENPMD_API_H
