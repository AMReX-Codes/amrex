#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>

#include <ParmParse.H>

#include <REAL.H>
#include <IntVect.H>
#include <Array.H>
#include <Amr.H>
#include <AmrLevel.H>
#include <Utility.H>
#include <Geometry.H>
#include <VisMF.H>
#include <Particles_F.H>
#include <RealBox.H>

#ifdef _OPENMP
#include <omp.h>
#endif

struct ParticleBase
{
    int     m_id;
    int     m_cpu;
    int     m_lev;
    int     m_grid;
    IntVect m_cell;
    Real    m_pos[BL_SPACEDIM];

    ParticleBase ()
        :
        m_id(-1),
        m_cpu(-1),
        m_lev(-1),
        m_grid(-1)
        {}

    static IntVect Index (const ParticleBase& p, int lev, const Amr* amr);

    static bool Where (ParticleBase& prt, const Amr* amr, bool update = false);

    static void Reset (ParticleBase& prt, const Amr* amr, bool update);

    static void PeriodicShift (ParticleBase& prt, const Amr* amr);

    static Real InterpDoit (const FArrayBox& fab, const IntVect& hiidx, const Real* frac, int comp);

    static void Interp (const ParticleBase& prt, const Amr* amr, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Amr* amr, const ParticleBase& p, Real* grav);

    static int MaxReaders ();
    //
    // Returns the next particle ID for this processor.
    // Particle IDs start at 1 and are never reused.
    // The pair, consisting of the ID and the CPU on which the particle is "born",
    // is a globally unique identifier for a particle.  The maximum of this value
    // across all processors must be checkpointed and then restored on restart
    // so that we don't reuse particle IDs.
    //
    static int NextID ();
    //
    // Reset on restart.
    //
    static void NextID (int nextid);
    //
    // Used by AssignDensity.
    //
    static bool CrseToFine (const BoxArray& cfba, const IntVect* cells, bool* which);

    static bool FineToCrse (const ParticleBase& p, const Real* plo, const Real* dx, const Box& vbx, const BoxArray& fba, const BoxArray& cba, bool* which, int* cgrid);

    static void AssignDensityCoeffs (const ParticleBase& p, const Real* plo, const Real* dx, Real* fracs,  IntVect* cells);
    //
    // Used for partitioning particles into valid & invalid groups.
    //
    class Invalid
    {
    public:
        bool operator () (const ParticleBase& p) const { return p.m_id <= 0; }
    };
};

std::ostream& operator<< (std::ostream& os, const ParticleBase& p);


template <int N>
struct Particle
    :
    public ParticleBase
{
    //
    // The amount of Real data we hold.
    //
    // In some cases this is:
    //
    // 0 - particle mass
    // 1 - x-velocity
    // 2 - y-velocity
    // 3 - z-velocity
    //
    Real m_data[N];
};

template <int N>
class ParticleContainer
{
public:
    //
    // The type of Particles we hold.
    //
    typedef Particle<N> ParticleType;

    ParticleContainer (Amr* amr)
        :
        m_verbose(1), m_amr(amr) { BL_ASSERT(amr != 0); }

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);

    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false);
    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass);
    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass, const Real shift[]);
    void InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass);

    Real sumParticleMass (const MultiFab& mf, int level) const;
    void sumParticleMomentum (const MultiFab& mf, int lev, Real* mom) const;

    void Increment (MultiFab& mf, int level);

    long IncrementWithTotal (MultiFab& mf, int level);

    void AssignDensity (MultiFab& mf, int level) const;

    void AssignDensity (PArray<MultiFab>& mf) const;

    void MultiplyParticleMass (int lev, Real mult);
 
    Real estTimestep (const MultiFab& grav_vector, Real a, int level, Real cfl) const;

    void Redistribute (bool where_already_called = false);

    bool OK () const;

    void ByteSpread () const;

    void MoveRandom ();

    void MoveRandom (int level);

    // **************************************************************************************************************** 
    //
    // If the particles move only with self-gravity from themselves and the gas, then 
    // we can move them according to the schemes below.
    // The gravitational force must be computed between the calls of the parts of the integration scheme.
    //
    // The following two functions form a PREDICTOR CORRECTOR scheme for integrating the motion of the particles
    // BE CAREFUL: This one uses a NGP interpolation, which is not consistent with the density assignment scheme!
    //
    void movePredict (const MultiFab& grav_vector, int level, Real timestep);
    void moveCorrect (const MultiFab& grav_vector_old, const MultiFab& grav_vector, int level, Real timestep);
    //
    // TODO: the methods should return a constraint on the timestep...
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   standard coordinates -- these rely on CELL-BASED gravity component
    //
    void moveKickDrift (const MultiFab& grav_vector, int level, Real timestep);
    void moveKick      (const MultiFab& grav_vector, int level, Real timestep);
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   standard coordinates -- these rely on EDGE-BASED gravity component
    //
    void moveKickDrift (PArray<MultiFab>& grav_vector, int level, Real timestep);
    void moveKick      (PArray<MultiFab>& grav_vector, int level, Real timestep);
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   comoving coordinates -- these rely on CELL-BASED gravity component
    //
    void moveKickDrift (const MultiFab& grav_vector, int level, Real timestep, Real a_old, Real a_half);
    void moveKick      (const MultiFab& grav_vector, int level, Real timestep, Real a_new, Real a_half);
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   comoving coordinates -- these rely on EDGE-BASED gravity component
    //
    void moveKickDrift (PArray<MultiFab>& grav_vector, int level, Real timestep, Real a_old, Real a_half);
    void moveKick      (PArray<MultiFab>& grav_vector, int level, Real timestep, Real a_new, Real a_half);
    //
    // after the moveKickDrift step the positions of the particles are advanced for a full timestep,
    // so this scheme should work in the overall algorithm...
    //
    // **************************************************************************************************************** 

    void AdvectWithUmac (const MultiFab* umac, int level, Real dt);

    void RemoveParticlesNotAtFinestLevel ();

    void Checkpoint (const std::string& dir, const std::string& name) const;

    void Restart (const std::string& dir, const std::string& file);

    void Timestamp (const std::string& file, const MultiFab& mf, int lev, Real time, const std::vector<int>& idx);

    void WriteAsciiFile (const std::string& file);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

protected:
    //
    // Helper function for Checkpoint().
    //
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Array<int>&    which,
                         Array<int>&    count,
                         Array<long>&   where) const;
    //
    // Helper functions for Restart().
    //
    void Restart_OneDotZero (const std::string& fullname,
                             std::ifstream&     HdrFile);

    void ReadParticles_OneDotZero (int            cnt,
                                   int            grd,
                                   int            lev,
                                   std::ifstream& ifs);
    //
    // We want to store the particles on a level by level and grid by grid basis.  This will
    // make accessing them and doing operations on them more memory efficient since most of our
    // operations on particles are done on a level by level basis or grid by grid basis.
    //
    typedef typename std::deque<ParticleType> PBox;
    //
    // A level of particles is stored in a map indexed by the grid number.
    //
    typedef typename std::map<int,PBox> PMap;
    //
    // The data.
    //
    int         m_verbose;
    Amr*        m_amr;
    Array<PMap> m_particles;
};

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <int N>
void
ParticleContainer<N>::ByteSpread () const
{
    long cnt = 0;

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator it = pmap.begin(), End = pmap.end(); it != End; ++it)
        {
            cnt += it->second.size();
        }
    }

    long mn = cnt*sizeof(ParticleType), mx = mn;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    ParallelDescriptor::ReduceLongMin(mn, IOProc);
    ParallelDescriptor::ReduceLongMax(mx, IOProc);
    ParallelDescriptor::ReduceLongSum(cnt,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "ParticleContainer<N> byte spread across MPI nodes: ["
                  << mn
                  << " ... "
                  << mx
                  << "] total particles: " << cnt << '\n';
    }
}

template <int N>
void
ParticleContainer<N>::InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep)
{
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= N);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const Real strttime = ParallelDescriptor::second();

    // Number of processes that read from the file
    int NReaders = ParticleBase::MaxReaders();

    // Number of chunks we break the redistribution loop into
    // This default is just for the case of NProcs = NReaders = 1.
    int NRedist = 1;

    if (NProcs <= 1024)
    {
       if (NReaders > 1)
          NRedist = 2;
    }
    else if (NProcs <= 4096)
    {
       NReaders = std::max(NReaders,128);
       NRedist = 4;
    }
    else if (NProcs <= 8192)
    {
       NReaders = std::max(NReaders,384);
       NRedist = 32;
    }
    else if (NProcs <= 16384)
    {
       NReaders = std::max(NReaders,512);
       NRedist = 48;
    }

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_amr->finestLevel()+1);

    IntVect lNrep(D_DECL(1,1,1));

    if (Nrep != 0)
        lNrep = *Nrep;

    long  howmany      = 0;
    long  howmany_read = 0;
    //
    // NReaders will read particles into a PBox.
    // Later they'll each put their chunk into m_particles and call Redistribute().
    //
    PBox nparticles;

    if (MyProc < NReaders)
    {
        //
        // Only the first NReaders MPI processes read particles.
        //
        std::ifstream ifs;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
 
        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(file.c_str(), std::ios::in);

        if (!ifs.good())
            BoxLib::FileOpenFailed(file);

        int cnt = 0;

        ifs >> cnt;

        ParticleType p, p_rep;

        const int Chunk = cnt / NReaders;
        //
        // Skip to my section of particles.
        //
        std::string line;

        ifs >> std::ws;  // Eat newline after "cnt" so we skip proper # of particles.

        for (int i = 0; i < MyProc*Chunk; i++)
        {
            std::getline(ifs,line);
        }

        if (!ifs.good())
        {
            std::string msg("ParticleContainer::InitFromAsciiFile(");
            msg += file;
            msg += ") failed @ 1";
            BoxLib::Error(msg.c_str());
        }

        int MyCnt = Chunk;

        if (MyProc == (NReaders - 1))
            //
            // We'll take the remainder.
            //
            MyCnt += cnt % NReaders;

        const Geometry& geom = m_amr->Geom(0);

        for (int i = 0; i < MyCnt; i++)
        {
            //
            // We don't read in m_id or m_cpu.  We'll set those later
            // in a manner to guarantee the global uniqueness of the pair.
            //
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                ifs >> p.m_pos[d];
            }

            //
            // If particle is right on a domain boundary then move it just inside the boundary.
            //
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                const Real domain_size = geom.ProbHi(d) - geom.ProbLo(d);

                if (p.m_pos[d] == geom.ProbLo(d)) p.m_pos[d] += 1.e-12 * domain_size;
                if (p.m_pos[d] == geom.ProbHi(d)) p.m_pos[d] -= 1.e-12 * domain_size;
            }

            for (int n = 0; n < extradata; n++)
            {
                ifs >> p.m_data[n];
            }

            if (!ifs.good())
            {
                std::string msg("ParticleContainer::InitFromAsciiFile(");
                msg += file;
                msg += ") failed @ 2";
                BoxLib::Error(msg.c_str());
            }

            if (!ParticleBase::Where(p,m_amr))
            {
                ParticleBase::PeriodicShift(p,m_amr);

                if (!ParticleBase::Where(p,m_amr))
                {
                    std::cout << "BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                    }

                    BoxLib::Abort("ParticleContainer<N>::InitFromAsciiFile(): invalid particle");
                }
            }

            p.m_id  = ParticleBase::NextID();
            p.m_cpu = MyProc;

            nparticles.push_back(p);

            howmany++;
            howmany_read++;

            Real domain_size[BL_SPACEDIM];
            for (int d=0; d<BL_SPACEDIM; ++d) {
                domain_size[d]  = (geom.ProbHi(d) - geom.ProbLo(d)) / float(lNrep[d]);
            }

            int rep[BL_SPACEDIM];
#if BL_SPACEDIM==3
            for (rep[2] = 1; rep[2] <= lNrep[2]; rep[2]++)
            {
#endif
                for (rep[1] = 1; rep[1] <= lNrep[1]; rep[1]++)
                {
                    for (rep[0] = 1; rep[0] <= lNrep[0]; rep[0]++) 
                    {
                        if (!(D_TERM( (rep[0] == 1), && (rep[1] == 1), && (rep[2] == 1) ) ) )
                        {
                            //
                            // Shift the position.
                            //
                            for (int d=0; d<BL_SPACEDIM; ++d) {
                                p_rep.m_pos[d] = p.m_pos[d] + float(rep[d]-1)*domain_size[d];
                            }
                            //
                            // Copy the mass and velocity.
                            //
                            for (int n = 0; n < extradata; n++)
                                p_rep.m_data[n] = p.m_data[n];

                            if (!ParticleBase::Where(p_rep,m_amr))
                            {
                                ParticleBase::PeriodicShift(p_rep,m_amr);
                                if (!ParticleBase::Where(p_rep,m_amr))
                                {
                                    std::cout << "BAD REPLICATED PARTICLE ID WOULD BE " << ParticleBase::NextID() << std::endl;
                                    BoxLib::Abort("ParticleContainer<N>::InitFromAsciiFile(): invalid replicated particle");
                                }
                            }
                            //
                            // Increment the particle ID.
                            //
                            p_rep.m_id  = ParticleBase::NextID();
                            //
                            // Assign to the same processor for now.
                            //
                            p_rep.m_cpu = MyProc;
   
                            nparticles.push_back(p_rep);
   
                            howmany++;
                        }
                    }
                }
#if BL_SPACEDIM==3
            }
#endif
        }
    }
    //
    // We've read in all the particles.
    // Now Redistribute() each chunk separately to minimize memory bloat.
    //
    int NRedist_chunk = NReaders / NRedist;

    for (int nr = 0; nr < NRedist; nr++)
    {
        if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Redistributing from processor " << nr*NRedist_chunk << " to " 
                                                          << (nr+1)*NRedist_chunk-1 << std::endl;
        for (int which = nr*NRedist_chunk; which < (nr+1)*NRedist_chunk; which++)
        {
            if (which == MyProc)
            {
                while (!nparticles.empty())
                {
                    ParticleType& p = nparticles.front();
    
                    m_particles[p.m_lev][p.m_grid].push_back(p);
    
                    nparticles.pop_front();
                }

                PBox().swap(nparticles);
            }
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute(true);
    }

    //
    // Take care of any leftover chunk
    //
    if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        if (NRedist*NRedist_chunk < NReaders)
            std::cout << "Redistributing from processor " << NRedist*NRedist_chunk << " to " 
                                                          << NReaders << std::endl;
    }
    for (int which = NRedist*NRedist_chunk; which < NReaders; which++)
    {
        if (which == MyProc)
        {
            while (!nparticles.empty())
            {
                ParticleType& p = nparticles.front();

                m_particles[p.m_lev][p.m_grid].push_back(p);

                nparticles.pop_front();
            }

            PBox().swap(nparticles);
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute(true);
    }

    //
    // Add up all the particles read in on each processor to get the total number of particles.
    //
    if (m_verbose > 0)
    {
        const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

        long num_particles = howmany; 

        ParallelDescriptor::ReduceLongSum(num_particles, IOProcNumber);

        if (D_TERM(lNrep[0] == 1, && lNrep[1] == 1, && lNrep[2] == 1))
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Total number of particles: " << num_particles << std::endl;
        }
        else
        {
            long num_particles_read = howmany_read; 

            ParallelDescriptor::ReduceLongSum(num_particles_read, IOProcNumber);

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "Replication the domain with vector           ";
                for (int d=0; d<BL_SPACEDIM; ++d) {
                    std::cout << lNrep[d] << " ";
                }
                std::cout << '\n';
                std::cout << "Total number of particles read in          : " << num_particles_read << '\n';
                std::cout << "Total number of particles after replication: " << num_particles      << std::endl;
            }
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor()) std::cout << "InitFromAsciiFile() time: " << runtime << std::endl;
    }
}

template <int N>
void
ParticleContainer<N>::InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_amr->Geom(0);
    const Real*     dx       = geom.CellSize();

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_amr->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    const Real   len[BL_SPACEDIM] = { geom.ProbLength(0),geom.ProbLength(1),geom.ProbLength(2) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        const int   fab_ix = fab_hi[0] - fab_lo[0] + 1;
        const int   fab_jx = fab_hi[1] - fab_lo[1] + 1;
        const int   fab_kx = fab_hi[2] - fab_lo[2] + 1;

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+0.5)*dx[n] +
		                      disp[n] * len[n];
                        //
			// Set the velocities.
                        //
	                p.m_data[n+1] = disp[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_amr))
                    {
      		        ParticleBase::PeriodicShift(p,m_amr);

                        if (!ParticleBase::Where(p,m_amr))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		    }

	            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_amr->finestLevel());
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
	            m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass)
{
    Real shift[] = {0,0,0};
    InitCosmo(mf, vel_fac, n_part, particleMass, shift);
}

template <int N>
void
ParticleContainer<N>::InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass, const Real shift[])
{
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_amr->Geom(0);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_amr->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { geom.ProbLength(0),geom.ProbLength(1),geom.ProbLength(2) };
    //
    // Print the grids as a sanity check.
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        if (mfi.validbox().isEmpty())
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Empty box in InitCosmo ");
        }
        if (!geom.Domain().contains(mfi.validbox()))
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Box in InitCosmo not contained in domain");
        }
    }
    //
    // We will need one ghost cell, so check wether we have one.
    //
    if (mf.nGrow() < 1)
        BoxLib::Abort("ParticleContainer<N>::InitCosmo: mf needs at least one correctly filled ghost zone!");

    if ( !(n_part[0] == n_part[1] && n_part[1] == n_part[2]) )
    {
	    std::cout << '\n';
	    std::cout << '\n';
	    std::cout << "Your particle lattice will have different spacings in the spatial directions!" << '\n';
	    std::cout << "You might want to change the particle number or the algorithm... ;)" << '\n';
	    std::cout << '\n';
	    std::cout << std::endl;
    }
    //
    // Place the particles evenly spaced in the problem domain.
    // Not perfectly fast - but easy
    //
    Real         pos[BL_SPACEDIM];
    ParticleType p;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box&  box     = mfi.validbox();
        RealBox     gridloc = RealBox(box, geom.CellSize(), geom.ProbLo());
	const Real* xlo     = gridloc.lo();
	const Real* xhi     = gridloc.hi();

        for (int k = 0; k < n_part[2]; k++)
        {
	    for (int j = 0; j < n_part[1]; j++)
            {
	        for (int i = 0; i < n_part[0]; i++)
                {
		    bool    isInValidBox = true;
            	    IntVect indices(D_DECL(i, j, k));

		    for (int n = 0; n < BL_SPACEDIM; n++)
                    {
  		        pos[n] = geom.ProbLo(n)
		               + (indices[n] + 0.5)*len[n]/n_part[n]
			       + shift[n];
                        //
			// Make sure particle is not on a boundary...
                        //
			pos[n] += 1e-14 * (geom.ProbHi(n) - geom.ProbLo(n));

			isInValidBox = isInValidBox 
				     && (pos[n] > xlo[n]) 
				     && (pos[n] < xhi[n]);
		    }

		    if (isInValidBox)
                    {
                        D_TERM(p.m_pos[0] = pos[0];,
                               p.m_pos[1] = pos[1];,
                               p.m_pos[2] = pos[2];);
                        //
		        // Set the mass of the particle.
                        //
	                p.m_data[0] = particleMass;
	                p.m_id      = ParticleBase::NextID();
	                p.m_cpu     = MyProc;

	                if (!ParticleBase::Where(p,m_amr))
                        {
      		            ParticleBase::PeriodicShift(p,m_amr);

                            if (!ParticleBase::Where(p,m_amr))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_amr->finestLevel());
	                //
	                // Add it to the appropriate PBox at the appropriate level.
	                //
	                m_particles[p.m_lev][p.m_grid].push_back(p);
		    }
	        }
	    }
        }
    }

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with equidistant placement" << std::endl;
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << std::endl;
    }

    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    PMap& pmap = m_particles[0];
    //
    // Make sure, that mf and m_amr.boxArray(0) are based on the same boxarray.
    //
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid    = pmap_it->first;
        PBox&            pbox    = pmap_it->second;
        const int        n       = pbox.size();
        const FArrayBox& dfab    = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            Real disp[BL_SPACEDIM];
            //
	    // Do CIC interpolation onto the particle positions.
	    // For CIC we need one ghost cell!
            //
            ParticleBase::GetGravity(dfab, m_amr, p, disp);

            D_TERM(p.m_pos[0] += len[0]*disp[0];,
                   p.m_pos[1] += len[1]*disp[1];,
                   p.m_pos[2] += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.m_data[1] = vel_fac[0]*disp[0];,
                   p.m_data[2] = vel_fac[1]*disp[1];,
                   p.m_data[3] = vel_fac[2]*disp[2];);

            if (!ParticleBase::Where(p,m_amr))
            {
	        ParticleBase::PeriodicShift(p,m_amr);

                if (!ParticleBase::Where(p,m_amr))
                    BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
	    }

            ParticleBase::Reset(p,m_amr,true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with particle displacement" << std::endl;
    }

    if (m_verbose > 1)
    {
        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "InitCosmo() done time: " << runtime << std::endl;
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitRandom (long          icount,
                                  unsigned long iseed,
                                  Real          mass,
                                  bool          serialize)
{
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_amr != 0);

    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       NProcs   = ParallelDescriptor::NProcs();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_amr->Geom(0);

    Real r, len[BL_SPACEDIM] = { geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2) };

    BoxLib::InitRandom(iseed+MyProc);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_amr->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    if (serialize)
    {
        //
        // We'll let IOProc generate the particles so we get the same
        // positions no matter how many CPUs we have.  This is here
        // mainly for debugging purposes.  It's not really useful for
        // very large numbers of particles.
        //
        //
        Array<Real> pos(icount*BL_SPACEDIM);

        if (ParallelDescriptor::IOProcessor())
        {
            for (long j = 0; j < icount; j++)
            {
                for (int i = 0; i < BL_SPACEDIM; i++)
                {
                    do
                    {
                        r = BoxLib::Random();
                    }
                    while (r == 0 || r == 1);

                    pos[j*BL_SPACEDIM + i] = r;
                }
            }
        }

        ParallelDescriptor::Bcast(pos.dataPtr(), icount*BL_SPACEDIM, IOProc);

        int cnt = 0;

        for (long j = 0; j < icount; j++)
        {
            ParticleType p;

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                r = pos[j*BL_SPACEDIM + i];

                p.m_pos[i] = geom.ProbLo(i) + (r * len[i]);

                BL_ASSERT(p.m_pos[i] < geom.ProbHi(i));
            }

            if (N > 0)
            {
                p.m_data[0] = mass;

                for (int i = 1; i < N; i++)
                    //
                    // Just zero out the rest of the data for lack of a better value.
                    //
                    p.m_data[i] = 0;
            }

            if (!ParticleBase::Where(p,m_amr))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_amr->finestLevel());

            const int who = m_amr->getLevel(p.m_lev).get_new_data(0).DistributionMap()[p.m_grid];

            if (who == MyProc)
            {
                //
                // We own it. Add it to the appropriate PBox at the appropriate level.
                //
                p.m_id  = ParticleBase::NextID();
                p.m_cpu = MyProc;

                m_particles[p.m_lev][p.m_grid].push_back(p);

                cnt++;
            }
        }

        BL_ASSERT(OK());
    }
    else
    {
        //
        // We'll generate the particles in parallel.
        //
        // Each CPU will key off the given seed to get independent streams of random numbers.
        //
        long M = icount / NProcs;
        //
        // Processor 0 will get the slop.
        //
        if (MyProc == 0)
        {
            M += (icount % NProcs);
        }

        for (long icnt = 0; icnt < M; icnt++)
        {
            ParticleType p;

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                do
                {
                    r = BoxLib::Random();
                }
                while (r == 0 || r == 1);

                p.m_pos[i] = geom.ProbLo(i) + (r * len[i]);

                BL_ASSERT(p.m_pos[i] < geom.ProbHi(i));
            }

            if (N > 0)
            {
                p.m_data[0] = mass;

                for (int i = 1; i < N; i++)
                    //
                    // Just zero out the rest of the data for lack of a better value.
                    //
                    p.m_data[i] = 0;
            }

            p.m_id  = ParticleBase::NextID();
            p.m_cpu = MyProc;

            if (!ParticleBase::Where(p,m_amr))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_amr->finestLevel());
            //
            // Add it to the appropriate PBox at the appropriate level.
            //
            m_particles[p.m_lev][p.m_grid].push_back(p);
        }
        //
        // Let Redistribute() sort out where the particles belong.
        //
        Redistribute(true);
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::InitRandom() time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::MoveRandom ()
{
    //
    // Move particles randomly at all levels
    //
    for (int lev = 0; lev < m_particles.size(); lev++)
    {
       MoveRandom(lev);
    }
}

template <int N>
void
ParticleContainer<N>::MoveRandom (int lev)
{
    BL_ASSERT(OK());
    BL_ASSERT(m_amr != 0);
    // 
    // Move particles up to FRAC*CellSize distance in each coordinate direction.
    //
    const Real FRAC = 0.25;

    static bool first = true;

    static Array<BoxLib::mt19937> rn;

    if (first)
    {
        first = false;
        //
        // Build and initialize a random number generator per thread.
        //
        int tnum = 1;

#ifdef _OPENMP
        tnum = omp_get_max_threads();
#endif
        rn.resize(tnum);

        for (int i = 0; i < tnum; i++)
        {
            //
            // We want to give each thread across all MPI processes a unique non-zero seed.
            //
            const unsigned long seedbase = 1+tnum*ParallelDescriptor::MyProc();

            rn[i] = BoxLib::mt19937(seedbase+i);
        }
    }

    PMap&       pmap              = m_particles[lev];
    const Real* dx                = m_amr->Geom(lev).CellSize();
    const Real  dist[BL_SPACEDIM] = { D_DECL(FRAC*dx[0], FRAC*dx[1], FRAC*dx[2]) };

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        PBox&     pbox = pmap_it->second;
        const int n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                p.m_pos[i] += dist[i]*(2*rn[tid].d_value()-1);
            }

            ParticleBase::Reset(p,m_amr,true);
        }
    }

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::Increment (MultiFab& mf,
                                 int       lev) 
{
    IncrementWithTotal(mf,lev);
}

template <int N>
long
ParticleContainer<N>::IncrementWithTotal (MultiFab& mf,
                                          int       lev)
{
    BL_ASSERT(OK());

    if (m_particles.empty()) return 0;

    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const PMap& pmap = m_particles[lev];
  
    long num_particles_in_domain = 0;

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int   grid = pmap_it->first;
        const PBox& pbox = pmap_it->second;
        FArrayBox&  fab  = mf[grid];

        for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
        {
            if (it->m_id > 0)
            {
                BL_ASSERT(it->m_grid == grid);

                fab(it->m_cell) += 1;
                num_particles_in_domain += 1;
            }
        }
    }

    ParallelDescriptor::ReduceLongSum(num_particles_in_domain);

    return num_particles_in_domain;
}

template <int N>
Real
ParticleContainer<N>::estTimestep (const MultiFab& gv,
                                   Real            a,
                                   int             lev,
                                   Real            cfl) const
{
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real      strttime         = ParallelDescriptor::second();
    const Geometry& geom             = m_amr->Geom(lev);
    const Real*     dx               = geom.CellSize();
    const Real      adx[BL_SPACEDIM] = { D_DECL(a*dx[0],a*dx[1],a*dx[2]) };
    Real            dt               = 1e50;
    const PMap&     pmap             = m_particles[lev];
    int             tnum             = 1;

#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif

    Array<Real> ldt(tnum,1e50);

    long num_particles_at_level = 0;

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        const PBox&      pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = gv[grid];

        num_particles_at_level += n;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            const Real mag_vel_over_dx[BL_SPACEDIM] = { D_DECL(abs(p.m_data[1])/adx[0],
                                                               abs(p.m_data[2])/adx[1],
                                                               abs(p.m_data[3])/adx[2]) };
            Real max_mag_vel_over_dx = mag_vel_over_dx[0];

#if (BL_SPACEDIM > 1)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[1], max_mag_vel_over_dx);
#endif
#if (BL_SPACEDIM > 2)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[2], max_mag_vel_over_dx);
#endif
            Real dt_part = (max_mag_vel_over_dx > 0) ? (cfl / max_mag_vel_over_dx) : 1e50;

            const Real gval[BL_SPACEDIM] = { D_DECL(gfab(p.m_cell,0),
                                                    gfab(p.m_cell,1),
                                                    gfab(p.m_cell,2)) };

            const Real mag_grav = sqrt(D_TERM(gval[0]*gval[0],
                                            + gval[1]*gval[1],
                                            + gval[2]*gval[2]));
            if (mag_grav > 0)
                dt_part = std::min( dt_part, 1/sqrt(mag_grav/dx[0]) );

            int tid = 0;

#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ldt[tid] = std::min(dt_part, ldt[tid]);
        }
    }

    for (int i = 0; i < ldt.size(); i++)
        dt = std::min(dt, ldt[i]);

    ParallelDescriptor::ReduceRealMin(dt);
    //
    // Set dt negative if there are no particles at this level.
    //
    ParallelDescriptor::ReduceLongSum(num_particles_at_level);

    if (num_particles_at_level == 0) dt = -1.e50;

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::estTimestep() time: " << stoptime << '\n';
        }
    }

    return dt;
}

//
// Assumes mass is in m_data[0]!
//

template <int N>
Real
ParticleContainer<N>::sumParticleMass (const MultiFab& mf,
                                       int             lev) const
{
    BL_ASSERT(N >= 1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    Real msum = 0;

    const PMap& pmap = m_particles[lev];

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const PBox& pbox = pmap_it->second;

        for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
        {
            if (it->m_id > 0)
            {
                msum += it->m_data[0];
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(msum);

    return msum;
}

//
// Assumes mass is in m_data[0], vx in m_dat[1], ...!
// dim defines the cartesian direction in which the momentum is summed, x is 0, y is 1, ...
//

template <int N>
void
ParticleContainer<N>::sumParticleMomentum (const MultiFab& mf,
                                           int             lev,
                                           Real*           mom) const
{
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const PMap& pmap = m_particles[lev];

    D_TERM(mom[0] = 0;, mom[1] = 0;, mom[2] = 0;);

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const PBox& pbox = pmap_it->second;
        const int   n    = pbox.size();

        Real mom_0 = 0, mom_1 = 0, mom_2 = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:mom_0,mom_1,mom_2)
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                D_TERM(mom_0 += p.m_data[0] * p.m_data[1];,
                       mom_1 += p.m_data[0] * p.m_data[2];,
                       mom_2 += p.m_data[0] * p.m_data[3];);
            }
        }
        
        D_TERM(mom[0] += mom_0;, mom[1] += mom_1;, mom[2] += mom_2;);
    }

    ParallelDescriptor::ReduceRealSum(mom,BL_SPACEDIM);
}

//
// This is the old single-level version.
//

template <int N>
void
ParticleContainer<N>::AssignDensity (MultiFab& mf,
                                     int       lev) const
{
    BL_ASSERT(N >= 1);
    BL_ASSERT(lev == 0);
    BL_ASSERT(mf.boxArray() == m_amr->boxArray(lev));

    const Real strttime = ParallelDescriptor::second();

    const Geometry& gm   = m_amr->Geom(lev);
    const Real*     dx   = gm.CellSize();
    const PMap&     pmap = m_particles[lev];
    const int       n    = pmap.size();
    const int       M    = D_TERM(2,+2,+4);
    const Real*     plo  = gm.ProbLo();

    if (gm.isAnyPeriodic() && !gm.isAllPeriodic())
        BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        mf[mfi].setVal(0);
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    Array<int>         pgrd(n);
    Array<const PBox*> pbxs(n);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if (n > 1)
#endif
    for (int i = 0; i < n; i++)
    {
        const PBox& pbx = *pbxs[i];
        FArrayBox&  fab = mf[pgrd[i]];

        bool    which[M];
        Real    fracs[M];
        IntVect cells[M];

        for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end(); it != End; ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0)
                //
                // Skip invalid particles.
                //
                continue;

            ParticleBase::AssignDensityCoeffs(p, plo, dx, fracs, cells);
            //
            // If this is not fully periodic then we have to be careful that no
            // particles go within one cell of the domain boundary.
            //
            const IntVect& csect = cells[0];

            if (!gm.isAllPeriodic())
                if (!gm.Domain().contains(csect))
                    BoxLib::Error("AssignDensityDoit: if not periodic, all particles must stay away from the domain boundary");

            for (int i = 0; i < M; i++)
            {
                fab(cells[i]) += p.m_data[0] * fracs[i];
            }
        }
    }

    mf.SumBoundary();
    gm.SumPeriodicBoundary(mf);

    const Real cell_vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf.mult(1/cell_vol);

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensity(single-level) time: " << stoptime << '\n';
        }
    }
}

//
// This is the multi-level version.
//
// The PArray should be empty on input.
//
// The MultiFabs in the PArray will be Managed'd on return.
//
// There'll be finest_level+1 of them.
//

template <int N>
void
ParticleContainer<N>::AssignDensity (PArray<MultiFab>& mf) const
{
    BL_ASSERT(N >= 1);

    const Real  strttime = ParallelDescriptor::second();

    const int   lfinest = m_amr->finestLevel();
    const Real* plo     = m_amr->Geom(0).ProbLo();
    const int   M       = D_TERM(2,+2,+4);

    mf.resize(lfinest+1, PArrayManage);

    for (int lev = 0; lev <= lfinest; lev++)
    {
        mf.set(lev, new MultiFab(m_amr->boxArray(lev), 1, 1));

        for (MFIter mfi(mf[lev]); mfi.isValid(); ++mfi)
        {
            mf[lev][mfi].setVal(0);
        }
    }

    for (int lev = 0; lev <= lfinest; lev++)
    {
        const Geometry& gm   = m_amr->Geom(lev);
        const Real*     dx   = gm.CellSize();
        const PMap&     pmap = m_particles[lev];
        const int       n    = pmap.size();
        const BoxArray& ba   = mf[lev].boxArray();
        BoxArray        cfba;

        if (gm.isAnyPeriodic() && !gm.isAllPeriodic())
            BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");

        if (lev < lfinest)
        {
            //
            // Build a coarsened fine BoxArray for Crse->Fine issue.
            //
            cfba = m_amr->boxArray(lev+1);
            cfba.coarsen(m_amr->refRatio(lev));
            BL_ASSERT(mf[lev].boxArray().contains(cfba));
        }
        //
        // This is a little funky.  What in effect this'll do is force
        // each thread to work on a single (separate) grid at a time.  That
        // way no thread will step on any other.  If there's only one grid per CPU,
        // then oh well ....
        //
        Array<int>         pgrd(n);
        Array<const PBox*> pbxs(n);

        int j = 0;
        for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it, ++j)
        {
            pgrd[j] =   pmap_it->first;
            pbxs[j] = &(pmap_it->second);
        }

//#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic,1) if (n > 1)
//#endif
        for (int i = 0; i < n; i++)
        {
            const PBox& pbx = *pbxs[i];
            FArrayBox&  fab =  mf[lev][pgrd[i]];
            const Box&  vbx =  ba[pgrd[i]];

            int      cgrid[M],  fgrid[M];
            bool    cwhich[M], fwhich[M];
            Real     fracs[M];
            IntVect  cells[M];

            for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end(); it != End; ++it)
            {
                const ParticleType& p = *it;

                if (p.m_id <= 0)
                    //
                    // Skip invalid particles.
                    //
                    continue;
                //
                // Get "fracs" and "cells" for the particle "p" at this level.
                //
                ParticleBase::AssignDensityCoeffs(p, plo, dx, fracs, cells);
                //
                // If this is not fully periodic then we have to be careful that no
                // particles go within one cell of the domain boundary.
                //
                if (!gm.isAllPeriodic())
                    if (!gm.Domain().contains(cells[0]))
                        BoxLib::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");

                const bool crsetofine = (lev < lfinest) ? ParticleBase::CrseToFine(cfba,cells,cwhich) : false;
                const bool finetocrse = (lev > 0      ) ? ParticleBase::FineToCrse(p,plo,dx,vbx,ba,m_amr->boxArray(lev-1),fwhich,cgrid) : false;
                //
                // Can crsetofine and finetocrse both be true simultaneously?
                //
                BL_ASSERT(!(crsetofine && finetocrse));

                if (!crsetofine && !finetocrse)
                {
                    //
                    // Just do it ...
                    //
                    for (int i = 0; i < M; i++)
                        fab(cells[i]) += p.m_data[0] * fracs[i];
                }

                if (finetocrse)
                {
                    //
                    // We have some fine->crse action.
                    //
                    // We first update only those fine ghost cells that cross
                    // fine->crse boundary using cfracs based on crse dx.
                    //
                    // This code assumes that ref ratio is 2 !!!
                    //
                    Real     sum_crse = 0, sum_fine = 0;;
                    Real     cfracs[M];
                    Real     cdx[BL_SPACEDIM] = { D_DECL(2*dx[0], 2*dx[1], 2*dx[2]) };
                    IntVect  ccells[M];

                    ParticleBase::AssignDensityCoeffs(p, plo, cdx, cfracs, ccells);

                    for (int i = 0; i < M; i++)
                    {
                        if (fwhich[i])
                        {
                            //
                            // Here we need to update the crse region.
                            // The coarse region is always going to be
                            // updated if we have a particle in a cell
                            // bordering a Fine->Crse boundary.  This
                            // is regardless of whether or not the particle
                            // is within dx/2 of the boundary.
                            //
                            // NOTE: This won't work in parallel.
                            //
                            BL_ASSERT(cgrid[i] >= 0);
                            BL_ASSERT(cgrid[i] < mf[lev-1].size());

                            FArrayBox& cfab = mf[lev-1][cgrid[i]];
                            //
                            // Here we need to update the crse region directly
                            // using cfracs and ccells.
                            //
                            cfab(ccells[i]) += p.m_data[0] * cfracs[i];
                            sum_crse        +=               cfracs[i];
                        }
                    }
                    //
                    // We've updated the Crse cells.  Now we have to update the fine
                    // cells in such a way that the total amount of mass we move
                    // around is precisely p.m_data[0]. In other words, the fractions
                    // we use at crse and fine have to sum to zero.  In the fine
                    // case, we have to account for the case where one or more of the
                    // cells indices is not in the valid region of the cell.
                    //
                    for (int i = 0; i < M; i++)
                    {
                        if (vbx.contains(cells[i]))
                        {
                            fwhich[i] = true;
                            sum_fine += fracs[i];
                        }
                        else
                        {
                            fwhich[i] = false;
                        }
                    }

                    BL_ASSERT(sum_fine > 0);

                    const Real mult = (1 - sum_crse) / sum_fine;

                    sum_fine = 0;

                    for (int i = 0; i < M; i++)
                    {
                        if (fwhich[i])
                        {
                            fab(cells[i]) += p.m_data[0] * fracs[i] * mult;
                            sum_fine      +=               fracs[i] * mult;
                        }
                    }

                    BL_ASSERT(std::abs(1-(sum_crse+sum_fine)) < 1.e-7);
                }

                if (crsetofine)
                {
                    Real sum_crse = 0, sum_fine = 0;

                    for (int i = 0; i < M; i++)
                    {
                        if (!cwhich[i])
                        {
                            fab(cells[i]) += p.m_data[0] * fracs[i];
                            sum_crse      +=               fracs[i];
                        }
                    }
                    //
                    // Next we have to find the Fine cells that need updating.
                    //
                    // NOTE: we may not own the Fine cells in parallel.
                    //
                    Real    ffracs[M];
                    IntVect fcells[M];

                    const Real*     fdx = m_amr->Geom(lev+1).CellSize();
                    const BoxArray& fba = mf[lev+1].boxArray();
                    //
                    // I want to use AssignDensityCoeffs() to find which fine
                    // cells need to be updated.  To do this I need to tweak
                    // the particle to within fine dx/2 of the C-F boundary.
                    //
                    ParticleBase pb = p;

                    RealBox rb(Box(pb.m_cell,pb.m_cell),dx,plo);

                    for (int i = 0; i < BL_SPACEDIM; i++)
                    {
                        if ((pb.m_pos[i] - rb.lo(i)) <= dx[i]/2)
                            //
                            // Put it close to the lo end so it sees fine grid.
                            //
                            pb.m_pos[i] = rb.lo(i) + fdx[i]/10;
                        else
                            //
                            // Put it close to the hi end so it sees fine grid.
                            //
                            pb.m_pos[i] = rb.hi(i) - fdx[i]/10;
                    }

                    ParticleBase::AssignDensityCoeffs(pb, plo, fdx, ffracs, fcells);

                    std::vector< std::pair<int,Box> > isects;

                    for (int i = 0; i < M; i++)
                    {
                        //
                        // Which grid at the fine level do we need to update?
                        //
                        isects = fba.intersections(Box(fcells[i],fcells[i]));

                        fwhich[i] = false;

                        if (!isects.empty())
                        {
                            BL_ASSERT(isects.size() == 1);
                            fwhich[i] = true;
                            sum_fine += ffracs[i];
                            fgrid[i]  = isects[0].first;
                        }
                    }

                    BL_ASSERT(sum_fine > 0);

                    const Real mult = (1 - sum_crse) / sum_fine;

                    sum_fine = 0;

                    for (int i = 0; i < M; i++)
                    {
                        if (fwhich[i])
                        {
                            BL_ASSERT(fgrid[i] >= 0);
                            BL_ASSERT(fgrid[i] < mf[lev+1].size());

                            FArrayBox& ffab  = mf[lev+1][fgrid[i]];

                            ffab(fcells[i]) += p.m_data[0] * ffracs[i] * mult;
                            sum_fine        +=               ffracs[i] * mult;
                        }
                    }

                    BL_ASSERT(std::abs(1-(sum_crse+sum_fine)) < 1.e-7);
                }
            }
        }

        mf[lev].SumBoundary();
        gm.SumPeriodicBoundary(mf[lev]);
    }

#if 0
    //
    // Next we have to sum down from fine to coarse only over
    // fine ghost cells that are not covered by fine valid region.
    //
    for (int i = lfinest; i > 0; i--)
    {
        //
        // First we make a temporary copy of the fine region
        // in which we set the valid region and any ghost cells
        // covered by valid region to zero.
        //
        MultiFab fmf(mf[i].boxArray(), 1, m_amr->MaxRefRatio(i-1));

        fmf.setVal(0);

        for (MFIter mfi(fmf); mfi.isValid(); ++mfi)
        {
            fmf[mfi].copy(mf[i][mfi]);

            fmf[mfi].setVal(0, mfi.validbox(), 0, 1);
        }

        fmf.FillBoundary();
        m_amr->Geom(i).FillPeriodicBoundary(fmf);
        //
        // Now sum down to a coarsen()d version of fmf.
        //
        // The only non-zero stuff in fmf is precisely what needs to be summed down.
        //
        BoxArray cba = fmf.boxArray();

        cba.coarsen(m_amr->refRatio(i-1));

        BL_ASSERT(mf[i-1].boxArray().contains(cba));

        cba.grow(1);

        MultiFab cmf(cba, 1, 0);

        cmf.setVal(0);

        for (MFIter mfi(cmf); mfi.isValid(); ++mfi)
        {
            FArrayBox&       cfab  = cmf[mfi];
            const Box&       ovlp  = cfab.box();
            const FArrayBox& ffab  = fmf[mfi];
            IntVect          ratio = m_amr->refRatio(i-1);

            BL_FORT_PROC_CALL(PART_SUMMASSDOWN,part_summassdown)
                (BL_TO_FORTRAN(cfab), 
                 BL_TO_FORTRAN(ffab),
                 ovlp.loVect(),ovlp.hiVect(),
                 ratio.getVect());
        }

        fmf.clear();
        //
        // We now have to sum cmf into mf[i-1].
        //
        MultiFab tcmf(mf[i-1].boxArray(), 1, 0);

        tcmf.setVal(0);

        tcmf.copy(cmf,FabArrayBase::ADD);

        for (MFIter mfi(tcmf); mfi.isValid(); ++mfi)
        {
            mf[i-1][mfi].plus(tcmf[mfi], mfi.validbox(), 0, 0, 1);
        }
        //
        // Make a version of cmf that does not have a grown boxArray but does have a ghost cell.
        // We need this in order to do SumPeriodicBoundary.
        //
        cba = mf[i].boxArray();

        cba.coarsen(m_amr->refRatio(i-1));

        MultiFab new_cmf(cba, 1, 1, Fab_allocate);

        for (MFIter mfi(new_cmf); mfi.isValid(); ++mfi)
        {
            new_cmf[mfi].copy(cmf[mfi]);
        }
        //
        // We now shift new_cmf periodically and sum into mf[i-1].
        //
        m_amr->Geom(i-1).SumPeriodicBoundary(mf[i-1], new_cmf);
    }
#endif

    for (int i = 0; i <= lfinest; i++)
    {
        const Real* dx       = m_amr->Geom(i).CellSize();
        const Real  cell_vol = D_TERM(dx[0], *dx[1], *dx[2]);

        mf[i].mult(1/cell_vol);
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensity(multi-level) time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::MultiplyParticleMass (int lev, Real mult)
{
   BL_ASSERT(lev == 0);
   PMap& pmap = m_particles[lev];
   for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
   {
       PBox&     pbx = pmap_it->second;
       const int n    = pbx.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
       for (int i = 0; i < n; i++)
       {
          ParticleType& p = pbx[i];
          if (p.m_id > 0)
          {
              //
              // Note: m_data[0] is mass, ...
              //
              p.m_data[0] *= mult;
          }
       }
   }
}

template <int N>
void
ParticleContainer<N>::movePredict (const MultiFab& gv,
                                   int             lev,
                                   Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime = ParallelDescriptor::second();

    PMap& pmap = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = gv[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                BL_ASSERT(p.m_grid == grid);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                D_TERM(p.m_data[1] += dt * gfab(p.m_cell,0);,
                       p.m_data[2] += dt * gfab(p.m_cell,1);,
                       p.m_data[3] += dt * gfab(p.m_cell,2););

                D_TERM(p.m_pos[0]  += dt * p.m_data[1];,
                       p.m_pos[1]  += dt * p.m_data[2];,
                       p.m_pos[2]  += dt * p.m_data[3];);

                ParticleBase::Reset(p,m_amr,true);
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::movePredict() time: " << stoptime << '\n';
        }
    }

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::moveCorrect (const MultiFab& gv_old,
                                   const MultiFab& gv,
                                   int             lev,
                                   Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime = ParallelDescriptor::second();

    PMap& pmap = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid     = pmap_it->first;
        PBox&            pbox     = pmap_it->second;
        const int        n        = pbox.size();
        const FArrayBox& gfab     = gv[grid];
        const FArrayBox& gfab_old = gv_old[grid];
        const Real       half_dt  = 0.5 * dt;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                BL_ASSERT(p.m_grid == grid);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                D_TERM(p.m_pos[0]  -= half_dt * p.m_data[1];,
                       p.m_pos[1]  -= half_dt * p.m_data[2];,
                       p.m_pos[2]  -= half_dt * p.m_data[3];);
            
                D_TERM(p.m_data[1] += half_dt * ( gfab(p.m_cell,0) - gfab_old(p.m_cell,0) );,
                       p.m_data[2] += half_dt * ( gfab(p.m_cell,1) - gfab_old(p.m_cell,1) );,
                       p.m_data[3] += half_dt * ( gfab(p.m_cell,2) - gfab_old(p.m_cell,2) ););

                D_TERM(p.m_pos[0]  += half_dt * p.m_data[1];,
                       p.m_pos[1]  += half_dt * p.m_data[2];,
                       p.m_pos[2]  += half_dt * p.m_data[3];);

                ParticleBase::Reset(p,m_amr,true);
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveCorrect() time: " << stoptime << '\n';
        }
    }

    Redistribute(true);
}

//
// This version takes as input the gravity vector at cell centers
//

template <int N>
void
ParticleContainer<N>::moveKickDrift (const MultiFab& grav_vector,
                                     int             lev,
                                     Real            dt)
{
    Real a_old  = 1.;
    Real a_half = 1.;
    moveKickDrift (grav_vector,lev,dt,a_old,a_half);
}

//
// This version takes as input the gravity vector at cell centers
//

template <int N>
void
ParticleContainer<N>::moveKickDrift (const MultiFab& grav_vector,
                                     int             lev,
                                     Real            dt,
                                     Real            a_old,
                                     Real            a_half) 
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime      = ParallelDescriptor::second();
    const Real half_dt       = 0.5 * dt;
    const Real a_half_inv    = 1 / a_half;
    const Real dt_a_half_inv = dt * a_half_inv;
    PMap&      pmap          = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = grav_vector[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
            Real grav[BL_SPACEDIM];

            ParticleBase::GetGravity(gfab, m_amr, p, grav);
            //
            // First update (a u)^half = (a u)^old + dt/2 grav^old
            //
            D_TERM(p.m_data[1] *= a_old;,
                   p.m_data[2] *= a_old;,
                   p.m_data[3] *= a_old;);
            //
            // Add adot/a and gravitational updates.
            //
            D_TERM(p.m_data[1] += half_dt * grav[0];,
                   p.m_data[2] += half_dt * grav[1];,
                   p.m_data[3] += half_dt * grav[2];);

            D_TERM(p.m_data[1] *= a_half_inv;,
                   p.m_data[2] *= a_half_inv;,
                   p.m_data[3] *= a_half_inv;);
            //
            // Now update x^new = x^old + dt grav^half / a^half
            //

            D_TERM(p.m_pos[0] += dt_a_half_inv * p.m_data[1];,
                   p.m_pos[1] += dt_a_half_inv * p.m_data[2];,
                   p.m_pos[2] += dt_a_half_inv * p.m_data[3];);

            ParticleBase::Reset(p,m_amr,true);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKickDrift() time: " << stoptime << '\n';
        }
    }

    Redistribute(true);
}

//
// This version takes as input the gravity vector at cell centers
//

template <int N>
void
ParticleContainer<N>::moveKick (const MultiFab& grav_vector,
                                int             lev,
                                Real            dt)
{
    Real a_new  = 1.;
    Real a_half = 1.;
    moveKick (grav_vector,lev,dt,a_new,a_half);
}

//
// This version takes as input the gravity vector at cell centers
//

template <int N>
void
ParticleContainer<N>::moveKick (const MultiFab& grav_vector,
                                int             lev,
                                Real            dt,
                                Real            a_new,
                                Real            a_half) 
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime  = ParallelDescriptor::second();
    const Real half_dt   = 0.5 * dt;
    const Real a_new_inv = 1 / a_new;
    PMap&      pmap      = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = grav_vector[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                BL_ASSERT(p.m_grid == grid);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                Real grav[BL_SPACEDIM];

                ParticleBase::GetGravity(gfab, m_amr, p, grav);
                //
                // Define (a u)^new = (a u)^half + dt/2 grav^new
                //
                D_TERM(p.m_data[1] *= a_half;,
                       p.m_data[2] *= a_half;,
                       p.m_data[3] *= a_half;);

                D_TERM(p.m_data[1] += half_dt * grav[0];,
                       p.m_data[2] += half_dt * grav[1];,
                       p.m_data[3] += half_dt * grav[2];);

                D_TERM(p.m_data[1] *= a_new_inv;,
                       p.m_data[2] *= a_new_inv;,
                       p.m_data[3] *= a_new_inv;);
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKick() time: " << stoptime << '\n';
        }
    }
    //
    // No need for Redistribution(), we only change the velocity.
    //
}

//
// This version takes as input the normal gravity component on each face
//

template <int N>
void
ParticleContainer<N>::moveKickDrift (PArray<MultiFab>& grav_vector,
                                     int             lev,
                                     Real            dt)
{
    Real a_old  = 1.;
    Real a_half = 1.;
    moveKickDrift (grav_vector,lev,dt,a_old,a_half);
}

//
// This version takes as input the normal gravity component on each face
//

template <int N>
void
ParticleContainer<N>::moveKickDrift (PArray<MultiFab>& grav_vector,
                                     int             lev,
                                     Real            dt,
                                     Real            a_old,
                                     Real            a_half) 
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real      strttime      = ParallelDescriptor::second();
    const Geometry& geom          = m_amr->Geom(lev);
    const Real*     dx            = geom.CellSize();
    const Real*     plo           = geom.ProbLo();
    const Real      half_dt       = 0.5 * dt;
    const Real      a_half_inv    = 1 / a_half;
    const Real      dt_a_half_inv = dt * a_half_inv;
    PMap&           pmap          = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int grid = pmap_it->first;
        PBox&     pbox = pmap_it->second;
        const int n    = pbox.size();

        const FArrayBox* gfab[BL_SPACEDIM] = { D_DECL(&grav_vector[0][grid],&grav_vector[1][grid],&grav_vector[2][grid]) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            //
            // First update (a u)^half = (a u)^old + dt/2 grav^old
            //
            D_TERM(p.m_data[1] *= a_old;,
                   p.m_data[2] *= a_old;,
                   p.m_data[3] *= a_old;);

            const IntVect lo = p.m_cell;

            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                IntVect hi = lo;

                hi[d] += 1;

                Real delta = (p.m_pos[d] - plo[d]) / dx[d] - lo[d];

                if (delta > 1) delta = 1;
                if (delta < 0) delta = 0;

                const Real grav_lo = (*gfab[d])(lo);
                const Real grav_hi = (*gfab[d])(hi);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                p.m_data[1+d] += half_dt * (grav_lo + delta * (grav_hi - grav_lo));
            }

            D_TERM(p.m_data[1] *= a_half_inv;,
                   p.m_data[2] *= a_half_inv;,
                   p.m_data[3] *= a_half_inv;);
            //
            // Now update x^new = x^old + dt grav^half / a^half
            //

            D_TERM(p.m_pos[0] += dt_a_half_inv * p.m_data[1];,
                   p.m_pos[1] += dt_a_half_inv * p.m_data[2];,
                   p.m_pos[2] += dt_a_half_inv * p.m_data[3];);

            ParticleBase::Reset(p,m_amr,true);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKickDrift() time: " << stoptime << '\n';
        }
    }

    Redistribute(true);
}

//
// This version takes as input the normal gravity component on each face
//

template <int N>
void
ParticleContainer<N>::moveKick (PArray<MultiFab>& grav_vector,
                                 int             lev,
                                 Real            dt)
{
    Real a_new  = 1.;
    Real a_half = 1.;
    moveKick (grav_vector,lev,dt,a_new,a_half);
}

//
// This version takes as input the normal gravity component on each face
//

template <int N>
void
ParticleContainer<N>::moveKick (PArray<MultiFab>& grav_vector,
                                int               lev,
                                Real              dt,
                                Real              a_new,
                                Real              a_half) 
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real      strttime  = ParallelDescriptor::second();
    const Geometry& geom      = m_amr->Geom(lev);
    const Real*     dx        = geom.CellSize();
    const Real*     plo       = geom.ProbLo();
    const Real      half_dt   = 0.5 * dt;
    const Real      a_new_inv = 1 / a_new;
    PMap&           pmap      = m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int grid = pmap_it->first;
        PBox&     pbox = pmap_it->second;
        const int n    = pbox.size();

        const FArrayBox* gfab[BL_SPACEDIM] = { D_DECL(&grav_vector[0][grid],&grav_vector[1][grid],&grav_vector[2][grid]) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            //
            // Define (a u)^new = (a u)^half + dt/2 grav^new
            //
            D_TERM(p.m_data[1] *= a_half;,
                   p.m_data[2] *= a_half;,
                   p.m_data[3] *= a_half;);

            const IntVect lo = p.m_cell;

            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                IntVect hi = lo;

                hi[d] += 1;

                Real delta = (p.m_pos[d] - plo[d]) / dx[d] - lo[d];

                if (delta > 1) delta = 1;
                if (delta < 0) delta = 0;

                const Real grav_lo = (*gfab[d])(lo);
                const Real grav_hi = (*gfab[d])(hi);
                //
                // Note: m_data[0] is mass, 1 is v_x, ...
                //
                p.m_data[1+d] += half_dt * (grav_lo + delta * (grav_hi - grav_lo));
            }

            D_TERM(p.m_data[1] *= a_new_inv;,
                   p.m_data[2] *= a_new_inv;,
                   p.m_data[3] *= a_new_inv;);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKick() time: " << stoptime << '\n';
        }
    }
    //
    // No need for Redistribution(), we only change the velocity.
    //
}

template <int N>
void
ParticleContainer<N>::RemoveParticlesNotAtFinestLevel ()
{
    BL_ASSERT(m_amr->finestLevel()+1 == m_particles.size());

    int cnt = 0;

    for (int lev = 0; lev < m_amr->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox& pbox = pmap_it->second;

            for (typename PBox::iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
            {
                ParticleType& p = *it;

                if (p.m_id > 0)
                {
                    //
                    // Negate it so that it's no longer considered a valid particle.
                    // It'll be removed on the next Redistribute().
                    //
                    p.m_id = -p.m_id;
                    cnt++;
                }
            }
        }
    }
    //
    // Print how many particles removed on each processor.
    //
    if (m_verbose > 1)
    {
        for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        {
            if (ParallelDescriptor::MyProc() == i)
            {
                if (cnt > 0)
                {
                    std::cout << "Processor "
                              << i
                              << " removed "
                              << cnt
                              << " particles not in finest level" << std::endl;
                }
            }
            ParallelDescriptor::Barrier();
        }
    }
}

//
// Uses a predictor/corrector to advance particles using umac.
//
//

template <int N>
void
ParticleContainer<N>::AdvectWithUmac (const MultiFab* umac,
                                      int             lev,
                                      Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= 2*BL_SPACEDIM);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    D_TERM(BL_ASSERT(umac[0].nGrow() >= 1);,
           BL_ASSERT(umac[1].nGrow() >= 1);,
           BL_ASSERT(umac[2].nGrow() >= 1););

    D_TERM(BL_ASSERT(!umac[0].contains_nan());,
           BL_ASSERT(!umac[1].contains_nan());,
           BL_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_amr->Geom(lev);
    const Real*     dx       = geom.CellSize();
    const Real*     plo      = geom.ProbLo();

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            const int grid = pmap_it->first;
            PBox&     pbox = pmap_it->second;
            const int n    = pbox.size();

            const FArrayBox* fab[BL_SPACEDIM] = { D_DECL(&umac[0][grid],&umac[1][grid],&umac[2][grid]) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];

                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

                const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + 0.5,
                                                       (p.m_pos[1]-plo[1])/dx[1] + 0.5,
                                                       (p.m_pos[2]-plo[2])/dx[2] + 0.5) };

                const IntVect hi(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

                const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-hi[0], len[1]-hi[1], len[2]-hi[2]) };

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    IntVect ehi = hi;

                    ehi[d] = p.m_cell[d] + 1;

                    Real efrac[BL_SPACEDIM] = { D_DECL(frac[0], frac[1], frac[2]) };

                    efrac[d] = (p.m_pos[d]-plo[d])/dx[d] - p.m_cell[d];

                    for (int j = 0; j < BL_SPACEDIM; j++)
                    {
                        if (efrac[j] > 1) efrac[j] = 1;
                        if (efrac[j] < 0) efrac[j] = 0;
                    }

                    const Real vel = ParticleBase::InterpDoit(*fab[d], ehi, efrac, 0);

                    if (ipass == 0)
                    {
                        //
                        // Predictor:
                        //
                        // Save old position and the vel & predict new time location.
                        //
                        p.m_data[d]             = p.m_pos[d];
                        p.m_data[d+BL_SPACEDIM] = p.m_pos[d];

                        p.m_pos[d] += dt*vel;
                    }
                    else
                    {
                        //
                        // Corrector:
                        //
                        // Update to final time using the orig position and the new vel.
                        //
                        // Save the velocity for use in Timestamp().
                        //
                        p.m_pos[d]  = 0.5 * (p.m_data[d] + p.m_pos[d] + dt*vel);

                        p.m_data[d] = 0.5 * (p.m_data[d+BL_SPACEDIM] + vel);
                    }
                }

                if (!ParticleBase::Where(p,m_amr,true))
                {
                    //
                    // If the particle goes out of the domain we invalidate it.
                    // Even if we're periodic in that direction.  Too much trouble
                    // to think about the right thing to do right now.
                    //
                    std::cout << "Invalidating out-of-domain particle: " << p << '\n';

                    BL_ASSERT(p.m_id > 0);

                    p.m_id = -p.m_id;
                }
            }
        }

        Redistribute(true);
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AdvectWithUmac() time: " << stoptime << '\n';
        }
    }
}

//
// This redistributes valid particles and discards invalid ones.
//

template <int N>
void
ParticleContainer<N>::Redistribute (bool where_already_called)
{
    const int MyProc   = ParallelDescriptor::MyProc();
    Real      strttime = ParallelDescriptor::second();

    if (m_particles.size() < m_amr->finestLevel()+1)
    {
       if (ParallelDescriptor::IOProcessor())
           std::cout << "ParticleContainer<N>::Redistribute() resizing m_particles from " << 
           m_particles.size() << " to " << m_amr->finestLevel()+1 << std::endl;
       m_particles.resize(m_amr->finestLevel()+1);
    }
    //
    // The valid particles that we don't own.
    //
    PBox notowned;
    //
    // And who rightly owns those particles.
    //
    std::deque<int> owner;

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            int   grid = pmap_it->first;
            PBox& pbox = pmap_it->second;

            //
            // Remove any invalid particles from front of container.
            //
            while (!pbox.empty() && pbox.front().m_id <= 0)
                pbox.pop_front();

            for (typename PBox::iterator it = pbox.begin(), End = pbox.end(); it != End; )
            {
                ParticleType& p = *it;

                if (p.m_id > 0)
                {
                    if (!where_already_called)
                        if (!ParticleBase::Where(p,m_amr))
                            BoxLib::Abort("ParticleContainer<N>::Redistribute(): invalid particle");
                    //
                    // The owner of the particle is the CPU owning the finest grid
                    // in state data that contains the particle.
                    //
                    const int who = m_amr->getLevel(p.m_lev).get_new_data(0).DistributionMap()[p.m_grid];

                    if (who == MyProc)
                    {
                        if (p.m_lev != lev || p.m_grid != grid)
                        {
                            //
                            // We own it but must shift it to another place.
                            //
                            m_particles[p.m_lev][p.m_grid].push_back(p);
                            //
                            // Invalidate the particle so we can reclaim its space.
                            //
                            p.m_id = -p.m_id;
                        }
                    }
                    else
                    {
                        owner.push_back(who);
                        notowned.push_back(p);
                        //
                        // Invalidate the particle so we can reclaim its space.
                        //
                        p.m_id = -p.m_id;
                    }
                }

                if (p.m_id <= 0)
                {
                    if (it != pbox.begin())
                    {
                        BL_ASSERT(pbox.front().m_id > 0);
                        std::swap(pbox.front(),p);
                    }
                    ++it;
                    pbox.pop_front();
                }
                else
                {
                    ++it;
                }
            }
        }
        //
        // Remove any map entries for which the particle container is now empty.
        //
        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; )
        {
            if (pmap_it->second.empty())
                pmap.erase(pmap_it++);
            else
                ++pmap_it;
        }
    }

    if (m_particles.size() > m_amr->finestLevel()+1)
    {
        //
        // Looks like we lost an AmrLevel on a regrid.
        //
        if (ParallelDescriptor::IOProcessor())
            std::cout << "ParticleContainer<N>::Redistribute() resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << m_amr->finestLevel()+1 << std::endl;

        BL_ASSERT(m_particles.size() >= 2);

        BL_ASSERT(m_particles[m_particles.size()-1].empty());

        m_particles.resize(m_amr->finestLevel()+1);
    }

    const int NProcs = ParallelDescriptor::NProcs();

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (NProcs == 1)
            ByteSpread();

        if (ParallelDescriptor::IOProcessor())
            std::cout << "ParticleContainer<N>::Redistribute()   scalar time: " << stoptime << '\n';
    }

    if (NProcs == 1)
    {
        BL_ASSERT(OK());
        BL_ASSERT(notowned.empty());
        return;
    }

#if BL_USE_MPI

    strttime = ParallelDescriptor::second();
    //
    // We may now have particles that are rightfully owned by another CPU.
    //
    Array<int> Snds(NProcs,0);
    Array<int> Rcvs(NProcs,0);

    BL_ASSERT(notowned.size() == owner.size());

    for (std::deque<int>::const_iterator it = owner.begin(), End = owner.end(); it != End; ++it)
        Snds[*it]++;

    long maxsendcount = 0;
    for (int i = 0; i < NProcs; i++)
        maxsendcount += Snds[i];
    ParallelDescriptor::ReduceLongMax(maxsendcount);

    if (maxsendcount == 0)
    {
        //
        // There's no parallel work to do.
        //
        BL_ASSERT(OK());
        return;
    }

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    int NumRcvs = 0;
    for (int i = 0; i < NProcs; i++)
        NumRcvs += Rcvs[i];

    int NumSnds = 0;
    for (int i = 0; i < NProcs; i++)
        NumSnds += Snds[i];

    BL_ASSERT(notowned.size() == NumSnds);

    PBox nparticles;

    {
        //
        // The "int" data.
        //
        const int iChunkSize = 4 + BL_SPACEDIM;

        Array<int> recvdata (NumRcvs*iChunkSize);
        Array<int> senddata (NumSnds*iChunkSize);

        Array<int> sendcnts(NProcs), recvcnts(NProcs);

        Array<int> sdispls(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * iChunkSize;
            sendcnts[i] = Snds[i] * iChunkSize;

            if (i > 0)
            {
                offset [i] = offset [i-1] + sendcnts[i-1];
                rdispls[i] = rdispls[i-1] + recvcnts[i-1];
                sdispls[i] = sdispls[i-1] + sendcnts[i-1];
            }
        }

        int i = 0;

        for (typename PBox::const_iterator it = notowned.begin(), End = notowned.end(); it != End; ++it)
        {
            int& ioff = offset[owner[i]];

            BL_ASSERT(it->m_id > 0);

            senddata[ioff+0] = it->m_id;
            senddata[ioff+1] = it->m_cpu;

            senddata[ioff+2] = it->m_lev;
            senddata[ioff+3] = it->m_grid;

            D_TERM(senddata[ioff+4] = it->m_cell[0];,
                   senddata[ioff+5] = it->m_cell[1];,
                   senddata[ioff+6] = it->m_cell[2];);

            ioff += iChunkSize;

            i++;
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      NumRcvs == 0 ? 0 : recvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );
        Array<int>().swap(senddata);
        //
        // Unpack data into the new particles.
        //
        ParticleType p;

        int* rcvp = NumRcvs == 0 ? 0 : recvdata.dataPtr();

        for (int i = 0; i < NumRcvs; i++)
        {
            BL_ASSERT(rcvp != 0);

            p.m_id   = rcvp[0];
            p.m_cpu  = rcvp[1];

            p.m_lev  = rcvp[2];
            p.m_grid = rcvp[3];

            D_TERM(p.m_cell[0] = rcvp[4];,
                   p.m_cell[1] = rcvp[5];,
                   p.m_cell[2] = rcvp[6];);

            nparticles.push_back(p);

            rcvp += iChunkSize;
        }
    }

    BL_ASSERT(nparticles.size() == NumRcvs);

    {
        //
        // The "Real" data (m_pos & m_data).
        //
        const int rChunkSize = BL_SPACEDIM+N;

        Array<Real> recvdata (NumRcvs*rChunkSize);
        Array<Real> senddata (NumSnds*rChunkSize);

        Array<int> sendcnts(NProcs), recvcnts(NProcs);

        Array<int> sdispls(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * rChunkSize;
            sendcnts[i] = Snds[i] * rChunkSize;

            if (i > 0)
            {
                offset [i] = offset [i-1] + sendcnts[i-1];
                rdispls[i] = rdispls[i-1] + recvcnts[i-1];
                sdispls[i] = sdispls[i-1] + sendcnts[i-1];
            }
        }

        int i = 0;

        for (typename PBox::const_iterator it = notowned.begin(), End = notowned.end(); it != End; ++it)
        {
            int& ioff = offset[owner[i]];

            D_TERM(senddata[ioff+0] = it->m_pos[0];,
                   senddata[ioff+1] = it->m_pos[1];,
                   senddata[ioff+2] = it->m_pos[2];);

            ioff += BL_SPACEDIM;

            for (int j = 0; j < N; j++)
                senddata[ioff+j] = it->m_data[j];

            ioff += N;

            i++;
        }

        PBox().swap(notowned);

        std::deque<int>().swap(owner);

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      NumRcvs == 0 ? 0 : recvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      ParallelDescriptor::Communicator()) );
        Array<Real>().swap(senddata);
        //
        // Unpack data into the new particles.
        //
        Real* rcvp = NumRcvs == 0 ? 0 : recvdata.dataPtr();

        while (!nparticles.empty())
        {
            BL_ASSERT(rcvp != 0);

            ParticleType& p = nparticles.front();

            D_TERM(p.m_pos[0] = rcvp[0];,
                   p.m_pos[1] = rcvp[1];,
                   p.m_pos[2] = rcvp[2];);

            rcvp += BL_SPACEDIM;

            for (int j = 0; j < N; j++)
                p.m_data[j] = rcvp[j];

            rcvp += N;

            m_particles[p.m_lev][p.m_grid].push_back(p);

            nparticles.pop_front();
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        ByteSpread();

        if (ParallelDescriptor::IOProcessor())
            std::cout << "ParticleContainer<N>::Redistribute() parallel time: " << stoptime << '\n';
    }
#endif
}

template <int N>
bool
ParticleContainer<N>::OK () const
{
    //
    // Check that the integer data in each valid particle is what it should be.
    // This includes checking that particles are in the proper place in the particle
    // container based on what Where() says they should be.
    //
    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            const int   grid = pmap_it->first;
            const PBox& pbox = pmap_it->second;

            for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
            {
                //
                // Yes I want to make a copy of the particle.
                //
                ParticleType p = *it;

                if (p.m_id > 0)
                {
                    const int     llev  = p.m_lev;
                    const int     lgrid = p.m_grid;
                    const IntVect cell  = p.m_cell;

                    if (!ParticleBase::Where(p,m_amr))
                        return false;

                    if ((lev  != p.m_lev  || lev  != llev)  ||
                        (grid != p.m_grid || grid != lgrid) || cell != p.m_cell)
                    {
                        std::cout << "PARTICLE NUMBER " << p.m_id << '\n';

                        std::cout << "POS IS ";
                        for (int i = 0; i < BL_SPACEDIM; i++)
                            std::cout << p.m_pos[i] << ' ';

                        if (lev != p.m_lev || lev != llev)
                           std::cout << "BAD LEV  " << lev  << " " << p.m_lev << '\n';

                        if (grid != p.m_grid || grid != lgrid)
                           std::cout << "BAD GRID " << grid << " " << p.m_grid << '\n';

                        if (cell != p.m_cell)
                           std::cout << "BAD CELL " << cell << " " << p.m_cell << '\n';

                        return false;
                    }
                }
            }
        }
    }

    return true;
}

template <int N>
void
ParticleContainer<N>::Checkpoint (const std::string& dir,
                                  const std::string& name) const
{
    BL_ASSERT(OK());

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // We store the particles in a subdirectory of "dir".
    //
    std::string pdir = dir;

    if (!pdir.empty() && pdir[pdir.size()-1] != '/')
        pdir += '/';

    pdir += name;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(pdir, 0755))
            BoxLib::CreateDirectoryFailed(pdir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the I/O processor writes to the header file.
    //
    std::ofstream HdrFile;

    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator pmap_it = pmap.begin(), End = pmap.end(); pmap_it != End; ++pmap_it)
        {
            const PBox& pbox = pmap_it->second;

            for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
            {
                if (it->m_id > 0)
                {
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceLongSum(nparticles,IOProc);

    int maxnextid = ParticleBase::NextID();

    ParticleBase::NextID(maxnextid);

    ParallelDescriptor::ReduceIntMax(maxnextid,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;

        if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
            HdrFileName += '/';

        HdrFileName += "Header";

        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!HdrFile.good())
            BoxLib::FileOpenFailed(HdrFileName);
        //
        // First thing written is our Checkpoint/Restart version string.
        //
        HdrFile << ParticleBase::Version() << '\n';
        //
        // BL_SPACEDIM and N for sanity checking.
        //
        HdrFile << BL_SPACEDIM << '\n';

        HdrFile << N << '\n';
        //
        // The total number of particles.
        //
        HdrFile << nparticles << '\n';
        //
        // The value of nextid that we need to restore on restart.
        //
        HdrFile << maxnextid << '\n';
        //
        // Then the finest level of the AMR hierarchy.
        //
        HdrFile << m_amr->finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= m_amr->finestLevel(); lev++)
        {
            HdrFile << m_amr->boxArray(lev).size() << '\n';
        }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    const int nOutFiles = std::min(64,NProcs);

    for (int lev = 0; lev <= m_amr->finestLevel(); lev++)
    {
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;

        if (!LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/')
            LevelDir += '/';

        LevelDir = BoxLib::Concatenate(LevelDir + "Level_", lev, 1);

        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(LevelDir, 0755))
                BoxLib::CreateDirectoryFailed(LevelDir);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();

        const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Array<int>  which(state.size(),0);
        Array<int > count(state.size(),0);
        Array<long> where(state.size(),0);

        const int   FileNumber   = MyProc % nOutFiles;
        std::string FullFileName = LevelDir;

        FullFileName += '/';
        FullFileName += ParticleBase::DataPrefix();
        FullFileName += BoxLib::Concatenate("", FileNumber, 4);

        std::ofstream ParticleFile;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        ParticleFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        const int nSets = ((NProcs + (nOutFiles - 1)) / nOutFiles);
        const int mySet = (MyProc / nOutFiles);

        for (int iSet = 0; iSet < nSets; ++iSet)
        {
            if (mySet == iSet)
            {
                //
                // Write all the data at this level to the file.
                //
                if (iSet == 0)
                    //
                    // First set.
                    //
                    ParticleFile.open(FullFileName.c_str(),
                                      std::ios::out|std::ios::trunc|std::ios::binary);
                else
                {
                    ParticleFile.open(FullFileName.c_str(),
                                      std::ios::out|std::ios::app|std::ios::binary);
                    //
                    // Set to the end of the file.
                    //
                    ParticleFile.seekp(0, std::ios::end);
                }

                if (!ParticleFile.good())
                    BoxLib::FileOpenFailed(FullFileName);
                //
                // Write out all the valid particles we own at the specified level.
                // Do it grid block by grid block remembering the seek offset
                // for the start of writing of each block of data.
                //
                WriteParticles(lev, ParticleFile, FileNumber, which, count, where);

                ParticleFile.flush();

                ParticleFile.close();

                if (!ParticleFile.good())
                    BoxLib::Abort("ParticleContainer<N>::Checkpoint(): problem writing ParticleFile");

                int iBuff = 0, wakeUpPID = (MyProc + nOutFiles), tag = (MyProc % nOutFiles);

                if (wakeUpPID < NProcs)
                {
                    ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
                }
            }

            if (mySet == (iSet + 1))
            {
                //
                // Next set waits.
                //
                int iBuff, waitForPID = (MyProc - nOutFiles), tag = (MyProc % nOutFiles);

                ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
            }
        }

        ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProc);
        ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProc);
        ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            for (int j = 0; j < state.size(); j++)
            {
                //
                // We now write the which file, the particle count, and the
                // file offset into which the data for each grid was written,
                // to the header file.
                //
                HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
            }
            //
            // Unlink any zero-length data files.
            //
            Array<long> cnt(nOutFiles,0);

            for (int i = 0; i < count.size(); i++)
                cnt[which[i]] += count[i];

            for (int i = 0; i < cnt.size(); i++)
            {
                if (cnt[i] == 0)
                {
                    std::string FullFileName = LevelDir;

                    FullFileName += '/';
                    FullFileName += ParticleBase::DataPrefix();
                    FullFileName += BoxLib::Concatenate("", i, 4);

                    BoxLib::UnlinkFile(FullFileName.c_str());
                }
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            HdrFile.flush();

            HdrFile.close();

            if (!HdrFile.good())
                BoxLib::Abort("ParticleContainer<N>::Checkpoint(): problem writing HdrFile");

            std::cout << "ParticleContainer<N>::Checkpoint() time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::WriteParticles (int            lev,
                                      std::ofstream& ofs,
                                      int            fnum,
                                      Array<int>&    which,
                                      Array<int>&    count,
                                      Array<long>&   where) const
{
    const PMap&     pmap  = m_particles[lev];
    const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();
        //
        // Only write out valid particles.
        //
        int cnt = 0;

        typename PMap::const_iterator pmap_it = pmap.find(grid);

        if (pmap_it != pmap.end())
        {
            const PBox& pbox = pmap_it->second;

            for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
            {
                if (it->m_id > 0)
                    cnt++;
            }
        }

        which[grid] = fnum;
        count[grid] = cnt;
        where[grid] = VisMF::FileOffset(ofs);

        if (cnt == 0) continue;

        const PBox& pbox = pmap_it->second;

        {
            //
            // First write out the integer data in binary.
            // We do not need to write out the m_lev and m_grid
            // info since it's implicit in how the particles
            // are stored.  We can easily recreate them on restart.
            //
            const int iChunkSize = 2+BL_SPACEDIM;

            Array<int> istuff(cnt*iChunkSize);

            int* iptr = istuff.dataPtr();

            for (typename PBox::const_iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
            {
                if (it->m_id > 0)
                {
                    BL_ASSERT(it->m_lev == lev);
                    BL_ASSERT(it->m_grid == grid);

                    iptr[0] = it->m_id;
                    iptr[1] = it->m_cpu;

                    D_TERM(iptr[2] = it->m_cell[0];,
                           iptr[3] = it->m_cell[1];,
                           iptr[4] = it->m_cell[2];);

                    iptr += iChunkSize;
                }
            }

            ofs.write((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
        }

        {
            //
            // Then the Real data in binary.
            //
            const int rChunkSize = BL_SPACEDIM+N;

            Array<Real> rstuff(cnt*rChunkSize);

            Real* rptr = rstuff.dataPtr();

            for (typename PBox::const_iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
            {
                if (it->m_id > 0)
                {
                    D_TERM(rptr[0] = it->m_pos[0];,
                           rptr[1] = it->m_pos[1];,
                           rptr[2] = it->m_pos[2];);

                    for (int i = 0; i < N; i++)
                        rptr[BL_SPACEDIM+i] = it->m_data[i];

                    rptr += rChunkSize;
                }
            }

            ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(Real));
        }
    }
}

template <int N>
void
ParticleContainer<N>::Restart (const std::string& dir,
                               const std::string& file)
{
    BL_ASSERT(!dir.empty());
    BL_ASSERT(!file.empty());

    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();

    std::string fullname = dir;

    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';

    fullname += file;
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the IO processor reads the header file.
    //
    // It'll then broadcast() stuff of interest to all CPUs.
    //
    std::ifstream HdrFile;

    std::string HdrFileName = fullname;

    if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
        HdrFileName += '/';

    HdrFileName += "Header";

    HdrFile.open(HdrFileName.c_str(), std::ios::in);

    if (!HdrFile.good())
        BoxLib::FileOpenFailed(HdrFileName);
    //
    // First value should be the version string.
    //
    Array<char> vbuf(128);

    std::string version;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> version;

        BL_ASSERT(!version.empty());
        BL_ASSERT(vbuf.size() > version.size());

        for (int i = 0; i < version.size(); i++)
            vbuf[i] = version[i];

        vbuf[version.size()] = '\0';
    }
    ParallelDescriptor::Bcast(vbuf.dataPtr(), vbuf.size(), IOProc);

    version = vbuf.dataPtr();

    if (version == ParticleBase::Version())
    {
        Restart_OneDotZero(fullname,HdrFile);
    }
    else
    {
        std::string msg("ParticleContainer<N>::Restart(): unknown version string: ");
        msg += version;
        BoxLib::Abort(msg.c_str());
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Restart() time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::Restart_OneDotZero (const std::string& fullname,
                                          std::ifstream&     HdrFile)
{
    BL_ASSERT(!fullname.empty());

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Next value should be BL_SPACEDIM;
    //
    int dm;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> dm;

        if (dm != BL_SPACEDIM)
            BoxLib::Abort("ParticleContainer<N>::Restart(): dm != BL_SPACEDIM");
    }
    ParallelDescriptor::Bcast(&dm, 1, IOProc);
    //
    // Next value should be our "N".
    //
    int n;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> n;

        if (n != N)
            BoxLib::Abort("ParticleContainer<N>::Restart(): n != N");
    }
    ParallelDescriptor::Bcast(&n, 1, IOProc);

    long nparticles;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The total number of particles.
        //
        HdrFile >> nparticles;

        BL_ASSERT(nparticles >= 0);
    }
    ParallelDescriptor::Bcast(&nparticles, 1, IOProc);

    int maxnextid;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The value of nextid that we need to restore.
        //
        HdrFile >> maxnextid;

        BL_ASSERT(maxnextid > 0);
    }
    ParallelDescriptor::Bcast(&maxnextid, 1, IOProc);
    //
    // Don't forget to restore it!!!
    //
    ParticleBase::NextID(maxnextid);
    //
    // Then the finest level of the AMR hierarchy.
    //
    int finest_level;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> finest_level;

        BL_ASSERT(finest_level >= 0);
    }
    ParallelDescriptor::Bcast(&finest_level, 1, IOProc);
    //
    // Then the number of grids at each level.
    //
    Array<int> ngrids(finest_level+1);

    BL_ASSERT(finest_level == m_amr->finestLevel());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            HdrFile >> ngrids[lev];

            BL_ASSERT(ngrids[lev] > 0);
            BL_ASSERT(ngrids[lev] == m_amr->boxArray(lev).size());
        }
    }
    ParallelDescriptor::Bcast(ngrids.dataPtr(), ngrids.size(), IOProc);
    //
    // The rest of HdrFile consists of triples of the form:
    //
    //   which count offset
    //
    // One for each grid at each level from 0 -> finest_level.
    //
    // We rebuild the filename from which and level.
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //
        // Read in the which, count & offset info for this level.
        //
        Array<int>  which(ngrids[lev]);
        Array<int>  count(ngrids[lev]);
        Array<long> where(ngrids[lev]);

        if (ParallelDescriptor::IOProcessor())
        {
            for (int i = 0; i < ngrids[lev]; i++)
            {
                HdrFile >> which[i] >> count[i] >> where[i];
            }
        }
        ParallelDescriptor::Bcast(which.dataPtr(), which.size(), IOProc);
        ParallelDescriptor::Bcast(count.dataPtr(), count.size(), IOProc);
        ParallelDescriptor::Bcast(where.dataPtr(), where.size(), IOProc);

        m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

        m_particles.resize(m_amr->finestLevel()+1);

        const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const int grid = mfi.index();

            if (count[grid] <= 0) continue;
            //
            // The file names in the header file are relative.
            //
            std::string name = fullname;

            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';

            name += "Level_";
            name += BoxLib::Concatenate("", lev, 1);
            name += '/';
            name += ParticleBase::DataPrefix();
            name += BoxLib::Concatenate("", which[grid], 4);

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in);

            if (!ParticleFile.good())
                BoxLib::FileOpenFailed(name);

            ParticleFile.seekg(where[grid], std::ios::beg);

            ReadParticles_OneDotZero(count[grid],grid,lev,ParticleFile);

            ParticleFile.close();

            if (!ParticleFile.good())
                BoxLib::Abort("ParticleContainer<N>::Restart(): problem reading particles");
        }
    }

    BL_ASSERT(OK());
}

template <int N>
void
ParticleContainer<N>::ReadParticles_OneDotZero (int            cnt,
                                                int            grd,
                                                int            lev,
                                                std::ifstream& ifs)
{
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < m_particles.size());
    BL_ASSERT(lev >= 0 && lev <= m_amr->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_amr->boxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2+BL_SPACEDIM;

    Array<int> istuff(cnt*iChunkSize);

    ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the Real data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+N;

    Array<Real> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(Real));
    //
    // Now reassemble the particles.
    //
    int*         iptr = istuff.dataPtr();
    Real*        rptr = rstuff.dataPtr();
    PBox&        pbox = m_particles[lev][grd];
    ParticleType p;

    for (int i = 0; i < cnt; i++)
    {
        p.m_id   = iptr[0];
        p.m_cpu  = iptr[1];
        p.m_lev  = lev;
        p.m_grid = grd;

        BL_ASSERT(p.m_id > 0);

        BL_ASSERT(p.m_lev  == lev);
        BL_ASSERT(p.m_grid == grd);

        D_TERM(p.m_cell[0] = iptr[2];,
               p.m_cell[1] = iptr[3];,
               p.m_cell[2] = iptr[4];);

        iptr += iChunkSize;

        D_TERM(p.m_pos[0] = rptr[0];,
               p.m_pos[1] = rptr[1];,
               p.m_pos[2] = rptr[2];);

        for (int i = 0; i < N; i++)
            p.m_data[i] = rptr[BL_SPACEDIM+i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

template <int N>
void
ParticleContainer<N>::Timestamp (const std::string&      basename,
                                 const MultiFab&         mf,
                                 int                     lev,
                                 Real                    time,
                                 const std::vector<int>& indices)
{
    BL_ASSERT(lev >= 0);
    BL_ASSERT(time >= 0);
    BL_ASSERT(!basename.empty());
    BL_ASSERT(N >= 2*BL_SPACEDIM);
    BL_ASSERT(lev <= m_amr->finestLevel());

    const Real strttime = ParallelDescriptor::second();
    //
    // We'll spread the output over this many files.
    //
    // Should this be ParmParse'd in?
    //
    const int MaxOutFiles = 32;

    const int   MyProc    = ParallelDescriptor::MyProc();
    const int   NProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::min(MaxOutFiles, NProcs);
    const int   nSets     = ((NProcs + (nOutFiles - 1)) / nOutFiles);
    const int   mySet     = (MyProc / nOutFiles);

    for (int iSet = 0; iSet < nSets; ++iSet)
    {
        if (mySet == iSet)
        {
            //
            // Do we have any particles at this level that need writing?
            //
            bool gotwork = false;

            const PMap& pmap = m_particles[lev];

            for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
                 pmap_it != pmapEnd && !gotwork;
                 ++pmap_it)
            {
                for (typename PBox::const_iterator it = pmap_it->second.begin(), pboxEnd = pmap_it->second.end();
                     it != pboxEnd && !gotwork;
                     ++it)
                {
                    if (it->m_id > 0)
                    {
                        gotwork = true;
                    }
                }
            }

            if (gotwork)
            {
                std::string FileName = BoxLib::Concatenate(basename + '_', MyProc % nOutFiles, 2);

                std::ofstream TimeStampFile;

                VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

                TimeStampFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

                TimeStampFile.open(FileName.c_str(), std::ios::out|std::ios::app|std::ios::binary);

                TimeStampFile.setf(std::ios_base::scientific,std::ios_base::floatfield);

                TimeStampFile.precision(10);

                TimeStampFile.seekp(0, std::ios::end);

                if (!TimeStampFile.good())
                    BoxLib::FileOpenFailed(FileName);

                std::vector<Real> vals(indices.size());

                for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
                     pmap_it != pmapEnd;
                     ++pmap_it)
                {
                    const int        grid = pmap_it->first;
                    const PBox&      pbox = pmap_it->second;
                    const FArrayBox& fab  = mf[grid];

                    for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
                         it != pboxEnd;
                         ++it)
                    {
                        const ParticleType& p = *it;

                        if (p.m_id <= 0) continue;

                        BL_ASSERT(p.m_lev == lev);
                        BL_ASSERT(p.m_grid == grid);

                        TimeStampFile << p.m_id  << ' ' << p.m_cpu << ' ';

                        for (int i = 0; i < BL_SPACEDIM; i++)
                        {
                            TimeStampFile << p.m_pos[i] << ' ';
                        }

                        TimeStampFile << time;
                        //
                        // AdvectWithUmac stores the velocity in m_data ...
                        //
                        for (int i = 0; i < BL_SPACEDIM; i++)
                        {
                            TimeStampFile << ' ' << p.m_data[i];
                        }

                        if (!indices.empty())
                        {
                            ParticleBase::Interp(p,m_amr,fab,&indices[0],&vals[0],indices.size());

                            for (int i = 0; i < indices.size(); i++)
                            {
                                TimeStampFile << ' ' << vals[i];
                            }
                        }

                        TimeStampFile << '\n';
                    }
                }

                TimeStampFile.flush();
                TimeStampFile.close();
            }

            const int iBuff     = 0;
            const int wakeUpPID = (MyProc + nOutFiles);
            const int tag       = (MyProc % nOutFiles);

            if (wakeUpPID < NProcs)
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
        if (mySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int       iBuff;
            const int waitForPID = (MyProc - nOutFiles);
            const int tag        = (MyProc % nOutFiles);

            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Timestamp: lev: " << lev << " time: " << stoptime << '\n';
        }
    }
}
template <int N>
void
ParticleContainer<N>::WriteAsciiFile (const std::string& filename)
{
    BL_ASSERT(!filename.empty());

    const Real strttime = ParallelDescriptor::second();
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator pmap_it = pmap.begin(), End = pmap.end(); pmap_it != End; ++pmap_it)
        {
            const PBox& pbox = pmap_it->second;

            for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
            {
                if (it->m_id > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }
    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle count.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            BoxLib::FileOpenFailed(filename);

        File << nparticles << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            BoxLib::Abort("ParticleContainer<N>::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (MyProc == i)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            std::ofstream File;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                BoxLib::FileOpenFailed(filename);
            
            for (int lev = 0; lev < m_particles.size(); lev++)
            {
                const PMap& pmap = m_particles[lev];

                for (typename PMap::const_iterator pmap_it = pmap.begin(), End = pmap.end(); pmap_it != End; ++pmap_it)
                {
                    const PBox& pbox = pmap_it->second;

                    for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end(); it != pboxEnd; ++it)
                    {
                        if (it->m_id > 0)
                        {
                            D_TERM(File << it->m_pos[0] << ' ',
                                        << it->m_pos[1] << ' ',
                                        << it->m_pos[2] << ' ');

                            for (int i = 0; i < N; i++)
                            {
                                char ws = (i == N-1) ? '\n' : ' ';

                                File << it->m_data[i] << ws;
                            }
                        }
                    }
                }
            }

            File.flush();

            File.close();

            if (!File.good())
                BoxLib::Abort("ParticleContainer<N>::WriteAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::WriteAsciiFile() time: " << stoptime << '\n';
        }
    }
}

#endif /*_PARTICLES_H_*/
