#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>

#include <ParmParse.H>

#include <ParGDB.H>
#include <REAL.H>
#include <IntVect.H>
#include <Array.H>
#include <Utility.H>
#include <Geometry.H>
#include <VisMF.H>
#include <Particles_F.H>
#include <RealBox.H>
#include <BL_CXX11.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
    std::string aggregation_type   = "";
    int         aggregation_buffer = 1;
    const int   GhostParticleID    = 2000000000;
    const int   VirtualParticleID  = 1000000000;
}

struct ParticleBase
{
    //
    // The floating point type used for particle positions
    // and the "extra" data in the particles themselves.
    //
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif

    int      m_id;
    int      m_cpu;
    int      m_lev;
    int      m_grid;
    IntVect  m_cell;
    RealType m_pos[BL_SPACEDIM];

    ParticleBase ()
        :
        m_id(-1),
        m_cpu(-1),
        m_lev(-1),
        m_grid(-1)
        {}

    static IntVect Index (const ParticleBase& p, const Geometry& geom);
    //
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    //
    static bool Where (ParticleBase& prt, const ParGDBBase* gdb, int lev_min = 0, int finest_level = -1);
    //
    // Checks/sets whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    //
    static bool PeriodicWhere (ParticleBase& prt, const ParGDBBase* gdb, int lev_min = 0, int finest_level = -1);
    //
    // Checks/sets whether a particle is within its grid (including grow cells).
    //
    static bool RestrictedWhere (ParticleBase& p, const ParGDBBase* gdb, int ngrow);
    //
    // Checks/sets a particle's location on a specific level.
    // (Yes this is distict from the functionality provided above)
    //
    static bool SingleLevelWhere (ParticleBase& p, const ParGDBBase* gdb, int level);
    //
    // Updates a particle's location (Where), tries to periodic shift any particles
    // that have left the domain. May need work (see inline comments)
    //
    static void Reset (ParticleBase& prt, const ParGDBBase* gdb, bool update, bool verbose=true); 
    //
    // Returns true if the particle was shifted.
    //
    static bool PeriodicShift (ParticleBase& prt, const ParGDBBase* gdb);

    static Real InterpDoit (const FArrayBox& fab, const Real* fracs, const IntVect* cells, int comp);

    static Real InterpDoit (const FArrayBox& fab, const IntVect& hi, const Real* frac, int comp);

    static void Interp (const ParticleBase& prt, const Geometry& geom, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Geometry& geom, const ParticleBase& p, Real* grav);

    static int MaxReaders ();

    static long MaxParticlesPerRead ();
    //
    // Returns the next particle ID for this processor.
    // Particle IDs start at 1 and are never reused.
    // The pair, consisting of the ID and the CPU on which the particle is "born",
    // is a globally unique identifier for a particle.  The maximum of this value
    // across all processors must be checkpointed and then restored on restart
    // so that we don't reuse particle IDs.
    //
    static int NextID ();
    // This version can only be used inside omp critical.
    static int UnprotectedNextID ();
    //
    // Reset on restart.
    //
    static void NextID (int nextid);
    //
    // Used by AssignDensity.
    //
    static bool CrseToFine (const BoxArray&       cfba, 
                            const Array<IntVect>& cells, 
                            Array<IntVect>&       cfshifts, 
                            const Geometry&       gm, 
                            Array<int>&           which, 
                            Array<IntVect>&       pshifts);

    static bool FineToCrse (const ParticleBase&                p, 
                            int                                flev, 
                            const ParGDBBase*                  gdb, 
                            const Array<IntVect>&              fcells, 
                            const BoxArray&                    fvalid, 
                            const BoxArray&                    compfvalid_grown, 
                            Array<IntVect>&                    ccells, 
                            Array<Real>&                       cfracs, 
                            Array<int>&                        which, 
                            Array<int>&                        cgrid, 
                            Array<IntVect>&                    pshifts, 
                            std::vector< std::pair<int,Box> >& isects);

    static void FineCellsToUpdateFromCrse (const ParticleBase&                p, 
                                           int lev, const ParGDBBase*         gdb, 
                                           const IntVect&                     ccell,
                                           const IntVect&                     cshift, 
                                           Array<int>&                        fgrid, 
                                           Array<Real>&                       ffrac, 
                                           Array<IntVect>&                    fcells, 
                                           std::vector< std::pair<int,Box> >& isects);

    static void CIC_Fracs (const Real* frac, Real* fracs);

    static void CIC_Cells (const IntVect& hicell, IntVect* cells);
    //
    // Old, *-based CIC for use in Interp.
    //
    static void CIC_Cells_Fracs_Basic (const ParticleBase& p, const Real* plo, const Real* dx, Real* fracs,  IntVect* cells);
    //
    // Wraps the arbitrary dx function.
    //
    static int CIC_Cells_Fracs (const ParticleBase& p, 
                                const Real*         plo, 
                                const Real*         dx, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Does CIC computations for arbitrary particle/grid dx's.
    //
    static int CIC_Cells_Fracs (const ParticleBase& p, 
                                const Real*         plo, 
                                const Real*         dx_geom, 
                                const Real*         dx_part, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Useful for sorting particles into lexicographic order of their cell position.
    //
    class Compare
    {
    public:
        bool operator () (const ParticleBase& lhs,
                          const ParticleBase& rhs) const
        {
            return lhs.m_cell.lexLT(rhs.m_cell);
        }
    };
};

std::ostream& operator<< (std::ostream& os, const ParticleBase& p);

template <int N>
struct Particle
    :
    public ParticleBase
{
    //
    // The amount of floating point data we hold.
    //
    // In some cases this is:
    //
    // 0 - particle mass
    // 1 - x-velocity
    // 2 - y-velocity
    // 3 - z-velocity
    //
    RealType m_data[N];
};

//
// A concrete base class for ParticleContainer.
//
// This will allow us to then maintain, say, a vector of pointers to
// ParticleContainerBases that are actually pointers to instantiations of
// of different ParticleContainers.  For every function in ParticleContainer
// that you'd like to call using a pointer to a ParticleContainerBase you'll
// need to add a pure virtual function of that exact signature to 
// ParticleContainerBase.  Make sure that the signatures of virtual functions
// are exact including defaults arguments!
//

struct ParticleContainerBase
{
    virtual ~ParticleContainerBase ();

    virtual void AssignDensity (PArray<MultiFab>& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const = 0;

    virtual void AssignDensitySingleLevel      (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const = 0;
    virtual void AssignCellDensitySingleLevel  (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const = 0;
    virtual void AssignNodalDensitySingleLevel (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const = 0;

    virtual Real sumParticleMass (int level) const = 0;

    virtual void RemoveParticlesAtLevel (int level) = 0;

    virtual void Redistribute (bool where_already_called = false,
                               bool full_where           = false,
                               int  lev_min              = 0, 
                               int  nGrow                = 0) = 0;

    virtual void moveKickDrift (MultiFab& grav_vector, int level, Real timestep, 
                                Real a_old = 1.0, Real a_half = 1.0) = 0;
    virtual void moveKick      (MultiFab& grav_vector, int level, Real timestep, 
                                Real a_new = 1.0, Real a_half = 1.0, int start_comp_for_accel = -1) = 0;
    virtual void moveKickDrift (PArray<MultiFab>& grav_vector, int level, Real timestep, 
                                Real a_old = 1.0, Real a_half = 1.0) = 0;
    virtual void moveKick      (PArray<MultiFab>& grav_vector, int level, Real timestep, 
                                Real a_new = 1.0, Real a_half = 1.0) = 0;
};

template <int N>
class ParticleContainer
    :
    public ParticleContainerBase
{
public:
    //
    // The type of Particles we hold.
    //
    typedef Particle<N> ParticleType;
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

    ParticleContainer (ParGDBBase* gdb)
        :
        m_verbose(1), m_relativistic(0), m_csq(-1.), m_gdb(gdb), allow_particles_near_boundary(false) { }

    ParticleContainer (const Geometry            & geom, 
		       const DistributionMapping & dmap,
		       const BoxArray            & ba)
	:
        m_verbose(1), m_relativistic(0), m_csq(-1.),
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba)
    {
	m_gdb = & m_gdb_object;
    }

    ParticleContainer (const Array<Geometry>            & geom, 
		       const Array<DistributionMapping> & dmap,
		       const Array<BoxArray>            & ba,
		       const Array<int>                 & rr)
	:
        m_verbose(1), m_relativistic(0), m_csq(-1.),
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba,rr)
    {
	m_gdb = & m_gdb_object;
    }

    void SetParticleBoxArray (int lev,
			      const DistributionMapping& new_dmap,
			      const BoxArray           & new_ba)
    { m_gdb->SetParticleBoxArray(lev, new_dmap, new_ba); }

    const BoxArray& ParticleBoxArray (int lev) const 
	{ return m_gdb->ParticleBoxArray(lev); }

    const DistributionMapping& ParticleDistributionMap (int lev) const 
	{ return m_gdb->ParticleDistributionMap(lev); }

    const ParGDBBase* GetParGDB () const { return m_gdb; }

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);

    void InitFromBinaryFile (const std::string& file, int extradata);

    void InitFromBinaryMetaFile (const std::string& file, int extradata);

    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false, RealBox bx = RealBox());
    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass);
    void InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, 
                     const Real particleMass, const Real shift[]);
    void InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass);
    void InitCosmo1ppcMultiLevel(MultiFab& mf, const Real disp_fac[], const Real vel_fac[], const Real particleMass, 
                                 int disp_idx, int vel_idx, BoxArray &baWhereNot, int lev);

    void addOneParticle (int id_in, int cpu_in, 
                         std::vector<double>& xloc, std::vector<double>& attributes); 

    virtual Real sumParticleMass (int level) const BL_OVERRIDE;
    void sumParticleMomentum (int lev, Real* mom) const;

    void Increment (MultiFab& mf, int level);

    long IncrementWithTotal (MultiFab& mf, int level);

    virtual void AssignDensitySingleLevel      (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const BL_OVERRIDE;
    virtual void AssignCellDensitySingleLevel  (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const BL_OVERRIDE;
    virtual void AssignNodalDensitySingleLevel (MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const BL_OVERRIDE;

    virtual void AssignDensity (PArray<MultiFab>& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const BL_OVERRIDE;

    void AssignDensityAndVels (PArray<MultiFab>& mf, int lev_min = 0) const;

    void AssignDensityDoit (PArray<MultiFab>* mf, PMap& data, int ncomp, int lev_min = 0) const;

    void MultiplyParticleMass (int lev, Real mult);

    void GetParticleIDs        (Array<int> & part_ids);
    void GetParticleCPU        (Array<int> & part_cpu);
    void GetParticleLocations  (Array<Real>& part_locs);
    void GetParticleVelocities (Array<Real>& part_vels);
    void GetParticleData       (Array<Real>& part_data, int start_comp, int num_comp);

    void SetParticleLocations  (Array<Real>& part_data);
    void SetParticleVelocities (Array<Real>& part_data);

    // Set the flag that allows particles to live near the domain boundary and throw away 
    //     the part of their contribution in AssignDensity that is outside the domain.
    void SetAllowParticlesNearBoundary(bool value);
 
    Real estTimestep (MultiFab& grav_vector,         int level, Real cfl) const;
    Real estTimestep (MultiFab& grav_vector, Real a, int level, Real cfl) const;

    virtual void Redistribute (bool where_already_called = false,
			       bool full_where           = false,
			       int  lev_min              = 0, 
			       int  nGrow                = 0) BL_OVERRIDE;

    void RedistributeMPI (PMap& not_ours);
    //
    // OK checks that all particles are in the right places (for some value of right)
    //
    // These flags are used to do proper checking for subcycling particles
    // the default values are fine for non-subcycling methods
    //
    bool OK (bool full_where = false, int lev_min = 0 , int ngrow = 0, int finest_level = -1) const;

    void ByteSpread () const;
    //
    // Returns # of particles at specified the level.
    //
    // If "only_valid" is true it only counts valid particles.
    //
    long NumberOfParticlesAtLevel (int level, bool only_valid = true, bool only_local = false) const;
    Array<long> NumberOfParticlesInGrid  (int level, bool only_valid = true, bool only_local = false) const;
    //
    // Returns # of particles at all levels
    //
    // If "only_valid" is true it only counts valid particles.
    //
    long TotalNumberOfParticles (bool only_valid=true, bool only_local=false) const;

    void MoveRandom ();

    void MoveRandom (int level);

    // **********************************************************************************
    // Nyx Specifc Methods
    // **********************************************************************************
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
    //   comoving coordinates -- these rely on CELL-BASED gravity component
    //
    virtual void moveKickDrift (MultiFab& grav_vector, int level, Real timestep, 
				Real a_old = 1.0, Real a_half = 1.0) BL_OVERRIDE;
    virtual void moveKick      (MultiFab& grav_vector, int level, Real timestep, 
				Real a_new = 1.0, Real a_half = 1.0,
				int start_comp_for_accel = -1) BL_OVERRIDE;
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   comoving coordinates -- these rely on EDGE-BASED gravity component
    //
    virtual void moveKickDrift (PArray<MultiFab>& grav_vector, int level, Real timestep, 
				Real a_old = 1.0, Real a_half = 1.0) BL_OVERRIDE;
    virtual void moveKick      (PArray<MultiFab>& grav_vector, int level, Real timestep, 
				Real a_new = 1.0, Real a_half = 1.0) BL_OVERRIDE;
    //
    // after the moveKickDrift step the positions of the particles are advanced for a full timestep,
    // so this scheme should work in the overall algorithm...
    //
    //
    // The Following methods are for managing Nyx's Virtual and Ghost Particles.
    //
    // Removes all particles at a given level
    //
    virtual void RemoveParticlesAtLevel (int level) BL_OVERRIDE;
    //
    // Creates virtual particles for a given level that represent
    // in some capacity all particles at finer levels
    //
    void CreateVirtualParticles (int level, PBox& virts) const;
    // 
    // Create ghost particles for a given level that are copies of particles
    // near coarse->fine boundaries in level-1
    //
    void CreateGhostParticles (int level, int ngrow, PBox& ghosts) const;
    //
    // Add particles from a pbox to the grid at this level
    //
    void AddParticlesAtLevel (int level, PBox& virts, bool where_already_called = false);
    // **************************************************************************************************************** 

    void AdvectWithUmac (MultiFab* umac, int level, Real dt, int vcomp = 0);

    void AdvectWithUcc (const MultiFab& ucc, int level, Real dt, int vcomp = 0);

    void RemoveParticlesNotAtFinestLevel ();

    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint = true) const;

    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint = true);

    void WritePlotFile (const std::string& dir, const std::string& name) const;

    void Timestamp (const std::string& file, const MultiFab& mf, int lev, Real time, const std::vector<int>& idx, const int vcomp = 0);

    void WriteAsciiFile          (const std::string& file);
    void WriteCoarsenedAsciiFile (const std::string& file);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (Real csq) { m_csq = csq; }

    const PMap& GetParticles(int lev) const { return m_particles[lev]; }
    PMap& GetParticles(int lev) { return m_particles[lev]; }

protected:

    bool OnSameGrids (int level, const MultiFab& mf) const { return m_gdb->OnSameGrids(level, mf); }

    //
    // Helper function for Checkpoint() and WritePlotFile().
    //
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Array<int>&    which,
                         Array<int>&    count,
                         Array<long>&   where,
                         bool           is_checkpoint) const;
    //
    // Helper functions for Restart().
    //
    void Restart_Doit (const std::string& fullname,
                       std::ifstream&     HdrFile,
                       const std::string& how,
                       bool           is_checkpoint);

    void ReadParticles_DoublePrecision (int            cnt,
                                        int            grd,
                                        int            lev,
                                        bool           is_checkpoint,
                                        std::ifstream& ifs);

    void ReadParticles_SinglePrecision (int            cnt,
                                        int            grd,
                                        int            lev,
                                        bool           is_checkpoint,
                                        std::ifstream& ifs);
    //
    // The data.
    //
    int         m_verbose;
    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    Real        m_csq;          // the square of the speed of light -- used to compute relativistic effects


    ParGDBBase* m_gdb;
    bool allow_particles_near_boundary;
    ParGDB      m_gdb_object;
    Array<PMap> m_particles;

};

template <int N>
long
ParticleContainer<N>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    long nparticles = 0;
    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,only_local);
    return nparticles;
}

template <int N>
Array<long>
ParticleContainer<N>::NumberOfParticlesInGrid (int  lev, 
                                               bool only_valid,
                                               bool only_local) const
{
    int ngrids = m_gdb->ParticleBoxArray(lev).size();
    Array<long> nparticles(ngrids,0);

    if (lev >= 0 && lev < m_particles.size())
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator pmap_it = pmap.begin(), End = pmap.end();
             pmap_it != End;
             ++pmap_it)
        {
            int         pgrd = pmap_it->first;
            const PBox& pbox = pmap_it->second;

            nparticles[pgrd] = 0;

            if (only_valid)
            {
                for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
                     it != pboxEnd;
                     ++it)
                {
                    if (it->m_id > 0)
                        nparticles[pgrd]++;
                }
            }
            else
            {
                nparticles[pgrd] = pbox.size();
            }
        }

        if (!only_local)
            ParallelDescriptor::ReduceLongSum(&nparticles[0],ngrids);
    }

    return nparticles;
}

template <int N>
long
ParticleContainer<N>::NumberOfParticlesAtLevel (int  lev,
                                                bool only_valid,
                                                bool only_local) const
{
    long nparticles = 0;

    if (lev >= 0 && lev < m_particles.size())
    {
        const PMap& pmap = m_particles[lev];

        for (typename PMap::const_iterator pmap_it = pmap.begin(), End = pmap.end();
             pmap_it != End;
             ++pmap_it)
        {
            const PBox& pbox = pmap_it->second;

            if (only_valid)
            {
                for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
                     it != pboxEnd;
                     ++it)
                {
                    if (it->m_id > 0)
                        nparticles++;
                }
            }
            else
            {
                nparticles += pbox.size();
            }
        }
    }

    if (!only_local)
        ParallelDescriptor::ReduceLongSum(nparticles);

    return nparticles;
}

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

    long mn = cnt, mx = mn;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const std::size_t sz = sizeof(ParticleType);

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
    ParallelDescriptor::ReduceLongMin(mn, IOProc);
    ParallelDescriptor::ReduceLongMax(mx, IOProc);
    ParallelDescriptor::ReduceLongSum(cnt,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "ParticleContainer<N> byte spread across MPI nodes: ["
                  << mn*sz
		  << " (" << mn << ")"
                  << " ... "
                  << mx*sz
		  << " (" << mx << ")"
                  << "] total particles: (" << cnt << ")\n";
    }
#ifdef BL_LAZY
    });
#endif
}

template <int N>
void
ParticleContainer<N>::addOneParticle (int                  id_in,
                                      int                  cpu_in, 
                                      std::vector<double>& xloc, 
 				      std::vector<double>& attributes)
{
    BL_PROFILE("ParticleContainer<N>::addOneParticle()");
    if (m_particles.size() == 0)
    {
       m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
       m_particles.resize(m_gdb->finestLevel()+1);
    }

    ParticleType p;

    p.m_id  = id_in;
    p.m_cpu = cpu_in;

    if (p.m_id <= 0)
        BoxLib::Abort("Particle ID's must be > 0 in addOneParticle");
 
    if (ParallelDescriptor::MyProc() != p.m_cpu)
        BoxLib::Abort("cpu_in must equal MyProc() in addOneParticle");

    for (int i = 0; i < BL_SPACEDIM; i++)
       p.m_pos[i] = xloc[i];

    for (int i = 0; i < N; i++)
       p.m_data[i] = attributes[i];
    //
    // It's possible that particles are defined to live right on the boundary,
    // so if that's the case we move them ever so slightly inside the domain 
    // instead.
    //
    const Geometry& geom = m_gdb->Geom(0);

    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };
    //
    // If particle is right on a domain boundary then move it just inside the boundary.
    //
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (p.m_pos[d] <= geom.ProbLo(d)) p.m_pos[d] += Delta[d];
        if (p.m_pos[d] >= geom.ProbHi(d)) p.m_pos[d] -= Delta[d];
    }

    if (!ParticleBase::Where(p,m_gdb))
    {
        ParticleBase::PeriodicShift(p,m_gdb);

        if (!ParticleBase::Where(p,m_gdb))
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
            }
            BoxLib::Abort("ParticleContainer<N>::addOneParticle(): invalid particle");
        }
    }

    m_particles[p.m_lev][p.m_grid].push_back(p);

    // Note that we will need to call Redistribute once we are done adding particles this way.
    // The Where call above assigns a particle to a grid (p.m_grid) based on the particle location.
    // However, the particle may not currently live on the processor that owns that grid.
    // The Redistribute routine should ensure that the particle ends up on the right processor.
}

template <int N>
void
ParticleContainer<N>::InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep)
{
    BL_PROFILE("ParticleContainer<N>::InitFromAsciiFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= N);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const Real strttime = ParallelDescriptor::second();
    //
    // Number of processes that read from the file.  //
    int NReaders = ParticleBase::MaxReaders();
    //
    // Number of chunks we break the redistribution loop into
    // This default is just for the case of NProcs = NReaders = 1.
    //
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
    m_particles.resize(m_gdb->finestLevel()+1);

    IntVect lNrep(D_DECL(1,1,1));

    if (Nrep != 0)
        lNrep = *Nrep;

    long how_many      = 0;
    long how_many_read = 0;
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

        ifs >> cnt >> std::ws;

        ParticleType p, p_rep;

        const int Chunk = cnt / NReaders;

        for (int i = 0; i < MyProc*Chunk; i++)
        {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            ifs >> std::ws;  // Eat newline.
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

        const Geometry& geom = m_gdb->Geom(0);

        const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
        const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
        const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                                  Real(.125)*geom.CellSize(1),
                                                  Real(.125)*geom.CellSize(2)) };
        for (int i = 0; i < MyCnt; i++)
        {
            //
            // We don't read in m_id or m_cpu.  We'll set those later
            // in a manner to guarantee the global uniqueness of the pair.
            //
            D_TERM(ifs >> p.m_pos[0];,
                   ifs >> p.m_pos[1];,
                   ifs >> p.m_pos[2];);
            //
            // If particle is right on a domain boundary then move it just inside the boundary.
            //
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                if (p.m_pos[d] <= ProbLo[d]) p.m_pos[d] += Delta[d];
                if (p.m_pos[d] >= ProbHi[d]) p.m_pos[d] -= Delta[d];
            }

            for (int n = 0; n < extradata; n++)
            {
                ifs >> p.m_data[n];
            }

            if (!ifs.good())
            {
                std::string msg("ParticleContainer::InitFromAsciiFile(");
                msg += file; msg += ") failed @ 2";
                BoxLib::Error(msg.c_str());
            }

            if (!ParticleBase::Where(p,m_gdb))
            {
                ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
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

            how_many++;
            how_many_read++;

            const Real DomSize[BL_SPACEDIM] = { D_DECL((geom.ProbHi(0)-geom.ProbLo(0))/lNrep[0],
                                                       (geom.ProbHi(1)-geom.ProbLo(1))/lNrep[1],
                                                       (geom.ProbHi(2)-geom.ProbLo(2))/lNrep[2]) };
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
                            for (int d=0; d<BL_SPACEDIM; ++d)
                            {
                                p_rep.m_pos[d] = p.m_pos[d] + (rep[d]-1)*DomSize[d];
                            }
                            //
                            // Copy the mass and velocity.
                            //
                            for (int n = 0; n < extradata; n++)
                                p_rep.m_data[n] = p.m_data[n];

                            if (!ParticleBase::Where(p_rep,m_gdb))
                            {
                                ParticleBase::PeriodicShift(p_rep,m_gdb);
                                if (!ParticleBase::Where(p_rep,m_gdb))
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
   
                            how_many++;
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
                                                          << (nr+1)*NRedist_chunk-1 << '\n';
        for (int which = nr*NRedist_chunk; which < (nr+1)*NRedist_chunk; which++)
        {
            if (which == MyProc)
            {
                while (!nparticles.empty())
                {
                    ParticleType& p = nparticles.back();
    
                    m_particles[p.m_lev][p.m_grid].push_back(p);
    
                    nparticles.pop_back();
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
                                                          << NReaders << '\n';
    }
    for (int which = NRedist*NRedist_chunk; which < NReaders; which++)
    {
        if (which == MyProc)
        {
            while (!nparticles.empty())
            {
                ParticleType& p = nparticles.back();

                m_particles[p.m_lev][p.m_grid].push_back(p);

                nparticles.pop_back();
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

        long num_particles = how_many; 

        ParallelDescriptor::ReduceLongSum(num_particles, IOProcNumber);

        if (D_TERM(lNrep[0] == 1, && lNrep[1] == 1, && lNrep[2] == 1))
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Total number of particles: " << num_particles << '\n';
        }
        else
        {
            long num_particles_read = how_many_read; 

            ParallelDescriptor::ReduceLongSum(num_particles_read, IOProcNumber);

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "Replication the domain with vector           ";
                for (int d=0; d<BL_SPACEDIM; ++d) {
                    std::cout << lNrep[d] << " ";
                }
                std::cout << '\n';
                std::cout << "Total number of particles read in          : " << num_particles_read << '\n';
                std::cout << "Total number of particles after replication: " << num_particles      << '\n';
            }
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromAsciiFile() time: " << runtime << '\n';
    }
}

//
// The format of a binary particle init file:
//
// NP -- The number of particles in the file.   A "long".
// DM -- Our dimension.  Either 1, 2, or 3.     A "int".
// NX -- The amount of "extra" data.            A "int".
// NP*(DM+NX) native floating-point numbers.    A "float" or "double".
//
// Note that there is nothing separating all these values.
// They're packed into the binary file like sardines.
//

template <int N>
void
ParticleContainer<N>::InitFromBinaryFile (const std::string& file,
                                          int                extradata)
{
    BL_PROFILE("ParticleContainer<N>::InitFromBinaryFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= N);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // The number of MPI processes that read from the file.
    // We limit this to a rather small number since there's a limit
    // to the number of independent I/O channels on most filesystems.
    //
    const int NReaders = ParticleBase::MaxReaders();

    BL_ASSERT(NReaders <= NProcs);
    //
    // How many particles each NReaders reads before redistributing.
    //
    const long NPartPerRedist = ParticleBase::MaxParticlesPerRead();

    if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Reading with " << NReaders << " readers\n";
        std::cout << "Redistributing after every " << NPartPerRedist << " particles for each reader\n";
    }
    //
    // tmp_particles should mirror how m_particles is built.
    // At the end of this routine it'll be the new m_particles.
    //
    Array<PMap> tmp_particles;

    tmp_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
      m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    tmp_particles.resize(m_gdb->finestLevel()+1);
      m_particles.resize(m_gdb->finestLevel()+1);
    //
    // All the processors need to participate in Redistribute() though.
    //
    int  NX = 0;
    long NP = 0;

    std::ifstream ifs;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    //
    // The "set" of MPI processor numbers of the readers.  We want
    // to be able to easily check whether or not we're a reader.
    //
    std::set<int> readers;
    //
    // The same set but in ascending order.
    //
    Array<int> rprocs(NReaders);

    if (NReaders == NProcs)
    {
        //
        // Just set'm.
        //
        for (int i = 0; i < NProcs; i++)
            rprocs[i] = i;
    }
    else
    {
        //
        // The I/O Proc builds a set of NReader integers in the range: [0,NProcs-1].
        //
        // It then broadcast'm to all MPI procs.
        //
        // We want these to be as evenly distributed over the full set of
        // [0,NProcs-1] MPI processors as possible, so that when reading we
        // minimize the number of readers per Node, and hence can use more
        // of the available Node memory for reading.
        //
        if (ParallelDescriptor::IOProcessor())
        {
            do
            {
                int n = int(BoxLib::Random() * (NProcs-1));

                BL_ASSERT(n >= 0);
                BL_ASSERT(n < NProcs);

                readers.insert(n);
            }
            while (readers.size() < NReaders);

            BL_ASSERT(readers.size() == rprocs.size());

            int i = 0;

            for (std::set<int>::const_iterator it = readers.begin(), End = readers.end();
                 it != End;
                 ++it, ++i)
            {
                rprocs[i] = *it;
            }
        }

        ParallelDescriptor::Bcast(rprocs.dataPtr(), rprocs.size(), IOProc);
    }

    if (readers.empty())
    {
        //
        // Set readers for non I/O procs.
        //
        readers.insert(rprocs.begin(), rprocs.end());

        BL_ASSERT(readers.size() == rprocs.size());

        BL_ASSERT(rprocs.size() == NReaders);
    }

    int RealSizeInFile = 0;

    if (readers.find(MyProc) != readers.end())
    {
        int DM = 0;

        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(file.c_str(), std::ios::in|std::ios::binary);

        if (!ifs.good())
            BoxLib::FileOpenFailed(file);

        ifs.read((char*)&NP, sizeof(NP));
        ifs.read((char*)&DM, sizeof(DM));
        ifs.read((char*)&NX, sizeof(NX));
        //
        // NP MUST be positive!
        //
        if (NP <= 0)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): NP <= 0");
        //
        // DM must equal BL_SPACEDIM.
        //
        if (DM != BL_SPACEDIM)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): DM != BL_SPACEDIM");
        //
        // NX MUST be in [0,N].
        //
        if (NX < 0 || NX > N)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): NX < 0 || NX > N");
        //
        // Can't ask for more data than exists in the file!
        //
        if (extradata > NX)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): extradata > NX");
        //
        // Figure out whether we're dealing with floats or doubles.
        //
        // First get our current position.
        //
        const std::streamoff CURPOS = ifs.tellg();
        //
        // Seek to end of file.
        //
        ifs.seekg(0,std::ios::end);
        //
        // ENDPOS - CURPOS should bracket the particle data.
        //
        const std::streamoff ENDPOS = ifs.tellg();

        RealSizeInFile = (ENDPOS - CURPOS) / (NP*(DM+NX));

        BL_ASSERT(RealSizeInFile == sizeof(float) || RealSizeInFile == sizeof(double));
        //
        // Now set stream back to earlier position.
        //
        ifs.seekg(CURPOS, std::ios::beg);
        //
        // Skip to our place in the file.
        //
        int id = 0;
        for ( ; id < NReaders; id++)
            if (rprocs[id] == MyProc)
                break;

        BL_ASSERT(id >= 0 && id < NReaders);

        const std::streamoff NSKIP = id * (NP/NReaders) * (DM+NX) * RealSizeInFile;

        if (NSKIP > 0)
        {
            ifs.seekg(NSKIP, std::ios::cur);
        }

        if (!ifs.good())
        {
            std::string msg("ParticleContainer::InitFromBinaryFile(");
            msg += file;
            msg += ") failed @ 1";
            BoxLib::Error(msg.c_str());
        }
    }
    //
    // Everyone needs to know NP -- the number of particles in the file.
    //
    ParallelDescriptor::ReduceLongMax(NP);
    //
    // How many particles each reader gets to read.
    //
    long MyCnt = NP / NReaders;

    if (MyProc == rprocs[0])
        //
        // Give any remainder to the first reader.
        //
        MyCnt += NP % NReaders;

    const Geometry& geom = m_gdb->Geom(0);

    const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
    const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };

    long how_many_redists = NP / (NPartPerRedist*NReaders), how_many_read = 0;

    if (NP % (NPartPerRedist*NReaders)) how_many_redists++;

    Array<float>  fxtra, fignore;
    Array<double> dxtra, dignore;

    if (extradata > 0)
    {
        fxtra.resize(extradata);
        dxtra.resize(extradata);
    }

    if ((NX-extradata) > 0)
    {
        fignore.resize(NX-extradata);
        dignore.resize(NX-extradata);
    }

    for (int j = 0; j < how_many_redists; j++)
    {
        if (readers.find(MyProc) != readers.end())
        {
            ParticleType p;

            BL_ASSERT(MyCnt > how_many_read);

            const long NR = std::min((MyCnt-how_many_read), NPartPerRedist);

            for (long i = 0; i < NR; i++)
            {
                //
                // We don't read in m_id or m_cpu.  We'll set those later
                // in a manner to guarantee the global uniqueness of the pair.
                //
                if (RealSizeInFile == sizeof(float))
                {
                    float fpos[BL_SPACEDIM];

                    ifs.read((char*)&fpos[0], BL_SPACEDIM*sizeof(float));

                    D_TERM(p.m_pos[0] = fpos[0];,
                           p.m_pos[1] = fpos[1];,
                           p.m_pos[2] = fpos[2];);

                }
                else if (RealSizeInFile == sizeof(double))
                {
                    double dpos[BL_SPACEDIM];

                    ifs.read((char*)&dpos[0], BL_SPACEDIM*sizeof(double));

                    D_TERM(p.m_pos[0] = dpos[0];,
                           p.m_pos[1] = dpos[1];,
                           p.m_pos[2] = dpos[2];);
                }
                //
                // If particle is right on a domain boundary then move it just inside the boundary.
                //
                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    if (p.m_pos[d] <= ProbLo[d]) p.m_pos[d] += Delta[d];
                    if (p.m_pos[d] >= ProbHi[d]) p.m_pos[d] -= Delta[d];
                }
                //
                // Read in any "extradata".
                //
                if (extradata > 0)
                {
                    if (RealSizeInFile == sizeof(float))
                    {
                        ifs.read((char*)&fxtra[0], extradata*sizeof(float));

                        for (int i = 0; i < extradata; i++)
                            p.m_data[i] = fxtra[i];
                    }
                    else if (RealSizeInFile == sizeof(double))
                    {
                        ifs.read((char*)&dxtra[0], extradata*sizeof(double));

                        for (int i = 0; i < extradata; i++)
                            p.m_data[i] = dxtra[i];
                    }
                }
                //
                // Read any remaining data for this particle.
                //
                if ((NX-extradata) > 0)
                {
                    if (RealSizeInFile == sizeof(float))
                    {
                        ifs.read((char*)&fignore[0], (NX-extradata)*sizeof(float));
                    }
                    else if (RealSizeInFile == sizeof(double))
                    {
                        ifs.read((char*)&dignore[0], (NX-extradata)*sizeof(double));
                    }
                }

                if (!ifs.good())
                {
                    std::string msg("ParticleContainer::InitFromBinaryFile(");
                    msg += file;
                    msg += ") failed @ 2";
                    BoxLib::Error(msg.c_str());
                }

                if (!ParticleBase::Where(p,m_gdb))
                {
                    ParticleBase::PeriodicShift(p,m_gdb);

                    if (!ParticleBase::Where(p,m_gdb))
                    {
                        std::cout << "BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                        for (int d = 0; d < BL_SPACEDIM; d++)
                        {
                            std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                        }

                        BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): invalid particle");
                    }
                }

                p.m_id  = ParticleBase::NextID();
                p.m_cpu = MyProc;

                m_particles[p.m_lev][p.m_grid].push_back(p);
            }

            how_many_read += NR;
        }

        Redistribute(true);
        //
        // Move particles in m_particles into tmp_particles so that
        // we don't keep trying to redistribute particles that have
        // already been redistributed correctly.
        //
        for (int lev = 0; lev < m_particles.size(); lev++)
        {
            PMap& pmap     = m_particles[lev];
            PMap& tmp_pmap = tmp_particles[lev];

            for (typename PMap::iterator pmap_it = pmap.begin(), End = pmap.end();
                 pmap_it != End;
                 ++pmap_it)
            {
                PBox& pbox     = pmap_it->second;
                PBox& tmp_pbox = tmp_pmap[pmap_it->first];

                tmp_pbox.insert(tmp_pbox.end(), pbox.begin(), pbox.end());

                PBox().swap(pbox);
            }

            PMap().swap(pmap);
        }
    }
    //
    // Make tmp_particles the new m_particles.
    //
    tmp_particles.swap(m_particles);
    //
    // Add up all the particles read in to get the total number of particles.
    //
    if (m_verbose > 0)
    {
        long num_particles_read = how_many_read; 

        ParallelDescriptor::ReduceLongSum(num_particles_read, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "\nTotal number of particles: " << num_particles_read << '\n';
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryFile() time: " << runtime << '\n';
    }
}

//
// This function expects to read a file containing the pathnames of
// binary particles files needing to be read in for input.  It expects
// one file name per line.
//

template <int N>
void
ParticleContainer<N>::InitFromBinaryMetaFile (const std::string& metafile,
                                              int                extradata)
{
    BL_PROFILE("ParticleContainer<N>::InitFromBinaryMetaFile()");
    const Real strttime = ParallelDescriptor::second();

    std::ifstream ifs(metafile.c_str(), std::ios::in);

    std::string file;

    for (;;)
    {
        std::getline(ifs,file);

        if (!ifs.good()) break;

        if (m_verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryMetaFile: processing file: " << file << '\n';

        InitFromBinaryFile(file, extradata);
    }

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryMetaFile() time: " << runtime << '\n';
    }
}

template <int N>
void
ParticleContainer<N>::InitCosmo1ppcMultiLevel(MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                                              const Real particleMass, int disp_idx, int vel_idx, 
                                              BoxArray &baWhereNot, int lev)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Array<int> calls;

    calls.resize(m_gdb->initialBaLevels()+1);

    calls[lev]++;

    if (calls[lev] > 1) return;

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->initialBaLevels()+1);

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    Real         vel[BL_SPACEDIM];
    
    Real 	mean_disp[BL_SPACEDIM]={D_DECL(0,0,0)};


    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));
		    totalcount++;
		    if (baWhereNot.contains(indices)) 
		    {
                       continue;
		    }

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,disp_idx+n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
			//
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+Real(0.5))*dx[n];
			if(disp[n]*disp_fac[n]>dx[n]/2.0)
			  outcount[n]++;
			if(disp[n]*disp_fac[n]<-dx[n]/2.0)
			  outcountminus[n]++;
			mean_disp[n]+=fabs(disp[n]);
			//
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.m_pos[n] += disp[n] * disp_fac[n];

                        //
			// Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
	                p.m_data[n+1] = vel[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		    }
		    BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
		    //handle particles that ran out of this level into a finer one. 
		    if (baWhereNot.contains(p.m_cell))
		    {
		      outside_counter++;
		      ParticleType newp[8];
		      for (int i=0;i<8;i++)
		      {
			newp[i].m_data[0] = particleMass/8.0;
			newp[i].m_id = ParticleBase::NextID();
			newp[i].m_cpu     = MyProc;
			for (int dim=0;dim<BL_SPACEDIM;dim++)
			{
			  newp[i].m_pos[dim]=p.m_pos[dim]+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
			  newp[i].m_data[dim+1]=p.m_data[dim+1];
			  
			  
			}
		      
			if (!ParticleBase::Where(newp[i],m_gdb))
			{
			    ParticleBase::PeriodicShift(newp[i],m_gdb);
    
			    if (!ParticleBase::Where(newp[i],m_gdb))
				BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
			}
			m_particles[newp[i].m_lev][newp[i].m_grid].push_back(newp[i]);
		      }
		      
		    }
	            
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
		    else
		      m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

template <int N>
void
ParticleContainer<N>::InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppc()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(0);
    const Real*     dx       = geom.CellSize();

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    const Real   len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                             geom.ProbLength(1),
                                             geom.ProbLength(2)) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

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
                                      (indices[n]+Real(0.5))*dx[n] +
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
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		    }

	            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
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
    BL_PROFILE("ParticleContainer<N>::InitCosmo()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };
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
	    std::cout << '\n' << '\n';
	    std::cout << "Your particle lattice will have different spacings in the spatial directions!" << '\n';
	    std::cout << "You might want to change the particle number or the algorithm... ;)" << '\n';
	    std::cout << '\n' << '\n';
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
		               + (indices[n] + Real(0.5))*len[n]/n_part[n]
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

	                if (!ParticleBase::Where(p,m_gdb))
                        {
      		            ParticleBase::PeriodicShift(p,m_gdb);

                            if (!ParticleBase::Where(p,m_gdb))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
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
        std::cout << "Done with equidistant placement" << '\n';
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << '\n';
    }

    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    PMap& pmap = m_particles[0];
    //
    // Make sure, that mf and m_gdb.boxArray(0) are defined on the same boxarray.
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
            ParticleBase::GetGravity(dfab,m_gdb->Geom(p.m_lev),p,disp);

            D_TERM(p.m_pos[0] += len[0]*disp[0];,
                   p.m_pos[1] += len[1]*disp[1];,
                   p.m_pos[2] += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.m_data[1] = vel_fac[0]*disp[0];,
                   p.m_data[2] = vel_fac[1]*disp[1];,
                   p.m_data[3] = vel_fac[2]*disp[2];);

            if (!ParticleBase::Where(p,m_gdb))
            {
	        ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                    BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
	    }

            ParticleBase::Reset(p,m_gdb,true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with particle displacement" << '\n';
    }

    if (m_verbose > 1)
    {
        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "InitCosmo() done time: " << runtime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitRandom (long          icount,
                                  unsigned long iseed,
                                  Real          mass,
                                  bool          serialize,
                                  RealBox       containing_bx)
{
    BL_PROFILE("ParticleContainer<N>::InitRandom()");
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_gdb != 0);

    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       NProcs   = ParallelDescriptor::NProcs();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    Real r, x, len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };

    // We will enforce that the particles are within the containing_bx.  
    // If containing_bx is not passed in, it defaults to the full domain.
    if (!containing_bx.ok()) containing_bx = geom.ProbDomain();

    // containing_bx is assumed to lie within the domain.
    BL_ASSERT(geom.ProbDomain().contains(containing_bx));

    const Real* xlo = containing_bx.lo();
    const Real* xhi = containing_bx.hi();

    BoxLib::InitRandom(iseed+MyProc);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

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
        Array<ParticleBase::RealType> pos(icount*BL_SPACEDIM);

        if (ParallelDescriptor::IOProcessor())
        {
            for (long j = 0; j < icount; j++)
            {
                for (int i = 0; i < BL_SPACEDIM; i++)
                {
                    do
                    {
                        r = BoxLib::Random();
                        x = geom.ProbLo(i) + (r * len[i]);
                    }
                    while (x < xlo[i] || x > xhi[i]);

                    pos[j*BL_SPACEDIM + i] = x;
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
                p.m_pos[i] = pos[j*BL_SPACEDIM + i];
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

            if (!ParticleBase::Where(p,m_gdb))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());

            const int who = m_gdb->ParticleDistributionMap(p.m_lev)[p.m_grid];

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
                    x = geom.ProbLo(i) + (r * len[i]);
                }
                while (x < xlo[i] || x > xhi[i]);

                p.m_pos[i] = x;

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

            if (!ParticleBase::Where(p,m_gdb))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
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
    BL_PROFILE("ParticleContainer<N>::MoveRandom(lev)");
    BL_ASSERT(OK());
    BL_ASSERT(m_gdb != 0);
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
    const Real* dx                = m_gdb->Geom(lev).CellSize();
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

            ParticleBase::Reset(p,m_gdb,true);
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
    BL_PROFILE("ParticleContainer<N>::IncrementWithTotal(lev)");
    BL_ASSERT(OK());

    if (m_particles.empty()) return 0;

    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const PMap& pmap = m_particles[lev];
  
    long num_particles_in_domain = 0;

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf))
    {
        // If we are already working with the internal mf defined on the
        // particle_box_array, then we just work with this.
        mf_pointer = &mf;
    }
    else
    {
        // If mf is not defined on the particle_box_array, then we need
        // to make a temporary mf_pointer here and copy it into mf at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),mf.nComp(),mf.nGrow(),
				  m_gdb->ParticleDistributionMap(lev), Fab_allocate);
    }

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int   grid = pmap_it->first;
        const PBox& pbox = pmap_it->second;
        FArrayBox&  fab  = (*mf_pointer)[grid];

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

    // If mf is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf) 
    {
	mf.copy(*mf_pointer,0,0,mf.nComp());  
	delete mf_pointer;
    }

    ParallelDescriptor::ReduceLongSum(num_particles_in_domain);

    return num_particles_in_domain;
}

template <int N>
Real
ParticleContainer<N>::estTimestep (MultiFab& gv,
                                   int             lev,
                                   Real            cfl) const
{
    return estTimestep(gv,1.0,lev,cfl);
}

template <int N>
Real
ParticleContainer<N>::estTimestep (MultiFab& gv,
                                   Real            a,
                                   int             lev,
                                   Real            cfl) const
{
    BL_PROFILE("ParticleContainer<N>::estTimestep(lev)");
    Real            dt               = 1e50;
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);

    if (m_particles.size() == 0)
        return dt;

    const Real      strttime         = ParallelDescriptor::second();
    const Geometry& geom             = m_gdb->Geom(lev);
    const Real*     dx               = geom.CellSize();
    const Real      adx[BL_SPACEDIM] = { D_DECL(a*dx[0],a*dx[1],a*dx[2]) };
    const PMap&     pmap             = m_particles[lev];
    int             tnum             = 1;

#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif

    Array<Real> ldt(tnum,1e50);

    long num_particles_at_level = 0;

    MultiFab* gv_pointer;
    if (OnSameGrids(lev, gv))
    {
        gv_pointer = 0;
    }
    else 
    {
        gv_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), gv.nComp(), gv.nGrow(),
				  m_gdb->ParticleDistributionMap(lev), Fab_allocate);
        gv_pointer->copy(gv,0,0,BL_SPACEDIM);
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        const PBox&      pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = (gv_pointer) ? (*gv_pointer)[grid] : gv[grid];

        num_particles_at_level += n;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            const Real mag_vel_over_dx[BL_SPACEDIM] = { D_DECL(std::abs(p.m_data[1])/adx[0],
                                                               std::abs(p.m_data[2])/adx[1],
                                                               std::abs(p.m_data[3])/adx[2]) };

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
                dt_part = std::min( dt_part, 1/std::sqrt(mag_grav/dx[0]) );

            int tid = 0;

#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ldt[tid] = std::min(dt_part, ldt[tid]);
        }
    }

    if (gv_pointer) delete gv_pointer;

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
ParticleContainer<N>::sumParticleMass (int lev) const
{
    BL_PROFILE("ParticleContainer<N>::sumParticleMass(lev)");
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
ParticleContainer<N>::sumParticleMomentum (int   lev,
                                           Real* mom) const
{
    BL_PROFILE("ParticleContainer<N>::sumParticleMomentum()");
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
// This is the single-level version -- it takes either cell-centered or node-centered MF's
//
template <int N>
void
ParticleContainer<N>::AssignDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                int       lev,
                                                int       ncomp,
                                                int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer<N>::AssignDensitySingleLevel()");
    BL_ASSERT(N >= 1);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

#ifdef NEUTRINO_PARTICLES
    BL_ASSERT(m_csq > 0.);
#endif

    if (lev >= m_particles.size())
    {
        //
        // Don't do anything if there are no particles at this level.
        //
        return;
    }

    // Keep the same external interface to the applications, but if the
    if (mf_to_be_filled.is_nodal())
    {
        AssignNodalDensitySingleLevel(mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
    else
    {
        AssignCellDensitySingleLevel(mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
}

//
// This is the single-level version for cell-centered density
//
template <int N>
void
ParticleContainer<N>::AssignCellDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                    int       lev,
                                                    int       ncomp,
                                                    int       particle_lvl_offset) const
{
    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, mf_to_be_filled.nGrow(),
				  m_gdb->ParticleDistributionMap(lev), Fab_allocate);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
    //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
    //    adds the value from one grid's ghost cell to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       BoxLib::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
        BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");
    }

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    // TODO: implement tiling with OpenMP in this grid loop.
    Array<int>         pgrd(ngrids);
    Array<const PBox*> pbxs(ngrids);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
        FArrayBox&  fab = (*mf_pointer)[pgrd[j]];

        Array<Real>    fracs;
        Array<IntVect> cells;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(fracs,cells) shared(plo,dx,dx_particle,gm,fab,ncomp,pbx)
#endif
        for (typename PBox::const_iterator it = pbx.begin(); it < pbx.end(); ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0) {
	      continue;
	    }

            const int M = ParticleBase::CIC_Cells_Fracs(p, plo, dx, dx_particle, fracs, cells);
            //
            // If this is not fully periodic then we have to be careful that the
            // particle's support leaves the domain unless we specifically want to ignore
            // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
            // We test this by checking the low and high corners respectively.
            //
            if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                    BoxLib::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		}
	    }

            for (int i = 0; i < M; i++)
            {
                if ( ! fab.box().contains(cells[i])) {
		  continue;
		}

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
		  continue;
		}
                //
                // Sum up mass in first component.
                //
#ifdef NEUTRINO_PARTICLES
                if (m_relativistic)
                {
                    Real vsq = 0.0;
                    for (int n = 1; n < ncomp; n++) {
                       vsq += p.m_data[n] * p.m_data[n];
		    }
                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                }
                else 
#endif /* NEUTRINO_PARTICLES */
                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                }
                // 
                // Sum up momenta in next components.
                //
                for (int n = 1; n < ncomp; n++)
#ifdef _OPENMP
#pragma omp atomic
#endif
                   fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
            }
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());
    //
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    //
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

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
// This is the single-level version for nodal density
//
template <int N>
void
ParticleContainer<N>::AssignNodalDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                     int       lev,
                                                     int       ncomp,
                                                     int       particle_lvl_offset) const
{
    BoxLib::Error("AssignNodalDensitySingleLevel not implemented yet");
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
ParticleContainer<N>::AssignDensityAndVels (PArray<MultiFab>& mf,
                                            int               lev_min) const
{
    AssignDensity(mf, lev_min, BL_SPACEDIM+1);
}

template <int N>
void
ParticleContainer<N>::AssignDensity (PArray<MultiFab>& mf_to_be_filled, 
                                     int               lev_min,
                                     int               ncomp,
                                     int               finest_level) const
{
    BL_PROFILE("ParticleContainer<N>::AssignDensity()");
    BL_ASSERT(N >= 1);
    BL_ASSERT(N >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

#ifdef NEUTRINO_PARTICLES
    BL_ASSERT(m_csq > 0.0);
#endif

    if (finest_level == -1)
    {
        finest_level = m_gdb->finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level))
    {
        finest_level--;
    }
    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level, 
    // lev_index is the corresponding index for mf. 
    //
    PArray<MultiFab>* mf;
    PArray<MultiFab>  mf_part;

    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1-lev_min, PArrayManage);
    for (int lev = lev_min; lev <= finest_level; lev++)
    { 
        const int lev_index = lev - lev_min;
        mf_to_be_filled.set(lev_index, new MultiFab(m_gdb->boxArray(lev), ncomp, 1));

        for (MFIter mfi(mf_to_be_filled[lev_index]); mfi.isValid(); ++mfi) {
            mf_to_be_filled[lev_index][mfi].setVal(0);
	}
    }

    // Test whether the grid structure of the boxArray is the same
    //       as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++) {
        if (!OnSameGrids(lev, mf_to_be_filled[lev-lev_min])) {
	    all_grids_the_same = false;
	    break;
	}
    }

    if (all_grids_the_same)
    {
        mf = &mf_to_be_filled;
    } 
    else     
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1-lev_min, PArrayManage);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_part.set(lev_index, new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, 1,
						m_gdb->ParticleDistributionMap(lev)));
    
            for (MFIter mfi(mf_to_be_filled[lev_index]); mfi.isValid(); ++mfi) {
                mf_part[lev_index][mfi].setVal(0);
	    }
        }

        mf = &mf_part;
    }

    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignDensitySingleLevel((*mf)[0],0,ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0].copy((*mf)[0],0,0,ncomp);
	}
        return;
    }
    
    const bool sub_cycle = m_gdb->subCycle();
    //
    // This'll hold all the info I need for parallel.
    // // What I'll use: m_lev, m_grid, m_cell & m_data[0..ncomp-1].
    //
    // This is the "data" needed by other MPI procs.
    //
    PMap data;

    const Real stime = ParallelDescriptor::second();
    //
    // Minimum M required.
    //
    const int M = D_TERM(2,+2,+4);

    Array<int>     cgrid(M);
    Array<int>    cwhich(M),  fwhich(M);
    Array<Real>    fracs(M),  cfracs(M);
    Array<IntVect> cells(M),  ccells(M), cfshifts(M);

    ParticleType pb;
    //
    // I'm going to allocate these badboys here & pass'm into routines that use'm.
    // This should greatly cut down on memory allocation/deallocation.
    //
    Array<IntVect>                    pshifts(27);
    std::vector< std::pair<int,Box> > isects;
    Array<int>                        fgrid(M);
    Array<Real>                       ffracs(M);
    Array<IntVect>                    fcells;
    //
    // "fvalid" contains all the valid region of the MultiFab at this level, together
    // with any ghost cells lying outside the domain, that can be periodically shifted into the
    // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
    // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
    // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
    //
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const Geometry& gm        = m_gdb->Geom(lev);
        const Geometry& gm_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1) : gm;
        const Geometry& gm_coarse = (lev > 0) ? m_gdb->Geom(lev-1) : gm;
        const Box&      dm        = gm.Domain();
        const Real*     dx        = gm.CellSize();
        const Real*     plo       = gm.ProbLo();
        const Real*     dx_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1).CellSize() : dx;
        const Real*     dx_coarse = (lev > 0) ? m_gdb->Geom(lev-1).CellSize() : dx;
        const int       lev_index = lev - lev_min;
        const BoxArray& grids     = (*mf)[lev_index].boxArray();
        const int       dgrow     = (lev == 0) ? 1 : m_gdb->MaxRefRatio(lev-1);

        BoxArray compfvalid, compfvalid_grown, fvalid = (*mf)[lev_index].boxArray();
        //
        // Do we have Fine->Crse overlap on a periodic boundary?
        // We want to add all ghost cells that can be shifted into valid region.
        //
        BoxList valid;

        for (int i = 0; i < grids.size(); i++)
        {
            if (gm.isAnyPeriodic())
            {
                const Box& dest = BoxLib::grow(grids[i],dgrow);

                if ( ! dm.contains(dest))
                {
                    for (int j = 0; j < grids.size(); j++)
                    {
                        BL_ASSERT(dm.contains(grids[j]));

                        gm.periodicShift(dest, grids[j], pshifts);

                        for (int k = 0; k < pshifts.size(); k++)
                        {
                            const Box& sbx = grids[j] + pshifts[k];
                            const Box& dbx = dest & sbx;

                            BL_ASSERT(dbx.ok());

                            valid.push_back(dbx);
                        }
                    }
                }
            }
        }
        if (valid.isNotEmpty())
        {
            //
            // We've got some Fine->Crse periodic overlap.
            // Don't forget to add the valid boxes too.
            //
            for (int i = 0; i < grids.size(); i++) {
                valid.push_back(grids[i]);
	    }
            fvalid = BoxArray(valid);
            fvalid.removeOverlap();
        }
        //
        // If we're at a lev < finestLevel, this is the coarsened fine BoxArray.
        // We use this for figuring out Crse->Fine issues.
        //
        BoxArray ccba;
        if (lev > 0)
        {
            ccba = m_gdb->boxArray(lev);
            ccba.coarsen(m_gdb->refRatio(lev-1));
        }
        BoxArray cfba;
        if (lev < finest_level)
        {
            cfba = m_gdb->boxArray(lev+1);
            cfba.coarsen(m_gdb->refRatio(lev));

            BL_ASSERT((*mf)[lev_index].boxArray().contains(cfba));
        }
        //
        // This is cfba with any shifted ghost cells.
        //
        BoxArray cfvalid = cfba;

        if (lev < finest_level)
        {
            BoxList cvalid;

            const BoxArray& cgrids = (*mf)[lev_index].boxArray();

            for (int i = 0; i < cfba.size(); i++)
            {
                if (gm.isAnyPeriodic())
                {
                    const Box& dest = BoxLib::grow(cfba[i],(*mf)[lev_index].nGrow());

                    if ( ! dm.contains(dest))
                    {
                        for (int j = 0; j < cgrids.size(); j++)
                        {
                            BL_ASSERT(dm.contains(cgrids[j]));

                            gm.periodicShift(dest, cgrids[j], pshifts);

                            for (int k = 0; k < pshifts.size(); k++)
                            {
                                const Box& sbx = cfba[i] - pshifts[k];

                                cvalid.push_back(sbx);
                            }
                        }
                    }
                }
            }
            if (cvalid.isNotEmpty())
            {
                //
                // We've got some Fine->Crse periodic overlap.
                // Don't forget to add the valid boxes too.
                //
                for (int i = 0; i < cfba.size(); i++) {
                    cvalid.push_back(cfba[i]);
		}
                cfvalid = BoxArray(cvalid);
                cfvalid.removeOverlap();
            }
        }
        //
        // The "+1" is so we enclose the valid region together with any
        //  ghost cells that can be periodically shifted into valid.
        //
        compfvalid = BoxLib::complementIn(BoxLib::grow(dm,dgrow+1), fvalid);

        compfvalid_grown = compfvalid;
        compfvalid_grown.grow(1);
        compfvalid_grown.removeOverlap();
            
        if (gm.isAnyPeriodic() && ! gm.isAllPeriodic())
        {
            BoxLib::Error("AssignDensity: problem must be periodic in no or all directions");
        }
        //
        // If we're at a lev > 0, this is the coarsened BoxArray.
        // We use this for figuring out Fine->Crse issues.
        //
        BoxArray cba;
        if (lev > 0)
        {
            cba = m_gdb->boxArray(lev);
            cba.coarsen(m_gdb->refRatio(lev-1));
        }
        //
        // Do the grids at this level cover the full domain? If they do
        // there can be no Fine->Crse interactions at this level.
        //
        const bool GridsCoverDomain = fvalid.contains(m_gdb->Geom(lev).Domain());
        
        for (typename PMap::const_iterator pmap_it = m_particles[lev].begin(),
                 pmapEnd = m_particles[lev].end();
             pmap_it != pmapEnd;
             ++pmap_it)
        {
            const PBox& pbx = pmap_it->second;
            FArrayBox&  fab = (*mf)[lev_index][pmap_it->first];

            for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end();
                 it != End;
                 ++it)
            {
                const ParticleType& p = *it;

                if (p.m_id <= 0) {
		  continue;
		}
                //
                // Get "fracs" and "cells" for the particle "p" at this level.
                //
                const int M = ParticleBase::CIC_Cells_Fracs(p, plo, dx, fracs, cells);
                //
                // If this is not fully periodic then we have to be careful that no
                // particle's support leaves the domain. We test this by checking the low
                // and high corners respectively.
                //
                if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                    if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                        BoxLib::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		    }
		}
                //
                // This section differs based on whether we subcycle.
                // Without subcycling we use the "stretchy" support for particles.
                // With subcycling a particles support is strictly defined 
                // by its resident level.
                //
                if (sub_cycle)
                {
                    bool isFiner    = false;
                    bool isBoundary = false;
                    //
                    // First sum the mass in the valid region
                    //
                    for (int i = 0; i < M; i++)
                    {
                        if (cfvalid.contains(cells[i]))
                        {
                            //
                            // Some part of the particle's mass lies in a 
                            // finer region; we'll deal with it shortly.
                            //
                            isFiner    = true;
                            isBoundary = true;
                            continue;
                        }
                        if ( ! fvalid.contains(cells[i]))
                        {
                            //
                            // We're out of the valid region.
                            //
                            isBoundary = true;
                            continue;
                        }
                        //
                        // Sum up mass in first component.
                        //
#ifdef NEUTRINO_PARTICLES
                        if (m_relativistic)
                        {
                            Real vsq = 0.0;
                            for (int n = 1; n < ncomp; n++) {
                               vsq += p.m_data[n] * p.m_data[n];
			    }
                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                            fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                        }
                        else 
#endif
                        {
                            fab(cells[i],0) += p.m_data[0] * fracs[i];
                        }
                        //
                        // Sum up momenta in next components.
                        //

                        // If the domain is not periodic and we want to let particles
                        //    live near the boundary but "throw away" the contribution that 
                        //    does not fall into the domain ...
                        if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
			     ! gm.Domain().contains(cells[i]))
			{
			  continue;
			}

                        for (int n = 1; n < ncomp; n++) {
                            fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
			}
                    }
                    //
                    // Deal with mass that doesn't belong at this level.
                    // Here we assume proper nesting so that only one special case can
                    // be true for a given particle.
                    //
                    if (isBoundary)
                    {
                        if (isFiner)
                        {
                            BL_ASSERT(lev < finest_level);
                            //
                            // We're at a coarse->fine interface
                            //
                            // get fine cells/fracs
                            //
                            const int MF = ParticleBase::CIC_Cells_Fracs(p, plo, dx_fine ,dx, ffracs, fcells);

                            for (int j = 0; j < MF; j++)
                            {
                                //
                                // Make sure this fine cell is valid. Check for periodicity.
                                //
                                const Box bx(fcells[j],fcells[j]);
                                gm_fine.periodicShift(bx, gm_fine.Domain(), pshifts);
                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(pshifts.size() == 1);
                                    fcells[j] = fcells[j] - pshifts[0];
                                }
                                (*mf)[lev_index + 1].boxArray().intersections(Box(fcells[j],fcells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first; 
                                const int who  = (*mf)[lev_index+1].DistributionMap()[grid];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf)[lev_index+1][grid](fcells[j],0) += p.m_data[0] * ffracs[j] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        (*mf)[lev_index+1][grid](fcells[j],0) += p.m_data[0] * ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf)[lev_index+1][grid](fcells[j],n) += p.m_data[n] * p.m_data[0] * ffracs[j];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev+1;
                                    pb.m_grid = grid;
                                    pb.m_cell = fcells[j];
                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] *  ffracs[j] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        pb.m_data[0] = p.m_data[0] *  ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * ffracs[j];
				    }
                                    data[who].push_back(pb);
                                }
                            }
                        }
                        else if (lev_index > 0)
                        {
                            //
                            // We must be at a fine->coarse interface.
                            //
                            const int MC = ParticleBase::CIC_Cells_Fracs(p, plo, dx_coarse, dx, cfracs, ccells);
                            for (int j = 0; j < MC; j++)
                            {
                                //
                                // Make sure this coarse cell isn't in this level's valid region.
                                // This may not matter.
                                //
                                if (cba.contains(ccells[j]))
                                    continue;
                                //
                                // Check for periodicity.
                                //
                                const Box bx(ccells[j],ccells[j]);
                                gm_coarse.periodicShift(bx, gm_coarse.Domain(), pshifts);

                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(pshifts.size() == 1);
                                    ccells[j] = ccells[j] - pshifts[0]; 
                                }
                                //
                                // Find its resident grid.
                                //
                                (*mf)[lev_index - 1].boxArray().intersections(Box(ccells[j],ccells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first;
                                const int who  = (*mf)[lev_index-1].DistributionMap()[grid];
                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf)[lev_index-1][grid](ccells[j],0) += p.m_data[0] * cfracs[j] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        (*mf)[lev_index-1][grid](ccells[j],0) += p.m_data[0] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf)[lev_index-1][grid](ccells[j],n) += p.m_data[n] * p.m_data[0] * cfracs[j];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev-1;
                                    pb.m_grid = grid;
                                    pb.m_cell = ccells[j];
                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] * cfracs[j] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        pb.m_data[0] = p.m_data[0] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * cfracs[j];
				    }

                                    data[who].push_back(pb);
                                }
                            }
                        }
                        else
                        {
                            // The mass is below levels we care about. Ignore it.
                        }
                    }
                }
                else 
                {
                    bool AnyCrseToFine = false;
                    if (lev < finest_level) {
                        AnyCrseToFine = ParticleBase::CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
		    }
                    //
                    // lev_index > 0 means that we don't do F->C for lower levels
                    // This may mean that the mass fraction is off.
                    //
                    bool AnyFineToCrse = false;
                    if (lev_index > 0 && !GridsCoverDomain)
                        AnyFineToCrse = ParticleBase::FineToCrse(p,lev,m_gdb,cells,fvalid,compfvalid_grown,ccells,cfracs,fwhich,cgrid,pshifts,isects);

                    BL_ASSERT(!(AnyCrseToFine && AnyFineToCrse));

                    if ( ! AnyCrseToFine && ! AnyFineToCrse)
                    {
                        //
                        // By far the most common case.  Just do it!
                        //
                        for (int i = 0; i < M; i++)
                        {

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if (! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i]))
			    {
			      continue;
			    }
                            //
                            // Sum up mass in first component.
                            //
#ifdef NEUTRINO_PARTICLES
                            if (m_relativistic)
                            {
                                Real vsq = 0.0;
                                for (int n = 1; n < ncomp; n++) {
                                   vsq += p.m_data[n] * p.m_data[n];
				}
                                Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                            }
                            else 
#endif
                            {
                                fab(cells[i],0) += p.m_data[0] * fracs[i];
                            }
                            //
                            // Sum up momenta in next components.
                            //
                            for (int n = 1; n < ncomp; n++) {
                                fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
			    }
                        }
                    }
                    else if (AnyFineToCrse)
                    {
                        Real sum_crse = 0, sum_fine = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // We're at a Fine->Crse boundary.
                                //
                                BL_ASSERT(cgrid[i] >= 0);
                                BL_ASSERT(cgrid[i] < (*mf)[lev_index-1].size());
                                //
                                // Here we need to update the crse region.  The coarse
                                // region is always going to be updated if we have a
                                // particle in a cell bordering a Fine->Crse boundary.
                                //
                                const int who = (*mf)[lev_index-1].DistributionMap()[cgrid[i]];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    if ( ! (*mf)[lev_index-1][cgrid[i]].box().contains(ccells[i])) {
				      continue;
				    }

                                    // If the domain is not periodic and we want to let particles
                                    //    live near the boundary but "throw away" the contribution that 
                                    //    does not fall into the domain ...
                                    if (! gm_coarse.isAllPeriodic() && allow_particles_near_boundary &&
				        ! gm_coarse.Domain().contains(ccells[i]))
				    {
				      continue;
				    }

                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf)[lev_index-1][cgrid[i]](ccells[i],0) += p.m_data[0] * cfracs[i] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        (*mf)[lev_index-1][cgrid[i]](ccells[i],0) += p.m_data[0] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf)[lev_index-1][cgrid[i]](ccells[i],n) += p.m_data[n] * p.m_data[0] * cfracs[i];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev-1;
                                    pb.m_grid = cgrid[i];
                                    pb.m_cell = ccells[i];
                                    //
                                    // Sum up mass in first component.
                                    //
#ifdef NEUTRINO_PARTICLES
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] * cfracs[i] * gamma;
                                    }
                                    else 
#endif
                                    {
                                        pb.m_data[0] = p.m_data[0] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * cfracs[i];
				    }
                                    data[who].push_back(pb);
                                }

                                sum_crse += cfracs[i];
                            }
                        }
                        //
                        // We've updated the Crse cells.  Now we have to update the fine
                        // cells in such a way that the total amount of mass we move
                        // around is precisely p.m_data[0]. In other words, the fractions
                        // we use at crse and fine have to sum to zero.  In the fine
                        // case, we have to account for the case where one or more of the
                        // cell indices is not in the valid region of the box containing 
                        // the particle.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {
                            //
                            // Reusing "fwhich" to indicate fine cells that need massaging.
                            //
                            fwhich[i] = true;

                            if ( ! compfvalid_grown.contains(cells[i]))
                            {
                                //
                                // Go ahead and add the full correct amount to these cells.
                                // They can't touch a Fine->Crse boundary.
                                //
                                sum_fine += fracs[i];
                                //
                                // Sum up mass in first component.
                                //
#ifdef NEUTRINO_PARTICLES
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                                }
                                else 
#endif
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
				}
                                fwhich[i] = false;
                            }
                            else if (compfvalid.contains(cells[i]))
                            {
                                fwhich[i] = false;
                            }
                        }

                        const Real sum_so_far = sum_crse + sum_fine; 

                        BL_ASSERT(sum_so_far > 0);
                        BL_ASSERT(sum_so_far < 1);

                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {       
                            if (fwhich[i])
                                //
                                // Got to weight cells in this direction differently.
                                //
                                sum_fine += fracs[i];
                        }

                        const Real mult = (1 - sum_so_far) / sum_fine;
                        //
                        // Now add the weighted amount to the fine cells touching the c-f interface.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // Sum up mass in first component.
                                //
#ifdef NEUTRINO_PARTICLES
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * mult * gamma;
                                }
                                else 
#endif
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * mult;
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i] * mult;
				}

                                sum_fine += fracs[i] * mult;
                            }
                        }

                        BL_ASSERT(std::abs(1-(sum_fine+sum_so_far)) < 1.e-9);
                    }
                    else if (AnyCrseToFine)
                    {
                        Real sum = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (!cwhich[i])
                            {
                                // If the domain is not periodic and we want to let particles
                                //    live near the boundary but "throw away" the contribution that 
                                //    does not fall into the domain ...
                                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
				     ! gm.Domain().contains(ccells[i]))
				{
				  continue;
				}
                                //
                                // Sum up mass in first component.
                                //
#ifdef NEUTRINO_PARTICLES
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                                }
                                else 
#endif
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
				}

                                sum += fracs[i];
                            }
                            else
                            {
                                //
                                // We're at a Crse->Fine boundary.
                                //
                                ParticleBase::FineCellsToUpdateFromCrse(p,lev,m_gdb,cells[i],cfshifts[i],fgrid,ffracs,fcells,isects);

                                for (int j = 0; j < fcells.size(); j++)
                                {
                                    const int who = (*mf)[lev_index+1].DistributionMap()[fgrid[j]];

                                    if (who == ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
#ifdef NEUTRINO_PARTICLES
                                        if (m_relativistic)
                                        {
                                            Real vsq = 0.0;
                                            for (int n = 1; n < ncomp; n++) {
                                               vsq += p.m_data[n] * p.m_data[n];
					    }
                                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                            (*mf)[lev_index+1][fgrid[j]](fcells[j],0) += p.m_data[0] * fracs[i] * ffracs[j] * gamma;
                                        }
                                        else 
#endif
                                        {
                                            (*mf)[lev_index+1][fgrid[j]](fcells[j],0) += p.m_data[0] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            (*mf)[lev_index+1][fgrid[j]](fcells[j],n) += p.m_data[n] * p.m_data[0] * fracs[i] * ffracs[j];
					}
                                    }
                                    else
                                    {
                                        pb.m_lev  = lev+1;
                                        pb.m_grid = fgrid[j];
                                        pb.m_cell = fcells[j];
                                        //
                                        // Sum up mass in first component.
                                        //
#ifdef NEUTRINO_PARTICLES
                                        if (m_relativistic)
                                        {
                                            Real vsq = 0.0;
                                            for (int n = 1; n < ncomp; n++) {
                                               vsq += p.m_data[n] * p.m_data[n];
					    }
                                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                            pb.m_data[0] = p.m_data[0] * fracs[i] * ffracs[j] * gamma;
                                        }
                                        else 
#endif
                                        {
                                            pb.m_data[0] = p.m_data[0] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            pb.m_data[n] = p.m_data[n] * p.m_data[0] * fracs[i] * ffracs[j];
					}

                                        data[who].push_back(pb);
                                    }

                                    sum += fracs[i] * ffracs[j];
                                }
                            }
                        }

                        BL_ASSERT(std::abs(1-sum) < 1.e-9);
                    }
                }
            }
        }
    }
    //
    // Send any needed data to other MPI processes.
    // This "may" touch ghost cells so we want to do it before
    // the SumBoundary() stuff.
    //
    AssignDensityDoit(mf,data,ncomp,lev_min);

    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const int       lev_index = lev - lev_min;
        const Geometry& gm        = m_gdb->Geom(lev);
        const Real*     dx        = gm.CellSize();
        const Real      vol       = D_TERM(dx[0], *dx[1], *dx[2]);

        (*mf)[lev_index].SumBoundary(gm.periodicity());
        //
        // If ncomp > 1, first divide the momenta (component n) 
        // by the mass (component 0) in order to get velocities.
        // Be careful not to divide by zero.
        //
        for (int n = 1; n < ncomp; n++)
        {
            for (MFIter mfi((*mf)[lev_index]); mfi.isValid(); ++mfi)
            {
                (*mf)[lev_index][mfi].protected_divide((*mf)[lev_index][mfi],0,n,1);
            }
        }
        //
        // Only multiply the first component by (1/vol) because this converts mass
        // to density. If there are additional components (like velocity), we don't
        // want to divide those by volume.
        //
        (*mf)[lev_index].mult(1/vol,0,1);
    }

    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level,  
    // lev_index is the corresponding index for mf. 
    //
    // I believe that we don't need any information in ghost cells so we don't copy those.
    //
    if ( ! all_grids_the_same)
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_to_be_filled[lev_index].copy(mf_part[lev_index],0,0,1);
        }
    
    if (m_verbose > 1)
    {
        Real etime = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(etime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensity(multi-level) time: " << etime << '\n';
        }
    }
}

//
// Used by AssignDensity (PArray<MultiFab>& mf).
//
// Passes data needed by Crse->Fine or Fine->Crse to CPU that needs it.
//
// We store the data that needs to be sent in "data". Note that m_lev is the
// real particle level, while mf may start at a fine level (e.g. lvls 1 and 2).
// Consequently, we must subtract lev_min from m_lev to get the mf lev.
//
// We only use: m_lev, m_grid, m_cell & m_data[0..ncomp-1] from the particles.
//
//

template <int N>
void
ParticleContainer<N>::AssignDensityDoit (PArray<MultiFab>* mf,
                                         PMap&             data,
                                         int               ncomp,
                                         int               lev_min) const
{
    BL_PROFILE("ParticleContainer<N>::AssignDensityDoit()");
    BL_ASSERT(N >= ncomp);

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
        BL_ASSERT(data.empty());
        return;
    }

#if BL_USE_MPI
    //
    // We may have data that needs to be sent to another CPU.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

    int NumSnds = 0, NumRcvs = 0;

    for (typename PMap::const_iterator it = data.begin(), End = data.end(); it != End; ++it)
    {
        NumSnds        += it->second.size();
        Snds[it->first] = it->second.size();
    }

    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0) {
        //
        // There's no parallel work to do.
        //
        return;
    }

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++) {
        if (Snds[i] > 0) {
            SndCnts[i] = Snds[i];
	}
    }

    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
        }
    }
    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);
    //
    // The data we want to receive.
    //
    // We only use: m_lev, m_grid, m_cell & m_data[0..ncomp-1] from the particles.
    //
    const int iChunkSize = 2 + BL_SPACEDIM;
    const int rChunkSize = ncomp;

    Array<int>                    irecvdata (NumRcvs*iChunkSize);
    Array<ParticleBase::RealType> rrecvdata (NumRcvs*rChunkSize);

    Array<int>         index(2*RcvCnts.size());
    Array<MPI_Status>  stats(2*RcvCnts.size());
    Array<MPI_Request> rreqs(2*RcvCnts.size());

    const int SeqNumI = ParallelDescriptor::SeqNum();
    const int SeqNumR = ParallelDescriptor::SeqNum();
    //
    // Post the receives.
    //
    int idx = 0;
    for (IntIntMap::const_iterator it = RcvCnts.begin(), End = RcvCnts.end(); it != End; ++it, ++idx)
    {
        const int Who  = it->first;
        const int iCnt = it->second   * iChunkSize;
        const int rCnt = it->second   * rChunkSize;
        const int iIdx = rOffset[Who] * iChunkSize;
        const int rIdx = rOffset[Who] * rChunkSize;

        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        rreqs[2*idx+0] = ParallelDescriptor::Arecv(&irecvdata[iIdx],iCnt,Who,SeqNumI).req();
        rreqs[2*idx+1] = ParallelDescriptor::Arecv(&rrecvdata[rIdx],rCnt,Who,SeqNumR).req();
    }
    //
    // Send the data.
    //
    Array<int>                    isenddata;
    Array<ParticleBase::RealType> rsenddata;

    for (IntIntMap::const_iterator it = SndCnts.begin(), End = SndCnts.end(); it != End; ++it)
    {
        const int Who  = it->first;
        const int iCnt = it->second * iChunkSize;
        const int rCnt = it->second * rChunkSize;

        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        isenddata.resize(iCnt);
        rsenddata.resize(rCnt);

        PBox& pbox = data[Who];

        int ioff = 0, roff = 0;
        for (typename PBox::const_iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
        {
            const ParticleType& p = *it;

            isenddata[ioff+0] = p.m_lev  - lev_min;
            isenddata[ioff+1] = p.m_grid;

            D_TERM(isenddata[ioff+2] = p.m_cell[0];,
                   isenddata[ioff+3] = p.m_cell[1];,
                   isenddata[ioff+4] = p.m_cell[2];);

            ioff += iChunkSize;

            for (int n = 0; n < ncomp; n++) {
                rsenddata[roff+n] = p.m_data[n];
	    }

            roff += ncomp;
        }

        PBox().swap(pbox);

        ParallelDescriptor::Send(isenddata.dataPtr(),iCnt,Who,SeqNumI);
        ParallelDescriptor::Send(rsenddata.dataPtr(),rCnt,Who,SeqNumR);
    }
    //
    // Receive the data.
    //
    for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
    {
        ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
    }
    //
    // Now update "mf".
    //
    if (NumRcvs > 0)
    {
        const int*                    idata = irecvdata.dataPtr();
        const ParticleBase::RealType* rdata = rrecvdata.dataPtr();

        for (int i = 0; i < NumRcvs; i++)
        {
            const int     lev  = idata[0];
            const int     grd  = idata[1];
            const IntVect cell = IntVect(D_DECL(idata[2],idata[3],idata[4]));

            BL_ASSERT((*mf)[lev].DistributionMap()[grd] == MyProc);
            BL_ASSERT((*mf)[lev][grd].box().contains(cell));

            for (int n = 0; n < ncomp; n++) {
                (*mf)[lev][grd](cell,n) += rdata[n];
	    }

            idata += iChunkSize;
            rdata += rChunkSize;
        }
    }

#endif /*BL_USE_MPI*/
}

template <int N>
void
ParticleContainer<N>::GetParticleIDs (Array<int>& part_ids)
{
    BL_PROFILE("ParticleContainer<N>::GetParticleIDs()");
    //
    // This gives us the starting point into the part_ids array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++) {
        cnt += cnts[iproc];
    }

    std::cout << "PROC CNT " << ParallelDescriptor::MyProc() << " " << cnt << std::endl;

    // Each particle takes up 1 int so no need to multiply cnt by anything
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_ids.resize(npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&     pbx = pmap_it->second;
            const int n    = pbx.size();
    
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbx[i];
                if (p.m_id > 0)
                {
                    // Load the ID
                    part_ids[cnt] = p.m_id;

                    // Update counter
                    cnt++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(part_ids.dataPtr(),part_ids.size()); 
}

template <int N>
void
ParticleContainer<N>::GetParticleCPU (Array<int>& part_cpu)
{
    BL_PROFILE("ParticleContainer<N>::GetParticleCPU()");
    //
    // This gives us the starting point into the part_cpu array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up 1 int so no need to multiply cnt by anything
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_cpu.resize(npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&     pbx = pmap_it->second;
            const int n    = pbx.size();
    
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbx[i];
                if (p.m_id > 0)
                {
                    // Load the ID
                    part_cpu[cnt] = p.m_cpu;

                    // Update counter
                    cnt++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(part_cpu.dataPtr(),part_cpu.size()); 
}

template <int N>
void
ParticleContainer<N>::GetParticleLocations (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<N>::GetParticleLocations()");
    //
    // This gives us the starting point into the part_data array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up BL_SPACEDIM Reals
    cnt *= (BL_SPACEDIM);
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_data.resize(BL_SPACEDIM*npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&     pbx = pmap_it->second;
            const int n    = pbx.size();
    
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbx[i];
                if (p.m_id > 0)
                {
                    // Load positions
                    for (int d=0; d < BL_SPACEDIM; d++)
                        part_data[cnt+d] = p.m_pos[d];

                    // Update counter
                    cnt += BL_SPACEDIM;
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int N>
void
ParticleContainer<N>::GetParticleVelocities (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<N>::GetParticleVelocities()");
    // This assumes that the mass/charge is stored in the first position 
    //      in the particle data, followed by the velocity components
    int start_comp = 1;
    int   num_comp = BL_SPACEDIM;
    GetParticleData(part_data,1,BL_SPACEDIM);
}

template <int N>
void
ParticleContainer<N>::GetParticleData (Array<Real>& part_data, int start_comp, int num_comp)
{
    BL_PROFILE("ParticleContainer<N>::GetParticleData()");
    //
    // This gives us the starting point into the part_data array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

    //
    // Make sure we don't try to get more than we have
    //
    if (start_comp + num_comp > N)
        BoxLib::Error("Tried to grab too many components in GetParticleData!!");

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up num_comp Reals
    cnt*= num_comp;
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    part_data.resize(num_comp*npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&     pbx = pmap_it->second;
            const int n    = pbx.size();
    
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbx[i];
                if (p.m_id > 0)
                {
                    // Load particle data, whatever it is.
                    for (int d = 0; d < num_comp; d++)
                      part_data[cnt+d] = p.m_data[start_comp+d];

                    // Update counter
                    cnt += num_comp;
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int N>
void
ParticleContainer<N>::SetAllowParticlesNearBoundary (bool value)
{
    allow_particles_near_boundary = value; 
}

template <int N>
void
ParticleContainer<N>::SetParticleLocations (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<N>::SetParticleLocations()");
   // This gives us the starting point into the part_data array
   // If only one processor (or no MPI), then that's all we need
   int cnt = 0;

#if BL_USE_MPI
   Array<long> cnts(ParallelDescriptor::NProcs());

   // This returns the number of particles on this processor
   long lcnt = TotalNumberOfParticles(true,true);

   // This accumulates the "lcnt" values into "cnts"
   MPI_Gather(&lcnt,1,              
              ParallelDescriptor::Mpi_typemap<long>::type(),
              cnts.dataPtr(),
              1,
              ParallelDescriptor::Mpi_typemap<long>::type(),
              ParallelDescriptor::IOProcessorNumber(),
              ParallelDescriptor::Communicator());

   ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

   for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
       cnt += cnts[iproc];

   // Each particle takes up BL_SPACEDIM Reals
   cnt*= BL_SPACEDIM;
#endif

   // This is the total number of particles on *all* processors
   long npart = TotalNumberOfParticles(true,false);

   // Mass + locations
   if (part_data.size() != npart*BL_SPACEDIM)
       BoxLib::Abort("Sending in wrong size part_data to SetParticleLocations");

   for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
   {
       PMap& pmap = m_particles[lev];

       for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
       {
           PBox&     pbx = pmap_it->second;
           const int n    = pbx.size();
    
           for (int i = 0; i < n; i++)
           {
              ParticleType& p = pbx[i];
              if (p.m_id > 0)
              {
                  // Load positions
                  for (int d=0; d < BL_SPACEDIM; d++)
                     p.m_pos[d] = part_data[cnt+d];

                  // Update counter
                  cnt += BL_SPACEDIM;
              }
           }
       }
    }
}

template <int N>
void
ParticleContainer<N>::SetParticleVelocities (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<N>::SetParticleVelocities()");
   // This gives us the starting point into the part_data array
   // If only one processor (or no MPI), then that's all we need
   int cnt = 0;

#if BL_USE_MPI
   Array<long> cnts(ParallelDescriptor::NProcs());

   // This returns the number of particles on this processor
   long lcnt = TotalNumberOfParticles(true,true);

   // This accumulates the "lcnt" values into "cnts"
   MPI_Gather(&lcnt,1,              
              ParallelDescriptor::Mpi_typemap<long>::type(),
              cnts.dataPtr(),
              1,
              ParallelDescriptor::Mpi_typemap<long>::type(),
              ParallelDescriptor::IOProcessorNumber(),
              ParallelDescriptor::Communicator());

   ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

   for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
       cnt += cnts[iproc];

   // Each particle takes up (BL_SPACEDIM) Reals
   cnt*= (BL_SPACEDIM);
#endif

   // This is the total number of particles on *all* processors
   long npart = TotalNumberOfParticles(true,false);

   // Velocities
   if (part_data.size() != npart*(BL_SPACEDIM))
       BoxLib::Abort("Sending in wrong size part_data to SetParticleVelocities");

   for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
   {
       PMap& pmap = m_particles[lev];

       for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
       {
           PBox&     pbx = pmap_it->second;
           const int n    = pbx.size();
    
           for (int i = 0; i < n; i++)
           {
              ParticleType& p = pbx[i];
              if (p.m_id > 0)
              {
                  // Load velocities
                  for (int d=0; d < BL_SPACEDIM; d++)
                     p.m_data[d+1] = part_data[cnt+d];

                  // Update counter
                  cnt += BL_SPACEDIM;
              }
           }
       }
    }
}

template <int N>
void
ParticleContainer<N>::MultiplyParticleMass (int lev, Real mult)
{
    BL_PROFILE("ParticleContainer<N>::MultiplyParticleMass()");
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
    BL_PROFILE("ParticleContainer<N>::movePredict()");
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime = ParallelDescriptor::second();

    PMap& pmap = m_particles[lev];

    MultiFab* gv_pointer;
    if (OnSameGrids(lev,gv))
    {
        gv_pointer = 0;
    }
    else 
    {
        gv_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), 1, gv.nGrow(),
				  m_gdb->ParticleDistributionMap(lev), Fab_allocate);
        gv_pointer->copy(gv,0,0,BL_SPACEDIM);
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = (gv_pointer) ? (*gv_pointer)[grid] : gv[grid];

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

                ParticleBase::Reset(p,m_gdb,true);
            } 
        }
    }

    if (gv_pointer) delete gv_pointer;

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
    BL_PROFILE("ParticleContainer<N>::moveCorrect()");
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime = ParallelDescriptor::second();

    PMap& pmap = m_particles[lev];

    MultiFab* gv_pointer;
    MultiFab* gv_pointer_old;
    if (OnSameGrids(lev,gv))
    {
        gv_pointer     = 0;
        gv_pointer_old = 0;
    }
    else 
    {
        gv_pointer     = new MultiFab(m_gdb->ParticleBoxArray(lev),1,gv.nGrow(),
				      m_gdb->ParticleDistributionMap(lev), Fab_allocate);
        gv_pointer_old = new MultiFab(m_gdb->ParticleBoxArray(lev),1,gv_old.nGrow(),
				      m_gdb->ParticleDistributionMap(lev), Fab_allocate);
        gv_pointer->copy(gv,0,0,BL_SPACEDIM);
        gv_pointer_old->copy(gv_old,0,0,BL_SPACEDIM);
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
        gv_pointer_old->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }


    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid     = pmap_it->first;
        PBox&            pbox     = pmap_it->second;
        const int        n        = pbox.size();
        const FArrayBox& gfab     = (gv_pointer) ? (*gv_pointer)[grid] : gv[grid];
        const FArrayBox& gfab_old = (gv_pointer_old) ? (*gv_pointer_old)[grid] : gv_old[grid];
        const Real       half_dt  = Real(0.5) * dt;

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

                ParticleBase::Reset(p,m_gdb,true);
            }
        }
    }

    if (gv_pointer) delete gv_pointer;
    if (gv_pointer_old) delete gv_pointer_old;

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
ParticleContainer<N>::moveKickDrift (MultiFab& grav_vector,
                                     int             lev,
                                     Real            dt,
                                     Real            a_old,
                                     Real            a_half) 
{
    BL_PROFILE("ParticleContainer::moveKickDrift()");
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);
    BL_ASSERT(grav_vector.nGrow() >= 2);

    //If there are no particles at this level
    if (lev >= m_particles.size())
        return;

    const Real strttime      = ParallelDescriptor::second();
    const Real half_dt       = Real(0.5) * dt;
    const Real a_half_inv    = 1 / a_half;
    const Real dt_a_half_inv = dt * a_half_inv;
    PMap&      pmap          = m_particles[lev];

    MultiFab* gv_pointer;
    if (OnSameGrids(lev, grav_vector))
    {
        gv_pointer = &grav_vector;
    }
    else
    {
        gv_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),grav_vector.nComp(),grav_vector.nGrow(),
				  m_gdb->ParticleDistributionMap(lev),Fab_allocate);
        for (MFIter mfi(*gv_pointer); mfi.isValid(); ++mfi)
            gv_pointer->setVal(0.);
        gv_pointer->copy(grav_vector,0,0,grav_vector.nComp());
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = (*gv_pointer)[grid];

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

            ParticleBase::GetGravity(gfab, m_gdb->Geom(p.m_lev), p, grav);
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
        }
    }

    if (gv_pointer != &grav_vector) delete gv_pointer;

    if (lev > 0 && m_gdb->subCycle())
    {
        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&            pbox = pmap_it->second;
            const int        n    = pbox.size();
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.m_id <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!ParticleBase::RestrictedWhere(p,m_gdb, grav_vector.nGrow()-2))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.m_id == GhostParticleID)
                    {
                        p.m_id = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.m_id << std::endl;
                        BoxLib::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
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
}

//
// This version takes as input the acceleration vector at cell centers, and has the option of
//      returning the acceleration at the particle location in the data array, starting at
//      component start_comp_for_accel
//
template <int N>
void
ParticleContainer<N>::moveKick (MultiFab& grav_vector,
                                int             lev,
                                Real            dt,
                                Real            a_new,
                                Real            a_half, 
                                int             start_comp_for_accel)
{
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real strttime  = ParallelDescriptor::second();
    const Real half_dt   = Real(0.5) * dt;
    const Real a_new_inv = 1 / a_new;
    PMap&      pmap      = m_particles[lev];

    MultiFab* gv_pointer;
    if (OnSameGrids(lev,grav_vector))
    {
        gv_pointer = &grav_vector;
    }
    else 
    {
        gv_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),grav_vector.nComp(),grav_vector.nGrow(),
				  m_gdb->ParticleDistributionMap(lev),Fab_allocate);
        for (MFIter mfi(*gv_pointer); mfi.isValid(); ++mfi)
            gv_pointer->setVal(0.);
        gv_pointer->copy(grav_vector,0,0,grav_vector.nComp());
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = (*gv_pointer)[grid];

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

                ParticleBase::GetGravity(gfab, m_gdb->Geom(p.m_lev), p, grav);
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

                if (start_comp_for_accel > BL_SPACEDIM)
                {
                   D_TERM(p.m_data[start_comp_for_accel  ] = grav[0];,
                          p.m_data[start_comp_for_accel+1] = grav[1];,
                          p.m_data[start_comp_for_accel+2] = grav[2];);
                }
            }
        }
    }

    if (gv_pointer != &grav_vector) delete gv_pointer;

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
                                     Real            dt,
                                     Real            a_old,
                                     Real            a_half) 
{
    BL_PROFILE("ParticleContainer::moveKickDrift()");
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real      strttime      = ParallelDescriptor::second();
    const Geometry& geom          = m_gdb->Geom(lev);
    const Real*     dx            = geom.CellSize();
    const Real*     plo           = geom.ProbLo();
    const Real      half_dt       = Real(0.5) * dt;
    const Real      a_half_inv    = 1 / a_half;
    const Real      dt_a_half_inv = dt * a_half_inv;
    PMap&           pmap          = m_particles[lev];

    MultiFab* gv_pointer[BL_SPACEDIM];
    if (OnSameGrids(lev,grav_vector[0]))
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
           gv_pointer[i] = &grav_vector[i];
    }
    else 
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
           gv_pointer[i] = new MultiFab(m_gdb->ParticleBoxArray(lev),grav_vector[i].nComp(),
                                        grav_vector[i].nGrow(),
					m_gdb->ParticleDistributionMap(lev),Fab_allocate);
           gv_pointer[i]->copy(grav_vector[i],0,0,grav_vector[i].nComp());
           gv_pointer[i]->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
        }
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int grid = pmap_it->first;
        PBox&     pbox = pmap_it->second;
        const int n    = pbox.size();

        const FArrayBox* gfab[BL_SPACEDIM] = 
          { D_DECL(&(*gv_pointer)[0][grid],&(*gv_pointer)[1][grid],&(*gv_pointer)[2][grid]) };

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
        }
    }

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	if (gv_pointer[i] != &grav_vector[i]) delete gv_pointer[i];
    }

    if (lev > 0 && m_gdb->subCycle())
    {
        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&            pbox = pmap_it->second;
            const int        n    = pbox.size();
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.m_id <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!ParticleBase::RestrictedWhere(p,m_gdb, grav_vector[0].nGrow()-2))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.m_id == GhostParticleID)
                    {
                        p.m_id = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.m_id << std::endl;
                        BoxLib::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
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
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    const Real      strttime  = ParallelDescriptor::second();
    const Geometry& geom      = m_gdb->Geom(lev);
    const Real*     dx        = geom.CellSize();
    const Real*     plo       = geom.ProbLo();
    const Real      half_dt   = Real(0.5) * dt;
    const Real      a_new_inv = 1 / a_new;
    PMap&           pmap      = m_particles[lev];

    MultiFab* gv_pointer[BL_SPACEDIM];
    if (OnSameGrids(lev,grav_vector[0]))
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
           gv_pointer[i] = &grav_vector[i];
    }
    else 
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
           gv_pointer[i] = new MultiFab(m_gdb->ParticleBoxArray(lev),grav_vector[i].nComp(),
                                        grav_vector[i].nGrow(),
					m_gdb->ParticleDistributionMap(lev),Fab_allocate);
           gv_pointer[i]->copy(grav_vector[i],0,0,grav_vector[i].nComp());
           gv_pointer[i]->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
        }
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int grid = pmap_it->first;
        PBox&     pbox = pmap_it->second;
        const int n    = pbox.size();

        const FArrayBox* gfab[BL_SPACEDIM] = 
          { D_DECL(&(*gv_pointer)[0][grid],&(*gv_pointer)[1][grid],&(*gv_pointer)[2][grid]) };

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

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	if (gv_pointer[i] != &grav_vector[i]) delete gv_pointer[i];
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
    BL_PROFILE("ParticleContainer<N>::RemoveParticlesNotAtFinestLevel()");
    BL_ASSERT(m_gdb->finestLevel()+1 == m_particles.size());

    int cnt = 0;

    for (int lev = 0; lev < m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = m_particles[lev];

        if (!pmap.empty())
        {
            for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
                 pmap_it != pmapEnd;
                 ++pmap_it)
            {
                cnt += pmap_it->second.size();
            }

            PMap().swap(pmap);

            BL_ASSERT(pmap.empty());
        }
    }
    //
    // Print how many particles removed on each processor if any were removed.
    //
    if (m_verbose > 1)
    {
        int maxcnt = cnt;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceIntMax(maxcnt);

        if (maxcnt > 0)
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
#ifdef BL_LAZY
	});
#endif
    }
}

template <int N>
void
ParticleContainer<N>::RemoveParticlesAtLevel (int level)
{
    BL_PROFILE("ParticleContainer<N>::RemoveParticlesAtLevel()");
    if (level >= m_particles.size())
        return;

    if (!m_particles[level].empty())
    {
        PMap().swap(m_particles[level]);
    }

    BL_ASSERT(m_particles[level].empty());
}

template <int N>
void
ParticleContainer<N>::AddParticlesAtLevel (int   level,
                                           PBox& virts,
                                           bool  where_already_called)
{
    BL_PROFILE("ParticleContainer<N>::AddParticlesAtLevel()");
    if (m_particles.size() < level+1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AddParticlesAtLevel resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << level+1 << '\n';
        }
        m_particles.resize(level + 1);
    }

    const int MyProc = ParallelDescriptor::MyProc();
    //
    // The valid particles that we don't own.
    //
    PMap not_ours;

    while (!virts.empty())
    {
        ParticleType& p = virts.back();

        if (p.m_id > 0)
        {
            if (!where_already_called)
            {
                //
                // Put the particle in this level.
                //
                p.m_lev = level;

                if (!ParticleBase::SingleLevelWhere(p, m_gdb, level))
                    //
                    // Virtuals shouldn't be in Ghost cells.
                    //
                    BoxLib::Abort("ParticleContainer<N>::AddParticlesAtLevel(): Can't add outside of domain\n");
            }
            else
            {
                BL_ASSERT(p.m_lev == level);
            }

            const int who = m_gdb->ParticleDistributionMap(p.m_lev)[p.m_grid];

            if (who == MyProc)
            {
                m_particles[p.m_lev][p.m_grid].push_back(p);
            }
            else
            {
                not_ours[who].push_back(p);
            }
        }

        virts.pop_back();
    }

    if (ParallelDescriptor::NProcs() == 1)
    {
        BL_ASSERT(not_ours.empty());
    }
    else
    {
        RedistributeMPI(not_ours);
    }
}

template <int N>
void
ParticleContainer<N>::CreateVirtualParticles (int   level,
                                              PBox& virts) const
{
    BL_PROFILE("ParticleContainer<N>::CreateVirtualParticles()");
    BL_ASSERT(level > 0);
    BL_ASSERT(virts.empty());

    if (level >= m_particles.size())
        //
        // This level could exist and simply have no particles.
        //
        return;
    //
    // Read these from the parm file if we haven't done so yet.
    //
    if (aggregation_type == "")
    {
        ParmParse pp("particles");
        aggregation_type = "None";
        pp.query("aggregation_type",aggregation_type);
        aggregation_buffer = 2;
        pp.query("aggregation_buffer",aggregation_buffer);
    }
    //
    // Create a buffer so that particles near the cf border are not aggregated.
    //
    BoxArray buffer = BoxLib::complementIn(m_gdb->Geom(level).Domain(), m_gdb->boxArray(level));

    buffer.grow(aggregation_buffer);

    const PMap& pmap = m_particles[level];

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it)
    {
        const PBox& pbox = pmap_it->second;
        //
        // Map for use in Cell aggregation.
        //
        std::map<IntVect,ParticleType,IntVect::Compare> agg_map;

        typename std::map<IntVect,ParticleType,IntVect::Compare>::iterator agg_map_it;

        for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
             it != pboxEnd;
             ++it)
        {
            if (buffer.contains(it->m_cell))
            {
                //
                // It's in the no-aggregation buffer.
                //
                virts.push_back(*it);
                //
                // Set its id to indicate that it's a virt.
                //
                virts.back().m_id = VirtualParticleID;
            }
            else
            {
                if (aggregation_type == "None")
                {
                    //
                    // No aggregation.  Simply clone the particle.
                    //
                    virts.push_back(*it);
                    //
                    // Set its id to indicate that it's a virt.
                    //
                    virts.back().m_id = VirtualParticleID;
                }
                else if (aggregation_type == "Cell")
                {
                    //
                    // Note that Cell aggregation assumes that p.m_data[0] is mass and
                    // that all other components should be combined in a mass-weighted
                    // average.
                    //
                    agg_map_it = agg_map.find(it->m_cell);

                    if (agg_map_it == agg_map.end())
                    {
                        //
                        // Add the particle.
                        //
                        ParticleType p = *it;
                        //
                        // Set its id to indicate that it's a virt.
                        //
                        p.m_id = VirtualParticleID;
                        agg_map[p.m_cell] = p;
                    }
                    else
                    {
                        BL_ASSERT(agg_map_it != agg_map.end());
                        const ParticleType&  pnew       = *it;
                        ParticleType&        pold       = agg_map_it->second;
                        const Real           old_mass   = pold.m_data[0];
                        const Real           new_mass   = pnew.m_data[0];
                        const Real           total_mass = old_mass + new_mass;
                        //
                        // Set the position to the center of mass.
                        //
                        for (int i = 0; i < BL_SPACEDIM; i++)
                        {
                            pold.m_pos[i] = (old_mass*pold.m_pos[i] + new_mass*pnew.m_pos[i])/total_mass;
                        }
                        BL_ASSERT(ParticleBase::Index(pold,m_gdb->Geom(level)) == it->m_cell);
                        //
                        // Set the metadata (presumably velocity) to the mass-weighted average.
                        //
                        for (int i = 1; i < N; i++)
                        {
                            pold.m_data[i] = (old_mass*pold.m_data[i] + new_mass*pnew.m_data[i])/total_mass;
                        }
                        pold.m_data[0] = total_mass;
                    }
                }
                else if (aggregation_type == "Flow")
                {
                    BoxLib::Abort("Flow aggregation not implemented");
                }
                else 
                {
                    BoxLib::Abort("Unknown Particle Aggregation mode");
                }
            }
        }
        if (aggregation_type == "Cell")
        {
            //
            // Add the aggregated particles to the virtuals.
            //
            for (typename std::map<IntVect,ParticleType>::iterator agg_it = agg_map.begin(), aggEnd = agg_map.end(); agg_it != aggEnd; ++agg_it)
            {
                virts.push_back((*agg_it).second);
            }
        }
    }
}

template <int N>
void
ParticleContainer<N>::CreateGhostParticles (int   level,
                                            int   ngrow,
                                            PBox& ghosts) const
{
    BL_PROFILE("ParticleContainer<N>::CreateGhostParticles()");
    BL_ASSERT(ghosts.empty());
    BL_ASSERT(level < m_gdb->finestLevel());

    if (level >= m_particles.size())
        //
        // This level could exist and simply have no particles.
        //
        return;

    const BoxArray& fine = m_gdb->ParticleBoxArray(level + 1);
    
    std::vector< std::pair<int,Box> > isects;

    const PMap& pmap = m_particles[level];

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it)
    {
        const PBox& pbox = pmap_it->second;

        for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
             it != pboxEnd;
             ++it)
        {
            //
            // Find particle location on the finer level.
            //
            const IntVect& iv = ParticleBase::Index(*it,m_gdb->Geom(level+1));
            //
            // Is it in the grown finer level?
            //
            fine.intersections(Box(iv,iv),isects,false,ngrow);
            //
            // Here we add the particle to each potential grid.
            //
            for (int i = 0; i < isects.size(); i++)
            {
                //
                // Create a copy.
                //
                ParticleType p = *it;
                //
                // Set its id to indicate that it's a ghost.
                //
                p.m_id = GhostParticleID;
                //
                // Set its position.
                //
                p.m_lev  = level + 1;
                p.m_grid = isects[i].first;
                p.m_cell = iv;
                //
                // Store it in the PBox.
                //
                ghosts.push_back(p);
            }
        }
    }
}
    
//
// Uses midpoint method to advance particles using umac.
//

template <int N>
void
ParticleContainer<N>::AdvectWithUmac (MultiFab* umac,
                                      int       lev,
                                      Real      dt,
                                      int       vcomp)
{
    BL_PROFILE("ParticleContainer<N>::AdvectWithUmac()");
    BL_ASSERT(OK(true, lev, umac[0].nGrow()-1));
    BL_ASSERT(vcomp >= 0);
    BL_ASSERT(N >= vcomp + BL_SPACEDIM);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    D_TERM(BL_ASSERT(umac[0].nGrow() >= 1);,
           BL_ASSERT(umac[1].nGrow() >= 1);,
           BL_ASSERT(umac[2].nGrow() >= 1););

    D_TERM(BL_ASSERT(!umac[0].contains_nan());,
           BL_ASSERT(!umac[1].contains_nan());,
           BL_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();
    const Real*     plo      = geom.ProbLo();

    PArray<MultiFab> umac_pointer;
    // We assume that if umac[0]'s boxArray matches then the others will too...
    if (OnSameGrids(lev, umac[0]))
    {
	umac_pointer.resize(BL_SPACEDIM, PArrayNoManage);
        for (int i = 0; i < BL_SPACEDIM; i++)
	    umac_pointer.set(i, &umac[i]);
    }
    else
    {
	umac_pointer.resize(BL_SPACEDIM, PArrayManage);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
	    int ng = umac[i].nGrow();
	    
	    umac_pointer.set(i, new MultiFab(m_gdb->ParticleBoxArray(lev),
					     umac[i].nComp(),
					     ng,
					     m_gdb->ParticleDistributionMap(lev),
					     Fab_allocate,
					     IntVect::TheDimensionVector(i)));
	    umac_pointer[i].copy(umac[i],0,0,umac[i].nComp(),ng,ng);
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            const int grid = pmap_it->first;
            PBox&     pbox = pmap_it->second;
            const int n    = pbox.size();

            FArrayBox* fab[BL_SPACEDIM] = { D_DECL(&umac_pointer[0][grid],
                                                   &umac_pointer[1][grid],
                                                   &umac_pointer[2][grid]) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                
                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

                const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + Real(0.5),
                                                       (p.m_pos[1]-plo[1])/dx[1] + Real(0.5),
                                                       (p.m_pos[2]-plo[2])/dx[2] + Real(0.5)) };

                const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

                const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    IntVect ecell = cell;

                    ecell[d] = p.m_cell[d] + 1;

                    Real efrac[BL_SPACEDIM] = { D_DECL(frac[0], frac[1], frac[2]) };

                    efrac[d] = (p.m_pos[d]-plo[d])/dx[d] - p.m_cell[d];

                    for (int j = 0; j < BL_SPACEDIM; j++)
                    {
                        if (efrac[j] > 1) efrac[j] = 1;
                        if (efrac[j] < 0) efrac[j] = 0;
                    }

                    const Real vel = ParticleBase::InterpDoit(*fab[d], ecell, efrac, 0);

                    if (ipass == 0)
                    {
                        //
                        // Save old position and the vel & predict location at dt/2.
                        //
                        p.m_data[vcomp+d] = p.m_pos[d];
                        p.m_pos[d] += 0.5*dt*vel;
                    }
                    else
                    {
                        //
                        // Update to final time using the orig position and the vel at dt/2.
                        //
                        p.m_pos[d]  = p.m_data[vcomp+d] + dt*vel;
                        // Save the velocity for use in Timestamp().
			p.m_data[vcomp+d] = vel;
                    }
                }
                
                ParticleBase::RestrictedWhere(p,m_gdb, umac[0].nGrow()); 
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AdvectWithUmac() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}

//
// Uses midpoint method to advance particles using cell-centered velocity
//

template <int N>
void
ParticleContainer<N>::AdvectWithUcc (const MultiFab& Ucc,
				     int             lev,
				     Real            dt,
				     int             vcomp)
{
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(true, lev, Ucc.nGrow()-1));
    BL_ASSERT(vcomp >= 0);
    BL_ASSERT(N >= vcomp + BL_SPACEDIM);
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    BL_ASSERT(!Ucc.contains_nan());

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    BL_ASSERT(OnSameGrids(lev,Ucc));

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            const int grid = pmap_it->first;
            PBox&     pbox = pmap_it->second;
            const int n    = pbox.size();

	    const FArrayBox& fab = Ucc[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                
                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

		Real v[BL_SPACEDIM];

		ParticleBase::Interp(p, geom, fab, idx, v, BL_SPACEDIM);

		if (ipass == 0) {
		    //
		    // Save old position and the vel & predict location at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
			p.m_data[vcomp+d] = p.m_pos[d];
                        p.m_pos[d] += 0.5*dt*v[d];
                    }
		} else {
		    //
		    // Update to final time using the orig position and the vel at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
                        p.m_pos[d]  = p.m_data[vcomp+d] + dt*v[d];
                        // Save the velocity for use in Timestamp().
			p.m_data[vcomp+d] = v[d];
                    }
                }
                
                ParticleBase::RestrictedWhere(p,m_gdb, Ucc.nGrow()); 
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AdvectWithUcc() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}

//
// This redistributes valid particles and discards invalid ones.
//

template <int N>
void
ParticleContainer<N>::Redistribute (bool where_already_called,
                                    bool full_where,
                                    int  lev_min,
                                    int  nGrow)
{
    BL_PROFILE("ParticleContainer::Redistribute()");
    const int MyProc   = ParallelDescriptor::MyProc();
    Real      strttime = ParallelDescriptor::second();
    //
    // On startup there are cases where Redistribute() could be called
    // with a given finestLevel() where that AmrLevel has yet to be defined.
    //
    int theEffectiveFinestLevel = m_gdb->finestLevel();

    while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
        theEffectiveFinestLevel--;

    if (m_particles.size() < theEffectiveFinestLevel+1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Redistribute() resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << theEffectiveFinestLevel+1 << '\n';
        }
        m_particles.resize(theEffectiveFinestLevel+1);
    }
    //
    // The valid particles that we don't own.
    //
    PMap not_ours;

    for (int lev = lev_min; lev < m_particles.size(); lev++)
    {
        PMap& pmap = m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            const int grid = pmap_it->first;
            PBox&     pbox = pmap_it->second;

            for (typename PBox::iterator it = pbox.begin(), End = pbox.end(); it != End; )
            {
                ParticleType& p = *it;

                if (p.m_id > 0)
                {
                    if (!where_already_called)
                    {
                        if (!ParticleBase::Where(p,m_gdb, lev_min, theEffectiveFinestLevel))
                        {                                
                            if (full_where) // Lengthier checks for subcycling.
                            {
                                if (!ParticleBase::PeriodicWhere(p, m_gdb, lev_min, theEffectiveFinestLevel))
                                {
                                    if (lev_min != 0) // RestrictedWhere should be unnecessary at top level.
                                    {
                                        if (!ParticleBase::RestrictedWhere(p, m_gdb, nGrow))
                                            BoxLib::Abort("ParticleContainer<N>::Redistribute(): invalid particle at non-coarse step");
                                    }
                                    else
                                    {
                                        //
                                        // The particle has left the domain; invalidate it.
                                        // This typically only happens on a coarse timestep.
                                        //
                                        p.m_id = -p.m_id;
                                    }
                                }
                            }
                            else
                            {
                                std::cout << "Bad Particle: " << p << '\n';
                                BoxLib::Abort("ParticleContainer<N>::Redistribute(): invalid particle in basic check");
                            }
                        }
                    }

                    if (p.m_id > 0)
                    {
                        //
                        // The owner of the particle is the CPU owning the finest grid
                        // in state data that contains the particle.
                        //
                        const int who = m_gdb->ParticleDistributionMap(p.m_lev)[p.m_grid];

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
                            not_ours[who].push_back(p);
                            //
                            // Invalidate the particle so we can reclaim its space.
                            //
                            p.m_id = -p.m_id;
                        }
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
            {
                pmap.erase(pmap_it++);
            }
            else
            {
                ++pmap_it;
            }
        }
    }

    if (m_particles.size() > theEffectiveFinestLevel+1)
    {
        //
        // Looks like we lost an AmrLevel on a regrid.
        //
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Redistribute() resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << theEffectiveFinestLevel+1 << '\n';
        }
        BL_ASSERT(m_particles.size() >= 2);
        BL_ASSERT(m_particles[m_particles.size()-1].empty());

        m_particles.resize(theEffectiveFinestLevel+1);
    }

    if (ParallelDescriptor::NProcs() == 1)
    {
        BL_ASSERT(not_ours.empty());
    }
    else
    {
        RedistributeMPI(not_ours);
    }

    BL_ASSERT(OK(full_where, lev_min, nGrow, theEffectiveFinestLevel));

    if (m_verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ByteSpread();

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "ParticleContainer<N>::Redistribute() time: " << stoptime << "\n\n";
#ifdef BL_LAZY
	});
#endif
    }
}

template <int N>
void
ParticleContainer<N>::RedistributeMPI (PMap& not_ours)
{
    BL_PROFILE("ParticleContainer<N>::RedistributeMPI()");
#if BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    //
    // We may now have particles that are rightfully owned by another CPU.
    //
    Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

    int NumSnds = 0, NumRcvs = 0;

    for (typename PMap::const_iterator it = not_ours.begin(), End = not_ours.end(); it != End; ++it)
    {
        NumSnds        += it->second.size();
        Snds[it->first] = it->second.size();
    }

    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0)
        //
        // There's no parallel work to do.
        //
        return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++)
        if (Snds[i] > 0)
            SndCnts[i] = Snds[i];

    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
        }
    }
    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);
    //
    // We'll store the particles we're to receive in a PMap indexed by proc # of receiver.
    //
    PMap nparticles;

    for (IntIntMap::const_iterator it = RcvCnts.begin(), End = RcvCnts.end(); it != End; ++it)
    {
        nparticles[it->first].resize(it->second);
    }
    Array<int>         owner(RcvCnts.size());
    Array<int>         index(RcvCnts.size());
    Array<MPI_Status>  stats(RcvCnts.size());
    Array<MPI_Request> rreqs(RcvCnts.size());
    //
    // First send/recv the integer parts of the particles.
    //
    {
        const int SeqNum     = ParallelDescriptor::SeqNum();
        const int iChunkSize = 4 + BL_SPACEDIM;
        //
        // Allocate data for rcvs as one big chunk.
        //
        Array<int> recvdata(NumRcvs * iChunkSize);
        //
        // Post receives.
        //
        int idx = 0;
        for (IntIntMap::const_iterator it = RcvCnts.begin(), End = RcvCnts.end(); it != End; ++it, ++idx)
        {
            const int Who = it->first;
            const int Cnt = it->second   * iChunkSize;
            const int Idx = rOffset[Who] * iChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            owner[idx] = Who;
            rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
        //
        // Send the integer data.
        //
        Array<int> senddata;

        for (IntIntMap::const_iterator it = SndCnts.begin(), End = SndCnts.end(); it != End; ++it)
        {
            const int Who = it->first;
            const int Cnt = it->second * iChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            senddata.resize(Cnt);

            const PBox& pbox = not_ours[Who];

            int ioff = 0;
            for (typename PBox::const_iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
            {
                const ParticleType& p = *it;

                BL_ASSERT(p.m_id > 0);

                senddata[ioff+0] = p.m_id;
                senddata[ioff+1] = p.m_cpu;

                senddata[ioff+2] = p.m_lev;
                senddata[ioff+3] = p.m_grid;

                D_TERM(senddata[ioff+4] = p.m_cell[0];,
                       senddata[ioff+5] = p.m_cell[1];,
                       senddata[ioff+6] = p.m_cell[2];);

                ioff += iChunkSize;
            }

            ParallelDescriptor::Send(senddata.dataPtr(),Cnt,Who,SeqNum);
        }
        //
        // Free up this memory ...
        //
        Array<int>().swap(senddata);
        //
        // Now receive and unpack the integer data.
        //
        for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(rreqs, completed, index, stats);

            for (int k = 0; k < completed; k++)
            {
                const int  Who  = owner[index[k]];
                const int  Idx  = rOffset[Who] * iChunkSize;
                const int* rcvp = &recvdata[Idx];
                PBox&      pbox = nparticles[Who];

                BL_ASSERT(pbox.size() == RcvCnts[Who]);

                for (typename PBox::iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
                {
                    ParticleType& p = *it;

                    BL_ASSERT(rcvp != 0);

                    p.m_id   = rcvp[0];
                    p.m_cpu  = rcvp[1];

                    p.m_lev  = rcvp[2];
                    p.m_grid = rcvp[3];

                    D_TERM(p.m_cell[0] = rcvp[4];,
                           p.m_cell[1] = rcvp[5];,
                           p.m_cell[2] = rcvp[6];);

                    rcvp += iChunkSize;
                }
            }
        }
    }
    //
    // Next send/recv the Real parts of the particles.
    //
    {
        const int SeqNum     = ParallelDescriptor::SeqNum();
        const int rChunkSize = BL_SPACEDIM+N;
        //
        // Allocate data for rcvs as one big chunk.
        //
        Array<ParticleBase::RealType> recvdata(NumRcvs * rChunkSize);
        //
        // Post receives.
        //
        int idx = 0;
        for (IntIntMap::const_iterator it = RcvCnts.begin(), End = RcvCnts.end(); it != End; ++it, ++idx)
        {
            const int Who = it->first;
            const int Cnt = it->second   * rChunkSize;
            const int Idx = rOffset[Who] * rChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
        //
        // Send the Real data.
        //
        Array<ParticleBase::RealType> senddata;

        for (IntIntMap::const_iterator it = SndCnts.begin(), End = SndCnts.end(); it != End; ++it)
        {
            const int Who = it->first;
            const int Cnt = it->second * rChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());
            
            senddata.resize(Cnt);

            PBox& pbox = not_ours[Who];

            int ioff = 0;
            for (typename PBox::const_iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
            {
                const ParticleType& p = *it;

                BL_ASSERT(p.m_id > 0);

                D_TERM(senddata[ioff+0] = p.m_pos[0];,
                       senddata[ioff+1] = p.m_pos[1];,
                       senddata[ioff+2] = p.m_pos[2];);

                ioff += BL_SPACEDIM;

                for (int j = 0; j < N; j++)
                    senddata[ioff+j] = p.m_data[j];

                ioff += N;
            }

            PBox().swap(pbox);

            ParallelDescriptor::Send(senddata.dataPtr(),Cnt,Who,SeqNum);
        }
        //
        // Free up this memory ...
        //
        Array<ParticleBase::RealType>().swap(senddata);
        //
        // Now receive and unpack the Real data.
        //
        for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(rreqs, completed, index, stats);

            for (int k = 0; k < completed; k++)
            {
                const int                     Who  = owner[index[k]];
                const int                     Idx  = rOffset[Who] * rChunkSize;
                const ParticleBase::RealType* rcvp = &recvdata[Idx];
                PBox&                         pbox = nparticles[Who];

                BL_ASSERT(pbox.size() == RcvCnts[Who]);

                for (typename PBox::iterator it = pbox.begin(), End = pbox.end(); it != End; ++it)
                {
                    ParticleType& p = *it;

                    BL_ASSERT(rcvp != 0);

                    D_TERM(p.m_pos[0] = rcvp[0];,
                           p.m_pos[1] = rcvp[1];,
                           p.m_pos[2] = rcvp[2];);

                    rcvp += BL_SPACEDIM;

                    for (int j = 0; j < N; j++)
                        p.m_data[j] = rcvp[j];

                    rcvp += N;

                    m_particles[p.m_lev][p.m_grid].push_back(p);
                }

                PBox().swap(pbox);
            }
        }
    }
#endif /*BL_USE_MPI*/
}

template <int N>
bool
ParticleContainer<N>::OK (bool full_where,
                          int  lev_min,
                          int  ngrow,
                          int  finest_level) const
{
    BL_PROFILE("ParticleContainer<N>::OK()");
    if (finest_level == -1)
        finest_level = m_gdb->finestLevel();

    BL_ASSERT(finest_level <= m_gdb->finestLevel());
    //
    // Check that the integer data in each valid particle is what it should be.
    // This includes checking that particles are in the proper place in the particle
    // container based on what Where() says they should be.
    //
    // Particles are copied to avoid accidentally moving them with where.
    //
    for (int lev = lev_min; lev < m_particles.size(); lev++)
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

                    if (!ParticleBase::Where(p,m_gdb, lev_min, finest_level))
                    {
                        if (full_where)
                        {
                            if (!ParticleBase::PeriodicWhere(p, m_gdb, lev_min, finest_level)) 
                            {
                                if (!ParticleBase::RestrictedWhere(p, m_gdb, ngrow))
                                    return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
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
                                  const std::string& name,
                                  bool               is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<N>::Checkpoint()");
    BL_ASSERT(OK());

    BL_ASSERT(sizeof(ParticleBase::RealType) == 4 || sizeof(ParticleBase::RealType) == 8);

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
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
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
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        if (sizeof(ParticleBase::RealType) == 4)
        {
            HdrFile << ParticleBase::Version() << "_single" << '\n';
        }
        else
        {
            HdrFile << ParticleBase::Version() << "_double" << '\n';
        }
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
        HdrFile << m_gdb->finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
        {
            HdrFile << m_gdb->boxArray(lev).size() << '\n';
        }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const bool gotsome = (NumberOfParticlesAtLevel(lev) > 0);
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;

        if (gotsome)
        {
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
        }

	MultiFab state(m_gdb->boxArray(lev),1,0,Fab_noallocate);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Array<int>  which(state.size(),0);
        Array<int > count(state.size(),0);
        Array<long> where(state.size(),0);

        if (gotsome)
        {
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
                    WriteParticles(lev, ParticleFile, FileNumber, which, count, where, is_checkpoint);

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
        }

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

            if (gotsome)
            {
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
ParticleContainer<N>::WritePlotFile (const std::string& dir,
                                     const std::string& name) const
{
    BL_PROFILE("ParticleContainer<N>::WritePlotFile()");
    BL_ASSERT(OK());
    bool is_checkpoint = false;

    // For yt we need exactly the chk particle format so would need to set is_checkpoint = true
    // Anyway, it's not too bad to have particle ids on disk,
    // think of merger trees or backtracing of particles for nested ics
    // is_checkpoint = true; 
    Checkpoint(dir,name,is_checkpoint);
}

template <int N>
void
ParticleContainer<N>::WriteParticles (int            lev,
                                      std::ofstream& ofs,
                                      int            fnum,
                                      Array<int>&    which,
                                      Array<int>&    count,
                                      Array<long>&   where,
                                      bool           is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<N>::WriteParticles()");
    const PMap&     pmap  = m_particles[lev];

    MultiFab state(m_gdb->boxArray(lev),1,0,Fab_noallocate);

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
            for (typename PBox::const_iterator it = pmap_it->second.begin(), pboxEnd =  pmap_it->second.end();
                 it != pboxEnd;
                 ++it)
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

        if (is_checkpoint)
        {
            //
            // First write out the integer data in binary.
            // We do not need to write out the m_lev and m_grid
            // info since it's implicit in how the particles
            // are stored.  We can easily recreate them on restart.
            //
            const int iChunkSize = 2+BL_SPACEDIM;

#ifdef BL_LOWMEMPWRITE
	    int maxItemsToWrite(8192);
	    int cntBufSize(maxItemsToWrite * iChunkSize);
	    int nItems(cnt), nItemsToWrite(0);
            Array<int> istuff(cntBufSize);
            typename PBox::const_iterator it  = pbox.begin();
            typename PBox::const_iterator End = pbox.end();

	    while(nItems > 0) {
              int *iptr = istuff.dataPtr();
	      int itemCount(0);
              for( ; it != End && itemCount < maxItemsToWrite; ++it) {
                if(it->m_id > 0) {
                    BL_ASSERT(it->m_lev == lev);
                    BL_ASSERT(it->m_grid == grid);

                    iptr[0] = it->m_id;
                    iptr[1] = it->m_cpu;

                    D_TERM(iptr[2] = it->m_cell[0];,
                           iptr[3] = it->m_cell[1];,
                           iptr[4] = it->m_cell[2];);

                    iptr += iChunkSize;
                }
		++itemCount;
              }
	      nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
              ofs.write((char *) istuff.dataPtr(), nItemsToWrite * iChunkSize * sizeof(int));
	      nItems -= nItemsToWrite;
	    }
#else
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
#endif
        }


        //
        // Write the Real data in binary.
        //
        const int rChunkSize = BL_SPACEDIM+N;

#ifdef BL_LOWMEMPWRITE
	int maxItemsToWrite(8192);
	int cntBufSize(maxItemsToWrite * rChunkSize);
	int nItems(cnt), nItemsToWrite(0);
        Array<ParticleBase::RealType> rstuff(cntBufSize);
        typename PBox::const_iterator it = pbox.begin();
        typename PBox::const_iterator End = pbox.end();

	while(nItems > 0) {
          ParticleBase::RealType *rptr = rstuff.dataPtr();
	  int itemCount(0);
          for( ; it != End && itemCount < maxItemsToWrite; ++it) {
            if(it->m_id > 0) {
                D_TERM(rptr[0] = it->m_pos[0];,
                       rptr[1] = it->m_pos[1];,
                       rptr[2] = it->m_pos[2];);

                for (int i = 0; i < N; i++) {
                  rptr[BL_SPACEDIM+i] = it->m_data[i];
		}
                rptr += rChunkSize;
            }
	    ++itemCount;
          }

	  nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
          ofs.write((char *) rstuff.dataPtr(), nItemsToWrite * rChunkSize * sizeof(ParticleBase::RealType));
	  nItems -= nItemsToWrite;
	}
#else
        Array<ParticleBase::RealType> rstuff(cnt*rChunkSize);

        ParticleBase::RealType* rptr = rstuff.dataPtr();

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

        ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(ParticleBase::RealType));
#endif
    }
}

template <int N>
void
ParticleContainer<N>::Restart (const std::string& dir,
                               const std::string& file,
                               bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<N>::Restart()");
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
    //
    // What do our version strings mean?
    //
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // 
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    //
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    //
    version = vbuf.dataPtr();

    if (version.find("Version_One_Dot_Zero") != std::string::npos)
    {
        Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
    }
    else if (version.find("Version_One_Dot_One") != std::string::npos)
    {
        if (version.find("_single") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"single",is_checkpoint);
        }
        else if (version.find("_double") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
        }
        else
        {
            std::string msg("ParticleContainer<N>::Restart(): bad version string: ");
            msg += version;
            BoxLib::Error(version.c_str());
        }
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
ParticleContainer<N>::Restart_Doit (const std::string& fullname,
                                    std::ifstream&     HdrFile,
                                    const std::string& how,
                                    bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<N>::RestartDoit()");
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

    BL_ASSERT(finest_level == m_gdb->finestLevel());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            HdrFile >> ngrids[lev];

            BL_ASSERT(ngrids[lev] > 0);
            BL_ASSERT(ngrids[lev] == m_gdb->boxArray(lev).size());
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

        m_particles.resize(m_gdb->finestLevel()+1);

	MultiFab state(m_gdb->boxArray(lev),1,0,Fab_noallocate);

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

            if (how == "single")
            {
                ReadParticles_SinglePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else if (how == "double")
            {
                ReadParticles_DoublePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else
            {
                std::string msg("ParticleContainer<N>::Restart_Doit(): bad parameter: ");
                msg += how;
                BoxLib::Error(msg.c_str());
            }
                
            ParticleFile.close();

            if (!ParticleFile.good())
                BoxLib::Abort("ParticleContainer<N>::Restart_Doit(): problem reading particles");
        }
    }

    BL_ASSERT(OK());        
}

//
// This one stores real data as doubles.
//

template <int N>
void
ParticleContainer<N>::ReadParticles_DoublePrecision (int            cnt,
                                                     int            grd,
                                                     int            lev,
                                                     bool           is_checkpoint,
                                                     std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<N>::ReadParticles_DoublePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < m_particles.size());
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->boxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2+BL_SPACEDIM;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the double data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+N;

    Array<double> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(double));
    //
    // Now reassemble the particles.
    //
    int*            iptr = istuff.dataPtr();
    double*         rptr = rstuff.dataPtr();
    PBox&           pbox = m_particles[lev][grd];
    const Geometry& geom = m_gdb->Geom(0);

    const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
    const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };
    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleBase::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        if (is_checkpoint)
        {
            p.m_id   = iptr[0];
            p.m_cpu  = iptr[1];
        }
        else
        {
           if (!ParticleBase::Where(p,m_gdb))
            {
                ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                {
                    std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                    }

                    BoxLib::Abort("ParticleContainer<N>::ReadParticles_DoublePrecision(): invalid particle");
                }
            }

            p.m_id   = ParticleBase::NextID();
            p.m_cpu  = ParallelDescriptor::MyProc();
        }
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
        //
        // If we're reading in doubles and storing'm in floats we have
        // to make sure the particles stay in the domain.
        //
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            if (p.m_pos[d] <= ProbLo[d]) p.m_pos[d] += Delta[d];
            if (p.m_pos[d] >= ProbHi[d]) p.m_pos[d] -= Delta[d];
        }

        for (int i = 0; i < N; i++)
            p.m_data[i] = rptr[BL_SPACEDIM+i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

//
// This one stores real data as floats.
//

template <int N>
void
ParticleContainer<N>::ReadParticles_SinglePrecision (int            cnt,
                                                     int            grd,
                                                     int            lev,
                                                     bool           is_checkpoint,
                                                     std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<N>::ReadParticles_SinglePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < m_particles.size());
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->boxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2+BL_SPACEDIM;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the float data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+N;

    Array<float> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(float));
    //
    // Now reassemble the particles.
    //
    int*   iptr = istuff.dataPtr();
    float* rptr = rstuff.dataPtr();
    PBox&  pbox = m_particles[lev][grd];

    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleBase::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        p.m_id   = iptr[0];
        p.m_cpu  = iptr[1];

        if (is_checkpoint)
        {
            p.m_id   = iptr[0];
            p.m_cpu  = iptr[1];
        }
        else
        {
           if (!ParticleBase::Where(p,m_gdb))
            {
                ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                {
                    std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                    }

                    BoxLib::Abort("ParticleContainer<N>::ReadParticles_SinglePrecision(): invalid particle");
                }
            }

            p.m_id   = ParticleBase::NextID();
            p.m_cpu  = ParallelDescriptor::MyProc();
        }

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
                                 const std::vector<int>& indices,
                                 const int               vcomp)
{
    BL_PROFILE("ParticleContainer<N>::Timestamp()");
    //
    // basename -> base filename for the output file
    // mf       -> the multifab
    // lev      -> level to check for particles
    // time     -> simulation time (will be recorded in Timestamp file)
    // indices  -> indices into mf that we output
    // vcomp    -> the offset of velocity terms in the particle metadata
    //
    BL_ASSERT(lev >= 0);
    BL_ASSERT(time >= 0);
    BL_ASSERT(vcomp >= 0);
    BL_ASSERT(!basename.empty());
    BL_ASSERT(N >= vcomp + BL_SPACEDIM);
    BL_ASSERT(lev <= m_gdb->finestLevel());

    const Real strttime = ParallelDescriptor::second();

    const int   MyProc    = ParallelDescriptor::MyProc();
    const int   NProcs    = ParallelDescriptor::NProcs();
    // We'll spread the output over this many files.
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
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

                const int       M  = indices.size();
                const BoxArray& ba = mf.boxArray();

                std::vector<Real> vals(M);

                for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
                     pmap_it != pmapEnd;
                     ++pmap_it)
                {
                    const int        grid = pmap_it->first;
                    const PBox&      pbox = pmap_it->second;
                    const Box&       bx   = ba[grid];
                    const FArrayBox& fab  = mf[grid];

                    for (typename PBox::const_iterator it = pbox.begin(), pboxEnd = pbox.end();
                         it != pboxEnd;
                         ++it)
                    {
                        const ParticleType& p = *it;

                        if (p.m_id <= 0) continue;

                        const IntVect& iv = ParticleBase::Index(p, m_gdb->Geom(lev));

                        if (!bx.contains(iv) && !ba.contains(iv)) continue;

                        BL_ASSERT(p.m_lev == lev);
                        BL_ASSERT(p.m_grid == grid);

                        TimeStampFile << p.m_id  << ' ' << p.m_cpu << ' ';

                        D_TERM(TimeStampFile << p.m_pos[0] << ' ';,
                               TimeStampFile << p.m_pos[1] << ' ';,
                               TimeStampFile << p.m_pos[2] << ' ';);

                        TimeStampFile << time;
                        //
                        // AdvectWithUmac stores the velocity in m_data ...
                        //
                        D_TERM(TimeStampFile << ' ' << p.m_data[vcomp+0];,
                               TimeStampFile << ' ' << p.m_data[vcomp+1];,
                               TimeStampFile << ' ' << p.m_data[vcomp+2];);

                        if (M > 0)
                        {
                            ParticleBase::Interp(p,m_gdb->Geom(p.m_lev),fab,&indices[0],&vals[0],M);

                            for (int i = 0; i < M; i++)
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

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Timestamp: lev: " << lev << " time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
        });
#endif
    }
}

template <int N>
void
ParticleContainer<N>::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer<N>::WriteAsciiFile()");
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

template <int N>
void
ParticleContainer<N>::WriteCoarsenedAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer<N>::WriteCoarsenedAsciiFile()");
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
                // Only keep particles in even cells
                if (it->m_id > 0 &&
                   (it->m_cell[0])%2 == 0 && (it->m_cell[1])%2 == 0 && (it->m_cell[2])%2 == 0)
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
            BoxLib::Abort("ParticleContainer<N>::WriteCoarsenedAsciiFile(): problem writing file");
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
                        // Only keep particles in even cells
                        if (it->m_id > 0 &&
                           (it->m_cell[0])%2 == 0 && (it->m_cell[1])%2 == 0 && (it->m_cell[2])%2 == 0)
                        {
                            D_TERM(File << it->m_pos[0] << ' ',
                                        << it->m_pos[1] << ' ',
                                        << it->m_pos[2] << ' ');

                            // Multiply mass by 8 since we are only taking 1/8 of the total particles
                            //   and want to keep the mass in the domain the same.
                            File << (8.0 * it->m_data[0]) << ' ';

                            // Now write out the velocity of the particle chosen; no weighting here
                            for (int i = 1; i < N; i++)
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
                BoxLib::Abort("ParticleContainer<N>::WriteCoarsenedAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::WriteCoarsenedAsciiFile() time: " << stoptime << '\n';
        }
    }
}

inline
void
ParticleBase::CIC_Fracs (const Real* frac, Real* fracs)
{
    //
    // "frac"  should be dimensioned: Real frac[BL_SPACEDIM]
    //
    // "fracs" should be dimensioned: Real fracs[D_TERM(2,+2,+4)]
    //
#if (BL_SPACEDIM == 1)
    // High
    fracs[0] = frac[0];

    // Low
    fracs[1] = (1-frac[0]);

#elif (BL_SPACEDIM == 2)
    // HH
    fracs[0] = frac[0] * frac[1] ;
    
    // LH
    fracs[1] = (1-frac[0]) * frac[1];
    
    // LL
    fracs[2] = (1-frac[0]) * (1-frac[1]);
    
    // HL
    fracs[3] = frac[0] * (1-frac[1]);

#elif (BL_SPACEDIM == 3)
    // HHH
    fracs[0] = frac[0] * frac[1] * frac[2];

    // LHH
    fracs[1] = (1-frac[0]) * frac[1] * frac[2];

    // LLH
    fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
    
    // HLH
    fracs[3] = frac[0] * (1-frac[1]) * frac[2];

    // HHL
    fracs[4] = frac[0] * frac[1] * (1-frac[2]);
    
    // LHL
    fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);

    // LLL
    fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
    
    // HLL
    fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
#endif
}

inline
void
ParticleBase::CIC_Cells (const IntVect& hicell, IntVect* cells)
{
    //
    // "cells" should be dimensioned: IntVect cells[D_TERM(2,+2,+4)]
    //
    IntVect cell = hicell;

#if (BL_SPACEDIM == 1)
    // High
    cells[0] = cell;

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

#elif (BL_SPACEDIM == 2)
    // HH
    cells[0] = cell;
    
    // LH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    
    // LL
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HL
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

#elif (BL_SPACEDIM == 3)
    // HHH
    cells[0] = cell;

    // LHH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

    // LLH
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HLH
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

    cell = hicell;

    // HHL
    cell[2]  = cell[2] - 1;
    cells[4] = cell;
    
    // LHL
    cell[0]  = cell[0] - 1;
    cells[5] = cell;

    // LLL
    cell[1]  = cell[1] - 1;
    cells[6] = cell;
    
    // HLL
    cell[0]  = cell[0] + 1;
    cells[7] = cell;
#endif
}

inline
int
ParticleBase::CIC_Cells_Fracs (const ParticleBase& p,
                               const Real*         plo,
                               const Real*         dx,
                               Array<Real>&        fracs,
                               Array<IntVect>&     cells)
{
    return ParticleBase::CIC_Cells_Fracs(p,plo,dx,dx,fracs,cells);
}

#endif /*_PARTICLES_H_*/
