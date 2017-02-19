#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <cstring>
#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <memory>
#include <limits>
#include <utility>
#include <tuple>
#include <type_traits>

#include <AMReX_ParmParse.H>
#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Deposition_F.H>
#include <AMReX_RealBox.H>
#include <AMReX_Print.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {
  
namespace
{
    std::string   aggregation_type   = "";
    int           aggregation_buffer = 1;
    constexpr int GhostParticleID    = std::numeric_limits<int>::max();
    constexpr int VirtualParticleID  = std::numeric_limits<int>::max()-1;
    constexpr int LastParticleID     = std::numeric_limits<int>::max()-2;
}

//
// A struct used for communicating particle data accross processes
// during multi-level operations.
//
struct ParticleCommData
{
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif
  int     m_lev;
  int     m_grid;
  IntVect m_cell;
  RealType m_data[1 + BL_SPACEDIM];
};


//
// A struct used for storing a particle's position in the AMR hierarchy.
//
struct ParticleLocData
{
  int     m_lev  = -1;
  int     m_grid = -1;
  int     m_tile = -1;
  IntVect m_cell {D_DECL(-1,-1,-1)};
  Box     m_gridbox;
  Box     m_tilebox;
};


//
// The struct used to store particles.
//
template<int NR, int NI=0>
struct Particle
{
    //
    // The floating point type used for the particles.
    //
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif

    //
    // The real data. We always have SPACEDIM position coordinates,
    // and optionally we have NR additional real attributes.
    //
    union rm_t
    {
      RealType pos[BL_SPACEDIM];
      RealType arr[BL_SPACEDIM+NR];
    };

    //
    // The integer data. We always have id and cpu, and optionally we
    // have NI additional integer attributes.
    //
    union im_t
    {
      struct {int id; int cpu;};
      int arr[2+NI];
    };

    rm_t m_rdata;
    im_t m_idata;

    static int the_next_id;

    int&  id()       & {return m_idata.id;}
    int   id() const & {return m_idata.id;}
    int& cpu()       & {return m_idata.cpu;}
    int  cpu() const & {return m_idata.cpu;}

    RealType& pos(int index)       & {return m_rdata.pos[index];}
    RealType  pos(int index) const & {return m_rdata.pos[index];}

    RealType& rdata(int index)       & {return m_rdata.arr[BL_SPACEDIM + index];}
    RealType  rdata(int index) const & {return m_rdata.arr[BL_SPACEDIM + index];}

    int& idata(int index)       & {return m_idata.arr[2 + index];}
    int  idata(int index) const & {return m_idata.arr[2 + index];}

    static Real InterpDoit (const FArrayBox& fab, const Real* fracs, const IntVect* cells, int comp);

    static Real InterpDoit (const FArrayBox& fab, const IntVect& hi, const Real* frac, int comp);

    static void Interp (const Particle<NR, NI>& prt, const Geometry& geom, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Geometry& geom, const Particle<NR, NI>& p, Real* grav);

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

    static void CIC_Fracs (const Real* frac, Real* fracs);

    static void CIC_Cells (const IntVect& hicell, IntVect* cells);

    //
    // Old, *-based CIC for use in Interp.
    //
    static void CIC_Cells_Fracs_Basic (const Particle<NR, NI>& p, const Real* plo, const Real* dx, Real* fracs,  IntVect* cells);

    //
    // Wraps the arbitrary dx function.
    //
    static int CIC_Cells_Fracs (const Particle<NR, NI>& p, 
                                const Real*         plo, 
                                const Real*         dx, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Does CIC computations for arbitrary particle/grid dx's.
    //
    static int CIC_Cells_Fracs (const Particle<NR, NI>& p, 
                                const Real*         plo, 
                                const Real*         dx_geom, 
                                const Real*         dx_part, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
};

template <int NR, int NI>
class ArrayOfStructs {
public:
    using ParticleType = Particle<NR, NI>;
    using RealType     = typename ParticleType::RealType;

    Array<ParticleType>& operator()() { return m_data; }
    const Array<ParticleType>& operator()() const { return m_data; }

    int numParticles () const { return m_data.size(); }

    RealType*       data ()       { return &(m_data[0].m_rdata.arr[0]); }
    const RealType* data () const { return &(m_data[0].m_rdata.arr[0]); }

    std::pair<int,int> dataShape () const {
	return std::make_pair(SizeInReal, static_cast<int>(m_data.size()));
    }
    
    static constexpr int SizeInReal = sizeof(ParticleType) / sizeof(RealType);

private:
    Array<ParticleType> m_data;
};

template <int NA>
struct StructOfArrays {
  
    Array<Real>& operator[] (const int index) {
        return m_data[index];
    }
    
    const Array<Real>& operator[] (const int index) const {
        return m_data[index];
    }
    
    int numParticles () const { return NA>0 ? m_data[0].size() : 0; }

private:
    std::array<Array<Real>, NA> m_data;
};

template <int NR, int NI, int NA> class ParIter;
template <int NR, int NI, int NA> class ParConstIter;

template <int NR, int NI = 0, int NA = 0>
class ParticleContainer
{
    friend class ParIter<NR,NI,NA>;
    friend class ParConstIter<NR,NI,NA>;

public:
    //
    // The type of Particles we hold.
    //
    using ParticleType = Particle<NR, NI>;
    using RealType = typename Particle<NR, NI>::RealType;

    // This defines the array-of-structs datatype.
    using AoS = ArrayOfStructs<NR,NI>;

    // The struct-of-arrays data is stored as a double array where the first index is the
    // component number and the second is the particle index.
    using SoA = StructOfArrays<NA>;

    // A single level worth of particles is indexed first by grid id, 
    // then by tile id, for both SoA and AoS data.
    using AoSMap = std::map<std::pair<int, int>, AoS>;
    using SoAMap = std::map<std::pair<int, int>, SoA>;

    ParticleContainer ()
      : 
      communicate_comp(NA, true),
      m_verbose(0), 
      m_gdb(nullptr),
      allow_particles_near_boundary(false)
    {
        Initialize ();
    }

    ParticleContainer (ParGDBBase* gdb)
        :
        communicate_comp(NA, true),
	m_verbose(0),
	m_gdb(gdb), 
	allow_particles_near_boundary(false)
    {
        Initialize ();
        reserveData();
        resizeData();
    }
  
    ParticleContainer (const Geometry            & geom, 
		       const DistributionMapping & dmap,
		       const BoxArray            & ba)
	:
        communicate_comp(NA, true),
        m_verbose(0), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    ParticleContainer (const Array<Geometry>            & geom, 
		       const Array<DistributionMapping> & dmap,
		       const Array<BoxArray>            & ba,
		       const Array<int>                 & rr)
	:
        communicate_comp(NA, true),
        m_verbose(0), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba,rr)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    ~ParticleContainer () {}

    void Define (ParGDBBase* gdb)
    {
      m_gdb = gdb;
        reserveData();
      resizeData();
    }

    void Define (const Geometry            & geom, 
		 const DistributionMapping & dmap,
		 const BoxArray            & ba)
  {
    m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = &m_gdb_object;
    reserveData();
    resizeData();
  }

    void Define (const Array<Geometry>            & geom, 
		 const Array<DistributionMapping> & dmap,
		 const Array<BoxArray>            & ba,
		 const Array<int>                 & rr)
    {
	m_gdb_object = ParGDB(geom, dmap, ba, rr);
	m_gdb = &m_gdb_object;
        reserveData();
	resizeData();
    }

    void SetParticleBoxArray (int lev, const BoxArray& new_ba)
	{ 
	  m_gdb->SetParticleBoxArray(lev, new_ba);
	}

    void SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
	{ m_gdb->SetParticleDistributionMap(lev, new_dmap); }

    const BoxArray& ParticleBoxArray (int lev) const 
	{ return m_gdb->ParticleBoxArray(lev); }

    const DistributionMapping& ParticleDistributionMap (int lev) const 
	{ return m_gdb->ParticleDistributionMap(lev); }

    const ParGDBBase* GetParGDB () const { return m_gdb; }

    void reserveData();
    void resizeData();

    void RedefineDummyMF (int lev);

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);
    void InitFromBinaryFile (const std::string& file, int extradata);
    void InitFromBinaryMetaFile (const std::string& file, int extradata);
    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false, RealBox bx = RealBox());
    void InitOnePerCell (Real x_off, Real y_off, Real z_off,
                         Real particleMass, const MultiFab& particle_mf);
    void InitNRandomPerCell (int n_per_cell, Real particleMass, const MultiFab& particle_mf);

    void addOneParticle (int id_in, int cpu_in, 
                         std::vector<Real>& xloc, std::vector<Real>& attributes); 
    void addNParticles (int n, D_DECL(Array<Real>& x, Array<Real>& y, Array<Real>& z),
                        int na, Array<Real>& attributes);

    void GetParticleIDs        (Array<int> & part_ids);
    void GetParticleCPU        (Array<int> & part_cpu);
    void GetParticleLocations  (Array<Real>& part_locs);
    void GetParticleData       (Array<Real>& part_data, int start_comp, int num_comp);

    void SetParticleLocations  (Array<Real>& part_locs);

    void MoveRandom ();
    void MoveRandom (int level);

    void Increment (MultiFab& mf, int level);
    long IncrementWithTotal (MultiFab& mf, int level, bool local = false);

    // rho_index: rho index in rdata
    Real sumParticleMass (int rho_index, int level, bool local = false) const;

    // Set the flag that allows particles to live near the domain boundary and throw away 
    //     the part of their contribution in AssignDensity that is outside the domain.
    void SetAllowParticlesNearBoundary(bool value);
 
    void Redistribute (bool where_already_called = false,
                       int  lev_min              = 0, 
		       int  nGrow                = 0);
    //
    // OK checks that all particles are in the right places (for some value of right)
    //
    // These flags are used to do proper checking for subcycling particles
    // the default values are fine for non-subcycling methods
    //
    bool OK (int lev_min = 0 , int ngrow = 0, int finest_level = -1) const;

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

    //
    // The Following methods are for managing Virtual and Ghost Particles.
    //
    // Removes all particles at a given level
    //
    void RemoveParticlesAtLevel (int level);
    // 
    void RemoveParticlesNotAtFinestLevel ();
    //
    // Creates virtual particles for a given level that represent
    // in some capacity all particles at finer levels
    //
    void CreateVirtualParticles (int level, AoS& virts) const;
    // 
    // Create ghost particles for a given level that are copies of particles
    // near coarse->fine boundaries in level-1
    //
    void CreateGhostParticles (int level, int ngrow, AoS& ghosts) const;
    //
    // Add particles from a pbox to the grid at this level
    //
    void AddParticlesAtLevel (int level, AoS& virts);

    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint = true) const;

    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint = true);

    void WritePlotFile (const std::string& dir, const std::string& name) const;

    void WriteAsciiFile (const std::string& file);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    const AoSMap& GetParticles(int lev) const { return m_particles[lev]; }
    AoSMap& GetParticles(int lev) { return m_particles[lev]; }

    // 
    // Functions depending the layout of the data.  Use with caution.
    //
    // Multi-level version
    void AssignDensity (int rho_index, bool sub_cycle, Array<std::unique_ptr<MultiFab> >& mf, 
			int lev_min = 0, int ncomp = 1, int finest_level = -1) const;
    // Single-level version
    void AssignDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				   int ncomp=1, int particle_lvl_offset = 0) const;
    void AssignCellDensitySingleLevelFort (int rho_index, MultiFab& mf, int level,
					   int ncomp=1, int particle_lvl_offset = 0) const;
    void AssignCellDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    void NodalDepositionSingleLevel   (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    //
    void moveKick (MultiFab& acceleration, int level, Real timestep, 
		   Real a_new = 1.0, Real a_half = 1.0,
		   int start_comp_for_accel = -1);

    IntVect Index (const Particle<NR, NI>& p, int lev) const;

    //
    // Updates a particle's location (Where), tries to periodic shift any particles
    // that have left the domain. May need work (see inline comments)
    //

    void Reset (ParticleType& prt, bool update, bool verbose=true) const; 

    //
    // Returns true if the particle was shifted.
    //
    bool PeriodicShift (ParticleType& prt) const;

    //
    // Used by AssignDensity.
    //
    bool CrseToFine (const BoxArray&       cfba, 
		     const Array<IntVect>& cells, 
		     Array<IntVect>&       cfshifts, 
		     const Geometry&       gm, 
		     Array<int>&           which, 
		     Array<IntVect>&       pshifts) const;

    bool FineToCrse (const ParticleType&                p,
		     int                                flev,
		     const Array<IntVect>&              fcells,
		     const BoxArray&                    fvalid,
		     const BoxArray&                    compfvalid_grown,
		     Array<IntVect>&                    ccells,
		     Array<Real>&                       cfracs,
		     Array<int>&                        which,
		     Array<int>&                        cgrid,
		     Array<IntVect>&                    pshifts,
		     std::vector< std::pair<int,Box> >& isects) const;

    void FineCellsToUpdateFromCrse (const ParticleType&                p,
				    int                                lev,
				    const IntVect&                     ccell,
				    const IntVect&                     cshift,
				    Array<int>&                        fgrid,
				    Array<Real>&                       ffrac,
				    Array<IntVect>&                    fcells,
				    std::vector< std::pair<int,Box> >& isects) const;

    std::vector<bool> communicate_comp;
    static bool do_tiling;
    static IntVect tile_size;

protected:

    //
    // Checks a particle's location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    //
    bool Where (const ParticleType& prt, ParticleLocData& pld, 
		int lev_min = 0, int finest_level = -1) const;

    //
    // Checks whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    //
    bool EnforcePeriodicWhere (ParticleType& prt, ParticleLocData& pld,
                               int lev_min = 0, int finest_level = -1) const;

    //
    // Checks whether a particle is within its grid (including ngrow ghost cells).
    //
    bool SingleLevelGrownWhere (const ParticleType& p, ParticleLocData& pld, int lev, int ngrow) const;

    bool OnSameGrids (int level, const MultiFab& mf) const { return m_gdb->OnSameGrids(level, mf); }


    // Helper function for Checkpoint() and WritePlotFile().
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Array<int>&    which,
                         Array<int>&    count,
                         Array<long>&   where,
                         bool           is_checkpoint) const;

    template <class RTYPE>
    void ReadParticles (int            cnt,
			int            grd,
			int            lev,
			bool           is_checkpoint,
			std::ifstream& ifs);

    //
    // The member data.
    //
    int         m_verbose;
    ParGDBBase* m_gdb;
    bool        allow_particles_near_boundary;
    ParGDB      m_gdb_object;

    //
    // The array-of-structs data for the entire AMR hierarchy
    //
    Array<AoSMap> m_particles;

    //
    // The struct-of-arrays data for the entire AMR hierarchy
    //
    Array<SoAMap> m_arraydata;

    Array<std::unique_ptr<MultiFab> > m_dummy_mf;

private:
    void AssignDensityDoit (int rho_index, 
			    Array<std::unique_ptr<MultiFab> >& mf, 
			    std::map<int, Array<ParticleCommData> >& data,
			    int ncomp, int lev_min = 0) const;

    std::pair<long,long> StartIndexInGlobalArray () const;

    void RedistributeMPI (std::map<int, Array<char> >& not_ours,
			  int lev_min = 0, int finest_level = 0, int nGrow = 0);

    void locateParticle(ParticleType& p, ParticleLocData& pld, int lev_min,
			int theEffectiveFinestLevel, int nGrow) const;

    static int getTileIndex(const IntVect& iv, const Box& box, Box& tbx);

    static void Initialize ();
};

template <int NR, int NI = 0, int NA = 0>
class ParIter
    : public MFIter
{
public:
    using ContainerType = ParticleContainer<NR,NI,NA>;
    using ParticleType  = typename ContainerType::ParticleType;
    using AoS           = typename ContainerType::AoS;
    using SoA           = typename ContainerType::SoA;

    ParIter (ContainerType& pc, int level) 
        : 
        MFIter(*pc.m_dummy_mf[level],
               pc.do_tiling ? pc.tile_size : IntVect::TheZeroVector()),
        m_aos_map(pc.m_particles[level]),
        m_soa_map(pc.m_arraydata[level]),
        m_level(level) 
        {}

    AoS& GetParticleData () {
        auto r = m_aos_map.find(std::make_pair(index(),LocalTileIndex()));
        BL_ASSERT(r != m_aos_map.end());
        return r->second;
    }

    const AoS& GetParticleData () const {
        auto r = m_aos_map.find(std::make_pair(index(),LocalTileIndex()));
        BL_ASSERT(r != m_aos_map.end());
        return r->second;
    }

    SoA& GetArrayData () {
        auto r = m_soa_map.find(std::make_pair(index(),LocalTileIndex()));
        BL_ASSERT(r != m_soa_map.end());
        return r->second;
    }

    const SoA& GetArrayData () const {
        auto r = m_soa_map.find(std::make_pair(index(),LocalTileIndex()));
        BL_ASSERT(r != m_soa_map.end());
        return r->second;
    }

private:
    typename ContainerType::AoSMap m_aos_map;
    typename ContainerType::SoAMap m_soa_map;
    int m_level;
};

#include "AMReX_ParticleI.H"
#include "AMReX_ParticleInit.H"
#include "AMReX_ParticleContainerI.H"

}

#endif /*_PARTICLES_H_*/
