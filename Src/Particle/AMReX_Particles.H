#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

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

#include <AMReX_ParmParse.H>
#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Particles_F.H>
#include <AMReX_RealBox.H>

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

template<int NR, int NI> class PartIter;

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
  int     m_lev;
  int     m_grid;
  int     m_tile;
  IntVect m_cell;
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

    static IntVect Index (const Particle<NR, NI>& p, const Geometry& geom);

    //
    // Checks a particle's location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    //
    static bool Where (const Particle<NR, NI>& prt, const ParGDBBase* gdb, 
		       ParticleLocData& pld, int lev_min = 0, int finest_level = -1);

    //
    // Checks whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    //
    static bool PeriodicWhere (Particle<NR, NI>& prt, const ParGDBBase* gdb, 
			       ParticleLocData& pld, int lev_min = 0, int finest_level = -1);

    //
    // Checks whether a particle is within its grid (including ngrow ghost cells).
    //
    static bool RestrictedWhere (const Particle<NR, NI>& p, const ParGDBBase* gdb, 
				 ParticleLocData& pld, int ngrow);

    //
    // Checks a particle's location on a specific level.
    // (Yes this is distict from the functionality provided above)
    //
    static bool SingleLevelWhere (const Particle<NR, NI>& p, const ParGDBBase* gdb, 
				  ParticleLocData& pld, int level);

    //
    // Updates a particle's location (Where), tries to periodic shift any particles
    // that have left the domain. May need work (see inline comments)
    //

    static void Reset (Particle<NR, NI>& prt, const ParGDBBase* gdb, bool update, bool verbose=true); 

    //
    // Returns true if the particle was shifted.
    //
    static bool PeriodicShift (Particle<NR, NI>& prt, const ParGDBBase* gdb);

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

    //
    // Used by AssignDensity.
    //
    static bool CrseToFine (const BoxArray&       cfba, 
                            const Array<IntVect>& cells, 
                            Array<IntVect>&       cfshifts, 
                            const Geometry&       gm, 
                            Array<int>&           which, 
                            Array<IntVect>&       pshifts);

    static bool FineToCrse (const Particle<NR, NI>&            p, 
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

    static void FineCellsToUpdateFromCrse (const Particle<NR, NI>&            p, 
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

template <int NR, int NI = 0>
class ParticleContainer
{
public:
    //
    // The type of Particles we hold.
    //
    typedef Particle<NR, NI> ParticleType;

    // This defines the array-of-structs datatype.
    using AoS = Array<ParticleType>;

    // The struct-of-arrays data is stored as a double array where the first index is the
    // component number and the second is the particle index.
    using SoA = Array<Array<Real> >;

    // A single level worth of particles is indexed first by grid id, 
    // then by tile id, for both SoA and AoS data.
    using AoSMap = std::map<int, Array<AoS> >;
    using SoAMap = std::map<int, Array<SoA> >;

    ParticleContainer ()
      : m_verbose(1), m_gdb(nullptr), allow_particles_near_boundary(false), m_npartdata(0) 
  {
    communicate_comp.resize(m_npartdata, true);
  }

    ParticleContainer (ParGDBBase* gdb)
        :
      m_verbose(1), m_gdb(gdb), allow_particles_near_boundary(false), m_npartdata(0) 
  {
    communicate_comp.resize(m_npartdata, true);
  }

    ParticleContainer (const Geometry            & geom, 
		       const DistributionMapping & dmap,
		       const BoxArray            & ba,
		       int                         npartdata = 0)
	:
        m_verbose(1), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba),
	m_npartdata(npartdata)
    {
	m_gdb = & m_gdb_object;
	communicate_comp.resize(m_npartdata, true);
    }

    ParticleContainer (const Array<Geometry>            & geom, 
		       const Array<DistributionMapping> & dmap,
		       const Array<BoxArray>            & ba,
		       const Array<int>                 & rr,
		       int                                npartdata = 0)
	:
        m_verbose(1), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba,rr),
	m_npartdata(npartdata)
    {
      m_gdb = & m_gdb_object;
      communicate_comp.resize(m_npartdata, true);
    }

    ~ParticleContainer () {}

    void Define (ParGDBBase* gdb)
    {
      m_gdb = gdb;
      m_npartdata = 0;
      communicate_comp.resize(m_npartdata, true);
    }

    void Define (const Geometry            & geom, 
		 const DistributionMapping & dmap,
		 const BoxArray            & ba,
		 int                         npartdata = 0)
  {
    m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = &m_gdb_object;
    m_npartdata = npartdata;
    communicate_comp.resize(m_npartdata, true);
  }

    void Define (const Array<Geometry>            & geom, 
		 const Array<DistributionMapping> & dmap,
		 const Array<BoxArray>            & ba,
		 const Array<int>                 & rr,
		 int                                npartdata = 0)
    {
	m_gdb_object = ParGDB(geom, dmap, ba, rr);
	m_gdb = &m_gdb_object;
	m_npartdata = npartdata;
	communicate_comp.resize(m_npartdata, true);
    }

    void SetParticleBoxArray (int lev, const BoxArray& new_ba)
	{ m_gdb->SetParticleBoxArray(lev, new_ba); }

    void SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
	{ m_gdb->SetParticleDistributionMap(lev, new_dmap); }

    const BoxArray& ParticleBoxArray (int lev) const 
	{ return m_gdb->ParticleBoxArray(lev); }

    const DistributionMapping& ParticleDistributionMap (int lev) const 
	{ return m_gdb->ParticleDistributionMap(lev); }

    const ParGDBBase* GetParGDB () const { return m_gdb; }

    void allocateData();

    void locateParticle(ParticleType& p, ParticleLocData& pld, int lev_min,
			int theEffectiveFinestLevel, int nGrow);

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);
    void InitFromBinaryFile (const std::string& file, int extradata);
    void InitFromBinaryMetaFile (const std::string& file, int extradata);
    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false, RealBox bx = RealBox());
    void InitOnePerCell (Real x_off, Real y_off, Real z_off,
                         Real particleMass, MultiFab& particle_mf);
    void InitNRandomPerCell (int n_per_cell, Real particleMass, MultiFab& particle_mf);

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
    long IncrementWithTotal (MultiFab& mf, int level);

    // rho_index: rho index in rdata
    Real sumParticleMass (int rho_index, int level) const;

    // Set the flag that allows particles to live near the domain boundary and throw away 
    //     the part of their contribution in AssignDensity that is outside the domain.
    void SetAllowParticlesNearBoundary(bool value);
 
    void Redistribute (int  lev_min              = 0, 
		       int  nGrow                = 0);

  void RedistributeMPI (std::map<int, AoS>& not_ours,
			std::map<int, SoA>& comm_pdata);

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
    void AssignCellDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    void NodalDepositionSingleLevel   (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    //
    void moveKick (MultiFab& acceleration, int level, Real timestep, 
		   Real a_new = 1.0, Real a_half = 1.0,
		   int start_comp_for_accel = -1);

    std::vector<bool> communicate_comp;
    static bool do_tiling;
    static IntVect tile_size;

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
    Array<SoAMap> m_partdata;

    int m_npartdata;

    using PIter = PartIter<NR, NI>;

private:
    void AssignDensityDoit (int rho_index, 
			    Array<std::unique_ptr<MultiFab> >& mf, 
			    std::map<int, std::vector<ParticleCommData> >& data,
			    int ncomp, int lev_min = 0) const;
};


template<int NR, int NI=0>
class PartIter
{

public:
  PartIter (const ParticleContainer<NR, NI>& _mypc, int _level)
    : mypc(_mypc), 
      pmap(_mypc.GetParticles(_level)),
      box_array(_mypc.ParticleBoxArray(_level)),
      level(_level),
      do_tiling(_mypc.do_tiling), 
      tile_size(_mypc.tile_size)
  {
    Initialize();
  }
  
  bool isValid () {return current_index < num_boxes;}

  void operator++ () {
    ++current_tile_index;
    if (current_tile_index >= num_tiles[current_index]) {
      ++current_index;
      current_tile_index = 0;
    }
  }

  int gridIndex () const {return grid_ids[current_index];}

  int tileIndex () const {return current_tile_index;}

  int numTiles(int grid_id) {
    for (int i = 0; i < num_boxes; i++) {
      if (grid_ids[i] == grid_id)
	return num_tiles[i];
    }
    return -1;
  }

  Box gridBox () {return box_array[grid_ids[current_index]];}

  Box tileBox () {
    if (do_tiling) {
      int grid_id = grid_ids[current_index];
      const Box& box = box_array[grid_id];
      const IntVect& small = box.smallEnd();
      const IntVect& size  = box.size();
      IntVect ntiles(D_DECL((size[0]+tile_size[0]-1)/tile_size[0],
			    (size[1]+tile_size[1]-1)/tile_size[1],
			    (size[2]+tile_size[2]-1)/tile_size[2]));
      
      IntVect index = { D_DECL(0, 0, 0) };
      int linear_index = current_tile_index;
      for (int i = 0; i < BL_SPACEDIM; i++) {
	index[i] = linear_index % ntiles[i];
	linear_index /= ntiles[i];
      }
      
      IntVect tile_lo_corner = small + index*tile_size; 
      IntVect tile_hi_corner = tile_lo_corner + tile_size - IntVect::TheUnitVector(); 
      Box tile_box(tile_lo_corner, tile_hi_corner);
      tile_box &= box;
      return tile_box;
    }
    
    return box_array[grid_ids[current_index]];
  }

private:

  void Initialize()
  {
    current_index = 0;
    current_tile_index = 0;
    num_boxes = pmap.size();
    num_tiles.resize(num_boxes, 0);
    grid_ids.resize(num_boxes, 0);
    int index = 0;
    for (auto it = pmap.begin(); it != pmap.end(); ++it, ++index) {
      grid_ids[index] = it->first;
      if (do_tiling) {
	const Box& box = box_array[it->first];
	const IntVect& small = box.smallEnd();
	const IntVect& size  = box.size();
	IntVect ntiles(D_DECL((size[0]+tile_size[0]-1)/tile_size[0],
			      (size[1]+tile_size[1]-1)/tile_size[1],
			      (size[2]+tile_size[2]-1)/tile_size[2]));
	int ntottiles = D_TERM(ntiles[0],*ntiles[1],*ntiles[2]);
	num_tiles[index] = ntottiles;
      }
      else {
	num_tiles[index] = 1;
      }
    }
  }

  const ParticleContainer<NR, NI>& mypc;
  const typename ParticleContainer<NR, NI>::AoSMap& pmap;
  const BoxArray& box_array;

  int level;
  bool do_tiling;
  IntVect tile_size;

  int current_index;
  int current_tile_index;
  int num_boxes;
  Array<int> grid_ids;
  Array<int> num_tiles;
};


#include "AMReX_ParticleI.H"
#include "AMReX_ParticleInit.H"
#include "AMReX_ParticleContainerI.H"

}

#endif /*_PARTICLES_H_*/
