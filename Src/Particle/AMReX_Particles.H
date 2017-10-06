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
#include <random>

#include <AMReX_ParmParse.H>
#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_RealBox.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_NFiles.H>

#include <AMReX_Particles_F.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if CUDA
#include "AMReX_CudaManagedAllocator.H"
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
  IntVect m_cell {AMREX_D_DECL(-1,-1,-1)};
  Box     m_gridbox;
  Box     m_tilebox;
};


//
// The struct used to store particles.
//
template<int NReal, int NInt=0>
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
    // and optionally we have NReal additional real attributes.
    //
    union rm_t
    {
      RealType pos[BL_SPACEDIM];
      RealType arr[BL_SPACEDIM+NReal];
    };

    //
    // The integer data. We always have id and cpu, and optionally we
    // have NInt additional integer attributes.
    //
    union im_t
    {
      struct {int id; int cpu;};
      int arr[2+NInt];
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

    static void Interp (const Particle<NReal, NInt>& prt, const Geometry& geom, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Geometry& geom, const Particle<NReal, NInt>& p, Real* grav);

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
    static void CIC_Cells_Fracs_Basic (const Particle<NReal, NInt>& p,
                                       const Real*                  plo,
                                       const Real*                  dx,
                                       Real*                        fracs,
                                       IntVect*                     cells);

    //
    // Wraps the arbitrary dx function.
    //
    static int CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                const Real*                  plo,
                                const Real*                  dx,
                                Vector<Real>&                 fracs,
                                Vector<IntVect>&              cells);
    //
    // Does CIC computations for arbitrary particle/grid dx's.
    //
    static int CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                const Real*                  plo,
                                const Real*                  dx_geom,
                                const Real*                  dx_part,
                                Vector<Real>&                 fracs,
                                Vector<IntVect>&              cells);
};

template <int NReal, int NInt>
class ArrayOfStructs {
public:
    using ParticleType  = Particle<NReal, NInt>;
    using RealType      = typename ParticleType::RealType;
    using Iterator      = typename Vector<ParticleType>::iterator;
    using ConstIterator = typename Vector<ParticleType>::const_iterator;

#if CUDA
    using ParticleVector = Vector<ParticleType, CudaManagedAllocator<ParticleType> >;
#else
    using ParticleVector = Vector<ParticleType>;
#endif

    static constexpr int SizeInReal = sizeof(ParticleType) / sizeof(RealType);

    const ParticleVector& operator()() const { return m_data; }
          ParticleVector& operator()()       { return m_data; }

    std::size_t size () const { return m_data.size(); }
    int numParticles () const { return m_data.size(); }

    bool empty () const { return m_data.empty(); }

    const RealType* data () const { return &(m_data[0].m_rdata.arr[0]); }
    RealType*       data ()       { return &(m_data[0].m_rdata.arr[0]); }

    std::pair<int,int> dataShape () const {
	return std::make_pair(SizeInReal, static_cast<int>(m_data.size()));
    }

    void push_back (const ParticleType& p) { return m_data.push_back(p); }
    void pop_back() {m_data.pop_back(); }
    bool empty() {return m_data.empty(); }

    const ParticleType& back() const {return m_data.back(); }
    ParticleType      & back()       {return m_data.back(); }

    const ParticleType& operator[] (int i) const { return m_data[i]; }
    ParticleType      & operator[] (int i)       { return m_data[i]; }
    
    Iterator erase( ConstIterator first, ConstIterator second) { return m_data.erase(first, second); }

    template< class InputIt >
    void insert( Iterator pos, InputIt first, InputIt last ) { m_data.insert(pos, first, last); }

    typename Vector<ParticleType>::iterator begin () { return m_data.begin(); }
    typename Vector<ParticleType>::const_iterator begin () const { return m_data.begin(); }
    typename Vector<ParticleType>::const_iterator cbegin () const { return m_data.cbegin(); }
        
    typename ParticleVector::iterator end () { return m_data.end(); }
    typename ParticleVector::const_iterator end () const { return m_data.end(); }
    typename ParticleVector::const_iterator cend () const { return m_data.cend(); }

private:
    ParticleVector m_data;
};
template <int NReal, int NInt> constexpr int ArrayOfStructs<NReal, NInt>::SizeInReal;


template <int NReal, int NInt>
struct StructOfArrays {

#if CUDA
    using RealVector = Vector<Real, CudaManagedAllocator<Real>>;
    using  IntVector = Vector<int , CudaManagedAllocator<int >>;
#else
    using RealVector = Vector<Real>;
    using  IntVector = Vector<int >;
#endif

    std::array<RealVector, NReal>& GetRealData () { return m_rdata; }
    std::array< IntVector,   NInt>& GetIntData  () { return m_idata; }

    const std::array<RealVector, NReal>& GetRealData () const { return m_rdata; }
    const std::array< IntVector,  NInt>& GetIntData  () const { return m_idata; }

    RealVector& GetRealData (const int index) {
        BL_ASSERT(index < NReal);
        return m_rdata[index];
    }
    
    const RealVector& GetRealData (const int index) const {
        BL_ASSERT(index < NReal);
        return m_rdata[index];
    }

    IntVector& GetIntData (const int index) {
        BL_ASSERT(index < NInt);
        return m_idata[index];
    }
    
    const IntVector& GetIntData (const int index) const {
        BL_ASSERT(index < NInt);
        return m_idata[index];
    }

    int numParticles () const { 
        if (NReal > 0)
            return m_rdata[0].size();
        else if (NInt > 0) 
            return m_idata[0].size();
        else
            return 0;
    }

private:
    std::array<RealVector, NReal> m_rdata;
    std::array< IntVector,  NInt> m_idata;    
};


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTile
{
    using ParticleType = Particle<NStructReal, NStructInt>;
    using AoS = ArrayOfStructs<NStructReal, NStructInt>;
    using SoA = StructOfArrays<NArrayReal, NArrayInt>;

    AoS&       GetArrayOfStructs ()       { return m_aos_tile; }
    const AoS& GetArrayOfStructs () const { return m_aos_tile; }

    SoA&       GetStructOfArrays ()       { return m_soa_tile; }
    const SoA& GetStructOfArrays () const { return m_soa_tile; }

    bool empty () const { return m_aos_tile.empty(); }
    
    std::size_t size () const { return m_aos_tile.size(); }

    int numParticles () const { return m_aos_tile.numParticles(); }

    ///
    /// Add one particle to this tile.
    ///
    void push_back (const ParticleType& p) { m_aos_tile().push_back(p); }

    ///
    /// Add a Real value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_real (int comp, Real v) { 
        m_soa_tile.GetRealData(comp).push_back(v);
    }

    ///
    /// Add Real values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_real (const std::array<Real, NArrayReal>& v) { 
        for (int i = 0; i < NArrayReal; ++i) {
            m_soa_tile.GetRealData(i).push_back(v[i]);
        }
    }

    ///
    /// Add a range of Real values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, const Real* beg, const Real* end) {
        auto it = m_soa_tile.GetRealData(comp).end();
        m_soa_tile.GetRealData(comp).insert(it, beg, end);
    }

    ///
    /// Add npar copies of the Real value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, std::size_t npar, Real v) {
        auto new_size = m_soa_tile.GetRealData(comp).size() + npar;
        m_soa_tile.GetRealData(comp).resize(new_size, v);
    }

    ///
    /// Add an int value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_int (int comp, int v) { 
        m_soa_tile.GetIntData(comp).push_back(v);
    }
    
    ///
    /// Add int values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_int (const std::array<int, NArrayInt>& v) { 
        for (int i = 0; i < NArrayInt; ++i) {
            m_soa_tile.GetIntData(i).push_back(v[i]);
        }
    }
    
    ///
    /// Add a range of int values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, const int* beg, const int* end) {
        auto it = m_soa_tile.GetIntData(comp).end();
        m_soa_tile.GetIntData(comp).insert(it, beg, end);
    }
    
    ///
    /// Add npar copies of the int value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, std::size_t npar, int v) {
        auto new_size = m_soa_tile.GetIntData(comp).size() + npar;
        m_soa_tile.GetIntData(comp).resize(new_size, v);
    }

private:

    AoS m_aos_tile;
    SoA m_soa_tile;
};

///
/// This struct is used to pass initial data into the various Init methods
/// of the particle container. That data should be initialized in the order
/// real struct data, int struct data, real array data, int array data. 
/// If fewer components are specified than the template parameters specify for,
/// a given component, then the extra values will be set to zero. If more components 
/// are specified, it is a compile-time error.
///
/// Example usage:
///
///     ParticleInitType<0, 2, 4, 1> pdata = {{}, {7, 9}, {1.5, 2.5, 3.5, 4.5}, {11}}; 
///
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleInitType
    {
        std::array<double, NStructReal> real_struct_data;
        std::array<int,    NStructInt > int_struct_data;
        std::array<double, NArrayReal > real_array_data;
        std::array<int,    NArrayInt  > int_array_data;
    };

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt> class ParIterBase;

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class ParticleContainer
{
    friend class ParIterBase<true,NStructReal, NStructInt, NArrayReal, NArrayInt>;
    friend class ParIterBase<false,NStructReal, NStructInt, NArrayReal, NArrayInt>;

public:
    //
    // The type of Particles we hold.
    //
    using ParticleType = Particle<NStructReal, NStructInt>;
    using RealType = typename Particle<NStructReal, NStructInt>::RealType;

    using ParticleTileType = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleInitData = ParticleInitType<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    // A single level worth of particles is indexed (grid id, tile id)
    // for both SoA and AoS data.
    using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    using AoS = typename ParticleTileType::AoS;
    using SoA = typename ParticleTileType::SoA;

    using RealArray     = typename SoA::RealArray;
    using  IntArray     = typename SoA::IntArray;
    using ParticleArray = typename AoS::ParticleArray;


    ParticleContainer ()
      : 
      communicate_real_comp(NArrayReal, true),
      communicate_int_comp(NArrayInt, true),
      m_verbose(0), 
      m_gdb(nullptr),
      allow_particles_near_boundary(false)
    {
        Initialize ();
    }

    ParticleContainer (ParGDBBase* gdb)
        :
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
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
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
        m_verbose(0), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    ParticleContainer (const Vector<Geometry>            & geom, 
		       const Vector<DistributionMapping> & dmap,
		       const Vector<BoxArray>            & ba,
		       const Vector<int>                 & rr)
	:
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
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

    void Define (const Vector<Geometry>            & geom, 
		 const Vector<DistributionMapping> & dmap,
		 const Vector<BoxArray>            & ba,
		 const Vector<int>                 & rr)
    {
	m_gdb_object = ParGDB(geom, dmap, ba, rr);
	m_gdb = &m_gdb_object;
        reserveData();
	resizeData();
    }

    void SetParticleBoxArray (int lev, const BoxArray& new_ba)
	{ m_gdb->SetParticleBoxArray(lev, new_ba); }

    void SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
	{ m_gdb->SetParticleDistributionMap(lev, new_dmap); }

    const BoxArray& ParticleBoxArray (int lev) const 
	{ return m_gdb->ParticleBoxArray(lev); }

    const DistributionMapping& ParticleDistributionMap (int lev) const 
	{ return m_gdb->ParticleDistributionMap(lev); }

    const Geometry& Geom (int lev) const { return m_gdb->Geom(lev); }

    int finestLevel () const { return m_gdb->finestLevel(); }
    int maxLevel ()    const { return m_gdb->maxLevel(); }
    
    int numLocalTilesAtLevel (int lev) const { return m_particles[lev].size(); }

    const ParGDBBase* GetParGDB () const { return m_gdb; }

    void reserveData();
    void resizeData();

    void RedefineDummyMF (int lev);

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);
    void InitFromBinaryFile (const std::string& file, int extradata);
    void InitFromBinaryMetaFile (const std::string& file, int extradata);
    
    ///
    /// This initializes the particle container with icount randomly distributed 
    /// particles. If serialize is true, then the particles will all be generated
    /// on the IO Process, and the particle positions will be broadcast to all 
    /// other process. If serialize is false, then the particle positions will be
    /// randomly generated in parallel, which each process using the random seed
    /// iseed + MyProc. The particles can be constrained to lie within the RealBox
    /// bx, if so desired. The default is the full domain.
    ///
    void InitRandom (long icount, unsigned long iseed, 
                     const ParticleInitData& mass,
                     bool serialize = false, RealBox bx = RealBox());

    ///
    /// This initializes the container with icount randomly distributed particles
    /// per box, using the random seed iseed. All the particles have the same data
    /// and attributes, which are passed using the pdata struct.
    ///
    /// This routine is used when we want to replicate a box for a scaling study -- 
    /// within each box the distribution is random but the particle data is replicated 
    /// across all boxes in the container. The boxes are assumed to be those 
    /// on the coarsest level.
    ///
    void InitRandomPerBox (long icount, unsigned long iseed, const ParticleInitData& pdata);
    
    ///
    /// This initializes the particle container with one particle per cell, 
    /// where the other particle data and attributes are all contant. The 
    /// coarsest level is used to generate the particle positions. The particle
    /// variable values are passed in through the pdata struct. The parameters
    /// x_off, y_off, and z_off represent offsets between 0 and 1 that show
    /// where inside the cells to place the particles. 0.5 means cell centered.
    ///
    void InitOnePerCell (Real x_off, Real y_off, Real z_off,
                         const ParticleInitData& pdata);

    ///
    /// This initializes the particle container with n_per_cell randomly
    /// distributed particles per cell, where the other particle data and
    /// and attributes are all constant. The cells on the coarsest level 
    /// are used to generate the particle positions. The particle variable
    /// values are passed in through the pdata struct.
    ///
    void InitNRandomPerCell (int n_per_cell, const ParticleInitData& pdata);

    void GetParticleIDs        (Vector<int> & part_ids);
    void GetParticleCPU        (Vector<int> & part_cpu);
    void GetParticleLocations  (Vector<Real>& part_locs);
    void GetParticleData       (Vector<Real>& part_data, int start_comp, int num_comp);
    void GetArrayData          (Vector<Real>& part_data, int start_comp, int num_comp);

    void SetParticleLocations  (Vector<Real>& part_locs);

    void MoveRandom ();
    void MoveRandom (int level);

    void Increment (MultiFab& mf, int level);
    long IncrementWithTotal (MultiFab& mf, int level, bool local = false);

    // rho_index: rho index in rdata
    Real sumParticleMass (int rho_index, int level, bool local = false) const;

    // Set the flag that allows particles to live near the domain boundary and throw away 
    //     the part of their contribution in AssignDensity that is outside the domain.
    void SetAllowParticlesNearBoundary(bool value);
 
    void Redistribute (int lev_min = 0, int lev_max = -1, int nGrow = 0);
    //
    // OK checks that all particles are in the right places (for some value of right)
    //
    // These flags are used to do proper checking for subcycling particles
    // the default values are fine for non-subcycling methods
    //
    bool OK (int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

    void ByteSpread () const;
    //
    // Returns # of particles at specified the level.
    //
    // If "only_valid" is true it only counts valid particles.
    //
    long NumberOfParticlesAtLevel (int level, bool only_valid = true, bool only_local = false) const;
    Vector<long> NumberOfParticlesInGrid  (int level, bool only_valid = true, bool only_local = false) const;
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
    void AddParticlesAtLevel (AoS& particles, int level, int nGrow=0);

    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint = true,
                     const Vector<std::string>& real_comp_names = Vector<std::string>(),
                     const Vector<std::string>&  int_comp_names = Vector<std::string>()) const;

    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint = true);

    void WritePlotFile (const std::string& dir, const std::string& name, 
                        const Vector<std::string>& real_comp_names = Vector<std::string>(),
                        const Vector<std::string>&  int_comp_names = Vector<std::string>()) const;

    void WriteAsciiFile (const std::string& file);

    void WriteCoarsenedAsciiFile (const std::string& filename);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    const Vector<ParticleLevel>& GetParticles () const { return m_particles; }
    Vector      <ParticleLevel>& GetParticles ()       { return m_particles; }

    const ParticleLevel& GetParticles (int lev) const { return m_particles[lev]; }
    ParticleLevel      & GetParticles (int lev)       { return m_particles[lev]; }

    const ParticleTileType& ParticlesAt (int lev, int grid, int tile) const 
        { return m_particles[lev].at(std::make_pair(grid,tile)); }
    ParticleTileType&       ParticlesAt (int lev, int grid, int tile)
        { return m_particles[lev].at(std::make_pair(grid, tile)); }

    // 
    // Functions depending the layout of the data.  Use with caution.
    //
    // Multi-level version
    void AssignDensity (int rho_index, bool sub_cycle, Vector<std::unique_ptr<MultiFab> >& mf, 
			int lev_min = 0, int ncomp = 1, int finest_level = -1) const;
    // Single-level version
    void AssignDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				   int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignDensityFort (int rho_index,
                            Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
                            int lev_min, int ncomp, int finest_level) const;

    void InterpolateFort (Vector<std::unique_ptr<MultiFab> >& mesh_data, 
                          int lev_min, int lev_max);

    void InterpolateSingleLevelFort (MultiFab& mesh_data, int lev);

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

    IntVect Index (const Particle<NStructReal, NStructInt>& p, int lev) const;

    //
    // Updates a particle's location (Where), tries to periodic shift any particles
    // that have left the domain. May need work (see inline comments)
    //

    ParticleLocData Reset (ParticleType& prt, bool update, bool verbose=true) const; 

    //
    // Returns true if the particle was shifted.
    //
    bool PeriodicShift (ParticleType& prt) const;

    //
    // Used by AssignDensity.
    //
    bool CrseToFine (const BoxArray&       cfba, 
		     const Vector<IntVect>& cells, 
		     Vector<IntVect>&       cfshifts, 
		     const Geometry&       gm, 
		     Vector<int>&           which, 
		     Vector<IntVect>&       pshifts) const;

    bool FineToCrse (const ParticleType&                p,
		     int                                flev,
		     const Vector<IntVect>&              fcells,
		     const BoxArray&                    fvalid,
		     const BoxArray&                    compfvalid_grown,
		     Vector<IntVect>&                    ccells,
		     Vector<Real>&                       cfracs,
		     Vector<int>&                        which,
		     Vector<int>&                        cgrid,
		     Vector<IntVect>&                    pshifts,
		     std::vector< std::pair<int,Box> >& isects) const;

    void FineCellsToUpdateFromCrse (const ParticleType&                p,
				    int                                lev,
				    const IntVect&                     ccell,
				    const IntVect&                     cshift,
				    Vector<int>&                        fgrid,
				    Vector<Real>&                       ffrac,
				    Vector<IntVect>&                    fcells,
				    std::vector< std::pair<int,Box> >& isects) const;

    MFIter MakeMFIter (int lev) const {
        BL_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], do_tiling ? tile_size : IntVect::TheZeroVector());
    }

    // struct-of-array stuff
    std::vector<bool> communicate_real_comp;
    std::vector<bool> communicate_int_comp;
    int num_real_soa_comps, num_int_soa_comps;

    static bool do_tiling;
    static IntVect tile_size;
    
protected:

    //
    // Checks a particle's location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    // Only if lev_min == lev_max, nGrow can be > 0 (i.e., including
    // nGrow ghost cells).
    //
    bool Where (const ParticleType& prt, ParticleLocData& pld, 
		int lev_min = 0, int lev_max = -1, int nGrow=0) const;

    //
    // Checks whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    //
    bool EnforcePeriodicWhere (ParticleType& prt, ParticleLocData& pld,
                               int lev_min = 0, int lev_max = -1) const;

    bool OnSameGrids (int level, const MultiFab& mf) const { return m_gdb->OnSameGrids(level, mf); }


    // Helper function for Checkpoint() and WritePlotFile().
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Vector<int>&    which,
                         Vector<int>&    count,
                         Vector<long>&   where,
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
    
private:
    void AssignDensityDoit (int rho_index, 
			    Vector<std::unique_ptr<MultiFab> >& mf, 
			    std::map<int, Vector<ParticleCommData> >& data,
			    int ncomp, int lev_min = 0) const;

    std::pair<long,long> StartIndexInGlobalArray () const;

    void RedistributeMPI (std::map<int, Vector<char> >& not_ours,
			  int lev_min = 0, int lev_max = 0, int nGrow = 0);

    void locateParticle(ParticleType& p, ParticleLocData& pld,
                        int lev_min, int lev_max, int nGrow) const;

    static int getTileIndex(const IntVect& iv, const Box& box, Box& tbx);

    void Initialize ();

    size_t particle_size, superparticle_size;
    int num_real_comm_comps, num_int_comm_comps;
    Vector<ParticleLevel> m_particles;
    Vector<std::unique_ptr<MultiFab> > m_dummy_mf;
};


template <bool is_const, int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class ParIterBase
    : public MFIter
{
private:
    using PCType = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ContainerRef    = typename std::conditional<is_const, PCType const&, PCType&>::type;
    using ParticleTileRef = typename std::conditional
        <is_const, typename PCType::ParticleTileType const&, typename PCType::ParticleTileType &>::type;
    using ParticleTilePtr = typename std::conditional
        <is_const, typename PCType::ParticleTileType const*, typename PCType::ParticleTileType *>::type;
    using AoSRef          = typename std::conditional
        <is_const, typename PCType::AoS const&, typename PCType::AoS&>::type;
    using SoARef          = typename std::conditional
        <is_const, typename PCType::SoA const&, typename PCType::SoA&>::type;

public:
    ParIterBase (ContainerRef pc, int level);

    void operator++ () {
        ++m_pariter_index;
        currentIndex = m_valid_index[m_pariter_index];
    }
    
    ParticleTileRef GetParticleTile () const { return *m_particle_tiles[m_pariter_index]; }

    AoSRef GetArrayOfStructs () const { return GetParticleTile().GetArrayOfStructs(); }

    SoARef GetStructOfArrays () const { return GetParticleTile().GetStructOfArrays(); }


    void GetPosition (AMREX_D_DECL(Vector<Real>& x,
                                   Vector<Real>& y,
                                   Vector<Real>& z)) const;

    int numParticles () const { return GetArrayOfStructs().numParticles(); }    
protected:
    int m_level;
    int m_pariter_index;
    Vector<int> m_valid_index;
    Vector<ParticleTilePtr> m_particle_tiles;
};


template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class ParIter
    : public ParIterBase<false,NStructReal,NStructInt, NArrayReal, NArrayInt>
{
public:
    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;

    ParIter (ContainerType& pc, int level)
        : ParIterBase<false,NStructReal,NStructInt, NArrayReal, NArrayInt>(pc,level)
        {}

    void SetPosition (AMREX_D_DECL(const Vector<Real>& x,
                                   const Vector<Real>& y,
                                   const Vector<Real>& z)) const;
};

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class ParConstIter
    : public ParIterBase<true,NStructReal,NStructInt, NArrayReal, NArrayInt>
{
public:
    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;

    ParConstIter (ContainerType const& pc, int level)
        : ParIterBase<true,NStructReal,NStructInt, NArrayReal, NArrayInt>(pc,level)
        {}
};


#include "AMReX_ParticleI.H"
#include "AMReX_ParticleInit.H"
#include "AMReX_ParticleContainerI.H"
#include "AMReX_ParIterI.H"

}

#endif /*_PARTICLES_H_*/
