#ifndef AMREX_PARTICLES_H_
#define AMREX_PARTICLES_H_

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
#include <AMReX_NFiles.H>
#include <AMReX_VectorIO.H>
#include <AMReX_Particle_mod_K.H>
#include <AMReX_ParticleMPIUtil.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_Functors.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_ParticleCommunication.H>
#include <AMReX_ParticleLocator.H>
#include <AMReX_Scan.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef AMREX_USE_CUDA
#include <thrust/device_vector.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/retag.h>
#include <thrust/device_ptr.h>
#include <thrust/tuple.h>
#include <thrust/gather.h>
#endif

namespace amrex {

/**
* \brief A struct used for communicating particle data accross processes
*  during multi-level operations.
*/
struct ParticleCommData
{
  typedef ParticleReal RealType;
  int     m_lev;
  int     m_grid;
  IntVect m_cell;
  RealType m_data[1 + AMREX_SPACEDIM];
};

/**
* \brief A struct used for storing a particle's position in the AMR hierarchy.
*/
struct ParticleLocData
{
  int     m_lev  = -1;
  int     m_grid = -1;
  int     m_tile = -1;
  IntVect m_cell {AMREX_D_DECL(-1,-1,-1)};
  Box     m_gridbox;
  Box     m_tilebox;
  Box     m_grown_gridbox;
};

/**
* /brief A struct used to pass initial data into the various Init methods.
* This struct is used to pass initial data into the various Init methods
* of the particle container. That data should be initialized in the order
* real struct data, int struct data, real array data, int array data.
* If fewer components are specified than the template parameters specify for,
* a given component, then the extra values will be set to zero. If more components
* are specified, it is a compile-time error.
*
* Example usage:
*
*     ParticleInitType<0, 2, 4, 1> pdata = {{}, {7, 9}, {1.5, 2.5, 3.5, 4.5}, {11}};
*/
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleInitType
{
    std::array<double, NStructReal> real_struct_data;
    std::array<int,    NStructInt > int_struct_data;
    std::array<double, NArrayReal > real_array_data;
    std::array<int,    NArrayInt  > int_array_data;
};

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
class ParIterBase;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
class ParIter;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
class ParConstIter;
    
/**
 * \brief A distributed container for Particles sorted onto the levels, grids,
 * and tiles of a block-structured AMR hierarachy.
 *
 * The data layout on a single tile is determined by the value of the following
 * template parameters:
 *
 * \tparam NStructReal The number of extra Real components in the particle struct
 * \tparam NStructInt The number of extra integer components in the particle struct
 * \tparam NArrayReal The number of extra Real components stored in struct-of-array form
 * \tparam NArrayInt The number of extra integer components stored in struct-of-array form
 *
 */
template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class ParticleContainer
{
    friend class ParIterBase<true,NStructReal, NStructInt, NArrayReal, NArrayInt>;
    friend class ParIterBase<false,NStructReal, NStructInt, NArrayReal, NArrayInt>;

public:
    //! \brief The type of Particles we hold.
    using ParticleType = Particle<NStructReal, NStructInt>;
    //! \brief The type of the "SuperParticle" which stored all components in AoS form
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;
    //! \brief The type of the Real data 
    using RealType = typename Particle<NStructReal, NStructInt>::RealType;

#ifdef BL_SINGLE_PRECISION_PARTICLES
    RealDescriptor ParticleRealDescriptor = FPC::Native32RealDescriptor();
#else
    RealDescriptor ParticleRealDescriptor = FPC::Native64RealDescriptor();
#endif

    using ParticleContainerType = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleTileType = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleInitData = ParticleInitType<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    //! A single level worth of particles is indexed (grid id, tile id)
    //! for both SoA and AoS data.
    using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    using AoS = typename ParticleTileType::AoS;
    using SoA = typename ParticleTileType::SoA;

    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;
    using ParticleVector   = typename AoS::ParticleVector;
    using CharVector       = Gpu::ManagedDeviceVector<char>;
    using SendBuffer       = Gpu::PolymorphicDeviceVector<char>;
    using ParIterType      = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParConstIterType = ParConstIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    //! \brief Default constructor - construct an empty particle container that has no concept
    //!  of a level hierarchy. Must be properly initialized later.
    ParticleContainer ()
      :
      communicate_real_comp(NArrayReal, true),
      communicate_int_comp(NArrayInt, true),
      m_verbose(0),
      m_gdb(nullptr),
      m_runtime_comps_defined(false),
      m_num_runtime_real(0),
      m_num_runtime_int(0)        
    {
        Initialize ();
    }

    //! \brief Construct a particle container using a ParGDB object. The container will
    //! track changes in the grid structure of the ParGDB automatically.
    //!
    //! \param gdb A pointer to a ParGDBBase, which contains pointers to the Geometry,
    //! DistributionMapping, and BoxArray objects that define the AMR hierarchy. Usually,
    //! this is generated by an AmrCore or AmrLevel object.
    //! 
    ParticleContainer (ParGDBBase* gdb)
        :
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
	m_verbose(0),
	m_gdb(gdb),
        m_runtime_comps_defined(false),
        m_num_runtime_real(0),
        m_num_runtime_int(0)
    {
        Initialize ();
        reserveData();
        resizeData();
    }

    //! \brief Construct a particle container using a given Geometry, DistributionMapping,
    //! and BoxArray. Single level version.
    //!
    //! \param the Geometry object, which describes the problem domain
    //! \param A DistributionMapping, which describes how the boxes are distributed onto MPI tasks
    //! \param A BoxArray, which gives the set of grid boxes
    //!     
    ParticleContainer (const Geometry            & geom,
		       const DistributionMapping & dmap,
		       const BoxArray            & ba)
	:
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
        m_verbose(0),
	m_gdb_object(geom,dmap,ba),
        m_runtime_comps_defined(false),
        m_num_runtime_real(0),
        m_num_runtime_int(0)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    //! \brief Construct a particle container using a given Geometry, DistributionMapping,
    //! BoxArray and Vector of refinement ratios. Multi-level version.
    //!
    //! \param geom A Vector of Geometry objects, one for each level
    //! \param dmap A Vector of DistributionMappings, one for each level
    //! \param ba A Vector of BoxArrays, one for each level
    //! \param rr A Vector of integer refinement ratios, of size num_levels - 1. rr[n] gives the
    //! refinement ratio between levels n and n+1
    //!
    ParticleContainer (const Vector<Geometry>            & geom,
		       const Vector<DistributionMapping> & dmap,
		       const Vector<BoxArray>            & ba,
		       const Vector<int>                 & rr)
	:
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
        m_verbose(0),
	m_gdb_object(geom,dmap,ba,rr),
        m_runtime_comps_defined(false),
        m_num_runtime_real(0),
        m_num_runtime_int(0)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    //! \brief Same as the above, but accepts different refinement ratios in each direction.
    //!
    //! \param geom A Vector of Geometry objects, one for each level
    //! \param dmap A Vector of DistributionMappings, one for each level
    //! \param ba A Vector of BoxArrays, one for each level
    //! \param rr A Vector of IntVect refinement ratios, of size num_levels - 1. rr[n] gives the
    //! refinement ratio between levels n and n+1
    //!
    ParticleContainer (const Vector<Geometry>            & geom,
		       const Vector<DistributionMapping> & dmap,
		       const Vector<BoxArray>            & ba,
		       const Vector<IntVect>             & rr)
	:
        communicate_real_comp(NArrayReal, true),
        communicate_int_comp(NArrayInt, true),
        m_verbose(0),
	m_gdb_object(geom,dmap,ba, [&]() -> Vector<int> {
                Vector<int> ref_ratio;
                for (int i = 0; i < static_cast<int>(rr.size()); ++i)
                {
#if AMREX_SPACEDIM > 1
                    AMREX_ASSERT(rr[i][0] == rr[i][1]);
#endif
#if AMREX_SPACEDIM > 2
                    AMREX_ASSERT(rr[i][0] == rr[i][2]);
#endif
                    ref_ratio.push_back(rr[i][0]);
                }
                return ref_ratio;  }() ),
        m_runtime_comps_defined(false),
        m_num_runtime_real(0),
        m_num_runtime_int(0)
    {
        Initialize ();
        m_gdb = & m_gdb_object;
        reserveData();
        resizeData();
    }

    virtual ~ParticleContainer () {}

    ParticleContainer ( const ParticleContainer &) = delete;
    ParticleContainer& operator= ( const ParticleContainer & ) = delete; 
    
    ParticleContainer ( ParticleContainer && ) = default;
    ParticleContainer& operator= ( ParticleContainer && ) = default;


    //! \brief Define a default-constructed ParticleContainer using a ParGDB object.
    //! The container will track changes in the grid structure of the ParGDB automatically.
    //!
    //! \param gdb A pointer to a ParGDBBase, which contains pointers to the Geometry,
    //! DistributionMapping, and BoxArray objects that define the AMR hierarchy. Usually,
    //! this is generated by an AmrCore or AmrLevel object.
    //! 
    void Define (ParGDBBase* gdb)
    {
        m_gdb = gdb;
        reserveData();
        resizeData();
    }

    //! \brief Define a default-constructed ParticleContainer using a ParGDB object.
    //! Single-level version.
    //!
    //! \param the Geometry object, which describes the problem domain
    //! \param A DistributionMapping, which describes how the boxes are distributed onto MPI tasks
    //! \param A BoxArray, which gives the set of grid boxes
    //!     
    void Define (const Geometry            & geom,
		 const DistributionMapping & dmap,
		 const BoxArray            & ba)
    {
        m_gdb_object = ParGDB(geom, dmap, ba);
        m_gdb = &m_gdb_object;
        reserveData();
        resizeData();
    }

    //! \brief Define a default-constructed ParticleContainer using a ParGDB object.
    //! Multi-level version
    //!
    //! \param geom A Vector of Geometry objects, one for each level
    //! \param dmap A Vector of DistributionMappings, one for each level
    //! \param ba A Vector of BoxArrays, one for each level
    //! \param rr A Vector of integer refinement ratios, of size num_levels - 1. rr[n] gives the
    //! refinement ratio between levels n and n+1
    //!
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

    //! \brief Set the particle BoxArray. If the container was previously set to
    //! to track the AMR hierarchy of an AmrCore or AmrLevel object, that correspondence
    //! will be broken here.
    //!
    //! \param lev The level on which to set the BoxArray.
    //! \param new_ba The new BoxArray to use.
    //!
    void SetParticleBoxArray (int lev, const BoxArray& new_ba)
	{ m_gdb->SetParticleBoxArray(lev, new_ba); }

    //! \brief Set the particle DistributionMapping. If the container was previously set to
    //! to track the AMR hierarchy of an AmrCore or AmrLevel object, that correspondence
    //! will be broken here.
    //!
    //! \param lev The level on which to set the dm.
    //! \param new_dmap The new DistributionMapping to use.
    //!
    void SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
	{ m_gdb->SetParticleDistributionMap(lev, new_dmap); }

    //! \brief Get the BoxArray for a given level
    //!
    //! \param lev The level.
    //!
    const BoxArray& ParticleBoxArray (int lev) const
	{ return m_gdb->ParticleBoxArray(lev); }

    //! \brief Get the DistributionMapping for a given level
    //!
    //! \param lev The level.
    //!
    const DistributionMapping& ParticleDistributionMap (int lev) const
	{ return m_gdb->ParticleDistributionMap(lev); }

    //! \brief Get the Geometry for a given level
    //!
    //! \param lev The level.
    //!
    const Geometry& Geom (int lev) const { return m_gdb->Geom(lev); }

    //! \brief the finest level actually defined for the ParticleContainer
    int finestLevel () const { return m_gdb->finestLevel(); }

    //! \brief the finest allowed level in the ParticleContainer, whether it is defined or not.
    int maxLevel ()    const { return m_gdb->maxLevel(); }

    //! \brief the number of defined levels in the ParticleContainer    
    int numLevels()    const { return finestLevel() + 1; }

    //! \brief The total number of tiles on this rank on this level     
    int numLocalTilesAtLevel (int lev) const { return m_particles[lev].size(); }

    //! \brief Get the ParGDB object used to define this container (const version)    
    const ParGDBBase* GetParGDB () const { return m_gdb; }

    //! \brief Get the ParGDB object used to define this container
          ParGDBBase* GetParGDB ()       { return m_gdb; }
    
    void reserveData ();
    void resizeData ();

    void RedefineDummyMF (int lev);

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);

    void InitFromBinaryFile (const std::string& file, int extradata);

    void InitFromBinaryMetaFile (const std::string& file, int extradata);

    /**
    * \brief 
    * This initializes the particle container with icount randomly distributed
    * particles. If serialize is true, then the particles will all be generated
    * on the IO Process, and the particle positions will be broadcast to all
    * other process. If serialize is false, then the particle positions will be
    * randomly generated in parallel, which each process using the random seed
    * iseed + MyProc. The particles can be constrained to lie within the RealBox
    * bx, if so desired. The default is the full domain.
    *
    * \param icount
    * \param iseed
    * \param mass
    * \param serialize
    * \param bx
    */
    void InitRandom (long icount, unsigned long iseed,
                     const ParticleInitData& mass,
                     bool serialize = false, RealBox bx = RealBox());


    /**
    * \brief
    * This initializes the container with icount randomly distributed particles
    * per box, using the random seed iseed. All the particles have the same data
    * and attributes, which are passed using the pdata struct.
    * 
    * This routine is used when we want to replicate a box for a scaling study --
    * within each box the distribution is random but the particle data is replicated
    * across all boxes in the container. The boxes are assumed to be those
    * on the coarsest level.
    *
    * \param icount
    * \param iseed
    * \param pdata
    */
    void InitRandomPerBox (long icount, unsigned long iseed, const ParticleInitData& pdata);


    /**
    * \brief
    * This initializes the particle container with one particle per cell,
    * where the other particle data and attributes are all contant. The
    * coarsest level is used to generate the particle positions. The particle
    * variable values are passed in through the pdata struct. The parameters
    * x_off, y_off, and z_off represent offsets between 0 and 1 that show
    * where inside the cells to place the particles. 0.5 means cell centered.
    *
    * \param x_off
    * \param y_off
    * \param z_off
    * \param pdata
    */
    void InitOnePerCell (Real x_off, Real y_off, Real z_off,
                         const ParticleInitData& pdata);


    /**
    * \brief
    * This initializes the particle container with n_per_cell randomly
    * distributed particles per cell, where the other particle data and
    * and attributes are all constant. The cells on the coarsest level
    * are used to generate the particle positions. The particle variable
    * values are passed in through the pdata struct.
    *
    * \param n_per_cell
    * \param pdata
    */
    void InitNRandomPerCell (int n_per_cell, const ParticleInitData& pdata);

    void MoveRandom ();

    void MoveRandom (int level);

    void Increment (MultiFab& mf, int level);

    long IncrementWithTotal (MultiFab& mf, int level, bool local = false);

    /**
    * \brief rho_index: rho index in rdata
    *
    * \param rho_index
    * \param level
    * \param local
    */
    Real sumParticleMass (int rho_index, int level, bool local = false) const;

    /**
    * \brief Redistribute puts all the particles back in the right places (for some value of right)
    *
    * Assigns particles to the levels, grids, and tiles that contain their current positions. 
    * If periodic boundary conditions are used, those will be enforced here.
    *
    * If Redistribute is called with default arguments, all particles will be placed on the
    * finest level that covers their current positions. 
    * 
    * The lev_min, lev_max, and nGrow flags are used to do proper checking for subcycling particles.
    * The default values are fine for non-subcycling methods
    *
    * The local flag controls whether this is `local` or `global` Redistribute. 
    * In a local Redistribute, particles can only have moved a certain distance since the last
    * time Redistribute() was called. Thus, communication only needs to happen between neighboring
    * ranks. In a global Redistribute, the particles can potentially go from any rank to any rank.
    * This usually happens after initialiation or when doing dynamic load balancing. 
    *
    * \param lev_min
    * \param lev_max
    * \param nGrow
    * \param local
    */
    void Redistribute (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0);

    /**
     * \brief Sort the particles on each tile by cell, using Fortran ordering.
     */
    void SortParticlesByCell();

    /**
    * \brief OK checks that all particles are in the right places (for some value of right)
    *
    * These flags are used to do proper checking for subcycling particles
    * the default values are fine for non-subcycling methods
    *
    * \param lev_min
    * \param lev_max
    * \param nGrow
    */
    bool OK (int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

    void ByteSpread () const;

    /**
    * \brief Returns # of particles at specified the level.
    *
    * If "only_valid" is true it only counts valid particles.
    *
    * \param level
    * \param only_valid
    * \param only_local
    */
    long NumberOfParticlesAtLevel (int level, bool only_valid = true, bool only_local = false) const;

    Vector<long> NumberOfParticlesInGrid  (int level, bool only_valid = true, bool only_local = false) const;

    /**
    * \brief Returns # of particles at all levels
    *
    * If "only_valid" is true it only counts valid particles.
    *
    * \param only_valid
    * \param only_local
    */
    long TotalNumberOfParticles (bool only_valid=true, bool only_local=false) const;


    /**
    * \brief The Following methods are for managing Virtual and Ghost Particles.
    *
    * Removes all particles at a given level
    *
    * \param level
    */
    void RemoveParticlesAtLevel (int level);

    void RemoveParticlesNotAtFinestLevel ();

    /**
    * \brief Creates virtual particles for a given level that represent
    * in some capacity all particles at finer levels
    *
    * \param level
    * \param virts
    */
    void CreateVirtualParticles (int level, AoS& virts) const;

    /**
    * \brief Create ghost particles for a given level that are copies of particles
    * near coarse-\>fine boundaries in level-1
    *
    * \param level
    * \param ngrow
    * \param ghosts
    */
    void CreateGhostParticles (int level, int ngrow, AoS& ghosts) const;

    /**
    * \brief Add particles from a pbox to the grid at this level
    *
    * \param particles
    * \param level
    * \param nGrow
    */
    void AddParticlesAtLevel (AoS& particles, int level, int nGrow=0);


    /**
    * \brief Clear all the particles in this container. This does not free memory.
    */
    void clearParticles();


    /**
    * \brief Copy particles from other to this ParticleContainer. Will clear all the
    * particles from this container first. local controls whether or not to call
    * Redistribute() after adding the particles.
    *
    * \param other
    */
    void copyParticles(const ParticleContainerType& other, bool local=false);

    /**
    * \brief Add particles from other to this ParticleContainer. local controls 
    * whether or not to call Redistribute after adding the particles.
    *
    * \param other
    */
    void addParticles(const ParticleContainerType& other, bool local=false);

    /**
    * \brief Write a contiguous chunk of real particle data to an ostream.
    * 
    * \param data A pointer to the start of the buffer to write
    * \param size The number of elements to write
    * \param os The ostream into which to write the data
    * \param rd A RealDescriptor describing the type of the floating point data
    */    
    void WriteParticleRealData (void* data, size_t size,
                                std::ostream& os, const RealDescriptor& rd) const;
    
    /**
    * \brief Read a contiguous chunk of real particle data from an istream.
    * 
    * \param data A pointer to the start of the buffer into which to read
    * \param size The number of elements to read
    * \param os The istream from which to read the data
    * \param rd A RealDescriptor describing the type of the floating point data
    */    
    void ReadParticleRealData (void* data, size_t size,
                               std::istream& is, const RealDescriptor& rd);

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void Checkpoint (const std::string& dir, const std::string& name) const;

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     */
    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint,
                     const Vector<std::string>& real_comp_names = Vector<std::string>(),
                     const Vector<std::string>& int_comp_names = Vector<std::string>()) const;

     /**
      * \brief Writes particle data to disk in the AMReX native format. 
      * 
      * \param dir The base directory into which to write (i.e. "plt00000")
      * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
      * \param write_real_comp for each real component, whether or not we include that component in the file
      * \param write_int_comp for each integer component, whether or not we include that component in the file
      * \param real_comp_names for each real component, a name to label the data with
      * \param int_comp_names for each integer component, a name to label the data with      
     */
    void WriteBinaryParticleData (const std::string& dir,
                                  const std::string& name,
                                  const Vector<int>& write_real_comp,
                                  const Vector<int>& write_int_comp,    
                                  const Vector<std::string>& real_comp_names,
                                  const Vector<std::string>&  int_comp_names) const;
    
    void CheckpointPre ();

    void CheckpointPost ();

    /**
     *   \brief Restart from checkpoint
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void Restart (const std::string& dir, const std::string& file);

    /**
     *   \brief Older version, for backwards compatibility
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param is_checkpoint Whether the particle id and cpu are included in the file.
     */
    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint);
    
    /**
     *  \brief This version of WritePlotFile writes all components and assigns component names
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void WritePlotFile (const std::string& dir, const std::string& name) const;

    /**
     *  \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components. 
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with      
     */
    void WritePlotFile (const std::string& dir, const std::string& name,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>&  int_comp_names) const;

    /**
     *  This version of WritePlotFile writes all components and allows the user to specify
     *  the names of the components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     */
    void WritePlotFile (const std::string& dir, const std::string& name,
                        const Vector<std::string>& real_comp_names) const;
    
    /**
     *  This version of WritePlotFile assigns component names, but allows the user to pass
     *  in a vector of ints that toggle on / off the writing of specfic components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     */
    void WritePlotFile (const std::string& dir,
                        const std::string& name,
                        const Vector<int>& write_real_comp,
                        const Vector<int>& write_int_comp) const;
    
    /**
     * \brief This is the most general version of WritePlotFile, which takes component 
     * names and flags for whether to write each variable as components. Note that 
     * the user should pass in vectors containing names of all the components, whether
     * they are written or not.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     */
    void WritePlotFile (const std::string& dir,
                        const std::string& name,
                        const Vector<int>& write_real_comp,
                        const Vector<int>& write_int_comp,    
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>&  int_comp_names) const;
    
    void WritePlotFilePre ();

    void WritePlotFilePost ();

    void WriteAsciiFile (const std::string& file);

    void WriteCoarsenedAsciiFile (const std::string& filename);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    const Vector<ParticleLevel>& GetParticles () const { return m_particles; }
    Vector      <ParticleLevel>& GetParticles ()       { return m_particles; }

    const ParticleLevel& GetParticles (int lev) const { return m_particles[lev]; }

    ParticleLevel      & GetParticles (int lev)       { return m_particles[lev]; }

    const ParticleTileType& ParticlesAt (int lev, int grid, int tile) const
    { return m_particles[lev].at(std::make_pair(grid, tile)); }

    ParticleTileType&       ParticlesAt (int lev, int grid, int tile)
    { return m_particles[lev].at(std::make_pair(grid, tile)); }

    template <class Iterator>

    const ParticleTileType& ParticlesAt (int lev, const Iterator& iter) const
        { return ParticlesAt(lev, iter.index(), iter.LocalTileIndex()); }

    template <class Iterator>

    ParticleTileType&       ParticlesAt (int lev, const Iterator& iter)
        { return ParticlesAt(lev, iter.index(), iter.LocalTileIndex()); }

    ParticleTileType& DefineAndReturnParticleTile (int lev, int grid, int tile)
    {
        m_particles[lev][std::make_pair(grid, tile)].define(NumRuntimeRealComps(), NumRuntimeIntComps());
        return ParticlesAt(lev, grid, tile);
    }

    /**
    * \brief Functions depending the layout of the data.  Use with caution.
    *
    * \param rho_index
    * \param mf_to_be_filled
    * \param lev_min
    * \param ncomp
    * \param finest_level
    * \param ngrow
    */
    void AssignDensity (int rho_index,
                        Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled,
                        int lev_min, int ncomp, int finest_level, int ngrow=2) const;

    void Interpolate (Vector<std::unique_ptr<MultiFab> >& mesh_data,
                      int lev_min, int lev_max);

    void InterpolateSingleLevel (MultiFab& mesh_data, int lev);

    void AssignCellDensitySingleLevel (int rho_index, MultiFab& mf, int level,
                                       int ncomp=1, int particle_lvl_offset = 0) const;

    void moveKick (MultiFab& acceleration, int level, Real timestep,
		   Real a_new = 1.0, Real a_half = 1.0,
		   int start_comp_for_accel = -1);

    IntVect Index (const Particle<NStructReal, NStructInt>& p, int lev) const;


    /**
    * \brief Updates a particle's location (Where), tries to periodic shift any particles
    * that have left the domain. May need work (see inline comments)
    *
    * \param prt
    * \param update
    * \param verbose
    * \param pld
    */
    ParticleLocData Reset (ParticleType& prt, bool update, bool verbose=true,
                           ParticleLocData pld = ParticleLocData()) const;


    /**
    * \brief Returns true if the particle was shifted.
    *
    * \param prt
    */
    bool PeriodicShift (ParticleType& prt) const;

    MFIter MakeMFIter (int lev, const MFItInfo& info) const {
        BL_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], info);
    }

    MFIter MakeMFIter (int lev) const {
        BL_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], do_tiling ? tile_size : IntVect::TheZeroVector());
    }

    MFIter MakeMFIter (int lev, bool tile) const {
        BL_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], tile ? tile_size : IntVect::TheZeroVector());
    }


    //! struct-of-array stuff
    std::vector<bool> communicate_real_comp;
    std::vector<bool> communicate_int_comp;

    static bool do_tiling;
    static IntVect tile_size;

    void SetLevelDirectoriesCreated(bool tf) {
      levelDirectoriesCreated = tf;
    }

    void SetUsePrePost(bool tf) {
      usePrePost = tf;
    }
    bool GetUsePrePost() {
      return usePrePost;
    }

    void SetUseUnlink(bool tf) {
      doUnlink = tf;
    }

    bool GetUseUnlink() {
      return doUnlink;
    }

    void RedistributeCPU (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0);

    void RedistributeGPU (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0);

    bool OKCPU (int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

    bool OKGPU (int lev_min = 0, int lev_max = -1, int nGrow = 0) const;

#ifdef AMREX_USE_CUDA    
    void EnforcePeriodicGPU ();
#endif

    template <typename T,
              typename std::enable_if<std::is_same<T,bool>::value,int>::type=0>
    void AddRealComp (T communicate=true)
    {
        m_runtime_comps_defined = true;
        m_num_runtime_real++;
        communicate_real_comp.push_back(communicate);
        SetParticleSize();
    }

    template <typename T,
              typename std::enable_if<std::is_same<T,bool>::value,int>::type=0>
    void AddIntComp (T communicate=true)
    {
        m_runtime_comps_defined = true;
        m_num_runtime_int++;
        communicate_int_comp.push_back(communicate);
        SetParticleSize();
    }

    const ParticleBufferMap& BufferMap () const {return m_buffer_map;} 

    Vector<int> NeighborProcs(int ngrow) const         
    {
        return computeNeighborProcs(this->GetParGDB(), ngrow);
    };

    bool OnSameGrids (int level, const MultiFab& mf) const { return m_gdb->OnSameGrids(level, mf); }

protected:

    mutable amrex::Vector<int> neighbor_procs;

    /**
    * \brief Checks a particle's location on levels lev_min and higher.
    * Returns false if the particle does not exist on that level.
    * Only if lev_min == lev_max, nGrow can be \> 0 (i.e., including
    * nGrow ghost cells).
    *
    * \param prt
    * \param pld
    * \param lev_min
    * \param lev_max
    * \param nGrow
    * \param local_grid
    */
    bool Where (const ParticleType& prt, ParticleLocData& pld,
		int lev_min = 0, int lev_max = -1, int nGrow=0, int local_grid=-1) const;


    /**
    * \brief Checks whether the particle has crossed a periodic boundary in such a way
    * that it is on levels lev_min and higher.
    *
    * \param prt
    * \param pld
    * \param lev_min
    * \param lev_max
    * \param local_grid
    */
    bool EnforcePeriodicWhere (ParticleType& prt, ParticleLocData& pld,
                               int lev_min = 0, int lev_max = -1, int local_grid=-1) const;

    /**
    * \brief Helper function for Checkpoint() and WritePlotFile().
    *
    * \param level
    * \param ofs
    * \param fnum
    * \param which
    * \param count
    * \param where
    */
    void WriteParticles (int level, std::ofstream& ofs, int fnum,
                         Vector<int>& which, Vector<int>& count, Vector<long>& where,
                         const Vector<int>& write_real_comp, const Vector<int>& write_int_comp) const;

    template <class RTYPE>
    void ReadParticles (int cnt, int grd, int lev, std::ifstream& ifs, int finest_level_in_file);
    
    void SetParticleSize ();

    void BuildRedistributeMask(int lev, int nghost=1) const;
    mutable std::unique_ptr<iMultiFab> redistribute_mask_ptr;
    mutable int redistribute_mask_nghost = std::numeric_limits<int>::min();

    void defineBufferMap () const;
    mutable ParticleBufferMap m_buffer_map;

    //! The member data.
    int         m_verbose;
    ParGDBBase* m_gdb;
    ParGDB      m_gdb_object;


    //! ---- variables for i/o optimization saved for pre and post checkpoint
    bool         levelDirectoriesCreated;
    bool         usePrePost;
    bool         doUnlink;
    int maxnextidPrePost;
    mutable int nOutFilesPrePost;
    long nparticlesPrePost;
    Vector<long> nParticlesAtLevelPrePost;          //!< ---- [level]
    mutable Vector<Vector<int>>  whichPrePost;      //!< ---- [level]
    mutable Vector<Vector<int>>  countPrePost;      //!< ---- [level]
    mutable Vector<Vector<long>> wherePrePost;      //!< ---- [level]
    mutable std::string HdrFileNamePrePost;
    mutable Vector<std::string> filePrefixPrePost;

#ifdef AMREX_USE_CUDA
    ParticleVector m_aos_r;
    std::array<RealVector, NArrayReal> m_rdata_r;
    std::array<IntVector, NArrayInt>  m_idata_r;
    Gpu::ManagedDeviceVector<int> m_grids_r;
    Gpu::ManagedDeviceVector<int> m_levs_r;

    Gpu::ManagedDeviceVector<int> m_lo;
    Gpu::ManagedDeviceVector<int> m_hi;
    Gpu::ManagedDeviceVector<int> m_output;

    Gpu::ManagedDeviceVector<int> m_grids_tmp;
    Gpu::ManagedDeviceVector<int> m_levs_tmp;

    ParticleVector m_aos_to_redistribute;
    std::vector<RealVector> m_real_arrays_to_redistribute;
    std::vector<IntVector > m_int_arrays_to_redistribute;
    Gpu::ManagedDeviceVector<int> m_grids_to_redistribute;

    std::map<int, SendBuffer> m_not_ours;

    ParticleCopyOp redistribute_copy_op;
    ParticleCopyPlan redistribute_copy_plan;
    Gpu::DeviceVector<SuperParticleType> redistribute_snd_buffer;
    Gpu::DeviceVector<SuperParticleType> redistribute_rcv_buffer;

    Gpu::PinnedDeviceVector<SuperParticleType> pinned_snd_buffer;
    Gpu::PinnedDeviceVector<SuperParticleType> pinned_rcv_buffer;

    mutable AmrParticleLocator m_particle_locator;
    
#endif

    int NumRuntimeRealComps () const { return m_num_runtime_real; }
    int NumRuntimeIntComps  () const { return m_num_runtime_int;  } 
    
    int NumRealComps () const { return NArrayReal + NumRuntimeRealComps(); }
    int NumIntComps  () const { return NArrayInt  + NumRuntimeIntComps() ; } 

    
private:

    virtual void particlePostLocate(ParticleType& p, const ParticleLocData& pld,
                                    const int lev) {};

    virtual void correctCellVectors(int old_index, int new_index,
				    int grid, const ParticleType& p) {};

    void RedistributeMPI (std::map<int, Vector<char> >& not_ours,
			  int lev_min = 0, int lev_max = 0, int nGrow = 0, int local=0);

    void locateParticle(ParticleType& p, ParticleLocData& pld,
                        int lev_min, int lev_max, int nGrow, int local_grid=-1) const;

    void Initialize ();

    bool m_runtime_comps_defined;
    int m_num_runtime_real;
    int m_num_runtime_int;
    
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

    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;
    using ParticleVector   = typename ContainerType::ParticleVector;

    ParIterBase (ContainerRef pc, int level);

    ParIterBase (ContainerRef pc, int level, MFItInfo& info);

#ifdef _OPENMP
    void operator++ () {

        if (dynamic) {
#pragma omp atomic capture
            m_pariter_index = nextDynamicIndex++;
        } else {
            ++m_pariter_index;
        }
        currentIndex = m_valid_index[m_pariter_index];
    }
#else
    void operator++ () {
        ++m_pariter_index;
        currentIndex = m_valid_index[m_pariter_index];
#ifdef AMREX_USE_GPU

        Gpu::Device::setStreamIndex(currentIndex);
#endif
    }
#endif

    ParticleTileRef GetParticleTile () const { return *m_particle_tiles[m_pariter_index]; }

    AoSRef GetArrayOfStructs () const { return GetParticleTile().GetArrayOfStructs(); }

    SoARef GetStructOfArrays () const { return GetParticleTile().GetStructOfArrays(); }

    template <typename Container>

    void GetPosition (AMREX_D_DECL(Container& x,
                                   Container& y,
                                   Container& z)) const;

    int numParticles () const { return GetArrayOfStructs().numParticles(); }

    int GetLevel () const { return m_level; }

    std::pair<int, int> GetPairIndex () const { return std::make_pair(this->index(), this->LocalTileIndex()); }
    
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
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParIter (ContainerType& pc, int level)
        : ParIterBase<false,NStructReal,NStructInt, NArrayReal, NArrayInt>(pc,level)
        {}

    ParIter (ContainerType& pc, int level, MFItInfo& info)
        : ParIterBase<false,NStructReal,NStructInt,NArrayReal,NArrayInt>(pc,level,info)
        {}

    template <typename Container>

    void SetPosition (AMREX_D_DECL(const Container& x,
                                   const Container& y,
                                   const Container& z)) const;
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
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParConstIter (ContainerType const& pc, int level)
        : ParIterBase<true,NStructReal,NStructInt, NArrayReal, NArrayInt>(pc,level)
        {}

    ParConstIter (ContainerType const& pc, int level, MFItInfo& info)
        : ParIterBase<true,NStructReal,NStructInt,NArrayReal,NArrayInt>(pc,level,info)
        {}

};

#include "AMReX_ParticleInit.H"
#include "AMReX_ParticleContainerI.H"
#include "AMReX_ParIterI.H"
#include "AMReX_ParticleIO.H"

}

#endif /*_PARTICLES_H_*/
