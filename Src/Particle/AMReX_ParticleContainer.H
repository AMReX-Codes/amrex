#ifndef AMREX_PARTICLE_CONTAINER_H_
#define AMREX_PARTICLE_CONTAINER_H_
#include <AMReX_Config.H>

#include <AMReX_ParticleContainerBase.H>
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
#include <AMReX_GpuContainers.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_ParticleCommunication.H>
#include <AMReX_ParticleLocator.H>
#include <AMReX_Scan.H>
#include <AMReX_DenseBins.H>
#include <AMReX_SparseBins.H>
#include <AMReX_ParticleTransformation.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParIter.H>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#ifdef AMREX_USE_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace amrex {

#ifdef AMREX_USE_HDF5_ASYNC
    extern hid_t es_par_g;
#endif
/**
* \brief A struct used for communicating particle data across processes
*  during multi-level operations.
*/
struct ParticleCommData
{
    using RealType = ParticleReal;
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
* \brief A struct used to pass initial data into the various Init methods.
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

template <bool is_const, typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
class ParIterBase_impl;

/**
 * \brief A distributed container for Particles sorted onto the levels, grids,
 * and tiles of a block-structured AMR hierarchy.
 *
 * The data layout on a single tile is determined by the value of the following
 * template parameters:
 *
 * \tparam T_NStructReal The number of extra Real components in the particle struct
 * \tparam T_NStructInt The number of extra integer components in the particle struct
 * \tparam T_NArrayReal The number of extra Real components stored in struct-of-array form
 * \tparam T_NArrayInt The number of extra integer components stored in struct-of-array form
 *
 */
template <typename T_ParticleType, int T_NArrayReal, int T_NArrayInt,
          template<class> class Allocator=DefaultAllocator,
          class T_CellAssignor=DefaultAssignor>

class ParticleContainer_impl : public ParticleContainerBase
{
public:
    using ParticleType = T_ParticleType;
    using ConstParticleType = typename ParticleType::ConstType;
    using CellAssignor = T_CellAssignor;

    //! \brief Number of extra Real components in the particle struct
    static constexpr int NStructReal = ParticleType::NReal;
    //! \brief Number of extra integer components in the particle struct
    static constexpr int NStructInt = ParticleType::NInt;
    //! \brief Number of extra Real components stored in struct-of-array form
    static constexpr int NArrayReal = T_NArrayReal;
    //! \brief Number of extra integer components stored in struct-of-array form
    static constexpr int NArrayInt = T_NArrayInt;
        //! \brief The type of the "Particle"

private:
    friend class ParIterBase_impl<true,ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    friend class ParIterBase_impl<false,ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;

public:
    //! \brief The memory allocator in use.
    template <typename T>
    using AllocatorType = Allocator<T>;
    //! \brief The type of the "SuperParticle" which stored all components in AoS form
    using SuperParticleType = Particle<NStructReal+NArrayReal, NStructInt+NArrayInt>;
    //! \brief The type of the Real data
    using RealType = typename Particle<NStructReal, NStructInt>::RealType;

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
    RealDescriptor ParticleRealDescriptor = FPC::Native32RealDescriptor();
#else
    RealDescriptor ParticleRealDescriptor = FPC::Native64RealDescriptor();
#endif

    using ParticleContainerType = ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ParticleTileType = ParticleTile<ParticleType, NArrayReal, NArrayInt, Allocator>;
    using ParticleInitData = ParticleInitType<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    //! A single level worth of particles is indexed (grid id, tile id)
    //! for both SoA and AoS data.
    using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    using AoS = typename ParticleTileType::AoS;
    using SoA = typename ParticleTileType::SoA;

    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;
    using ParticleVector   = typename AoS::ParticleVector;
    using CharVector       = Gpu::DeviceVector<char>;
    using ParIterType      = ParIter_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ParConstIterType = ParConstIter_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;

    //! \brief Default constructor - construct an empty particle container that has no concept
    //!  of a level hierarchy. Must be properly initialized later.
    ParticleContainer_impl ()
      :
      ParticleContainerBase()
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
    ParticleContainer_impl (ParGDBBase* gdb)
        :
        ParticleContainerBase(gdb)
    {
        Initialize ();
        ParticleContainer_impl::reserveData();
        ParticleContainer_impl::resizeData();
    }

    //! \brief Construct a particle container using a given Geometry, DistributionMapping,
    //! and BoxArray. Single level version.
    //!
    //! \param the Geometry object, which describes the problem domain
    //! \param A DistributionMapping, which describes how the boxes are distributed onto MPI tasks
    //! \param A BoxArray, which gives the set of grid boxes
    //!
    ParticleContainer_impl (const Geometry            & geom,
                            const DistributionMapping & dmap,
                            const BoxArray            & ba)
        :
        ParticleContainerBase(geom, dmap, ba)
    {
        Initialize ();
        ParticleContainer_impl::reserveData();
        ParticleContainer_impl::resizeData();
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
    ParticleContainer_impl (const Vector<Geometry>            & geom,
                            const Vector<DistributionMapping> & dmap,
                            const Vector<BoxArray>            & ba,
                            const Vector<int>                 & rr)
        :
        ParticleContainerBase(geom, dmap, ba, rr)
    {
        Initialize ();
        ParticleContainer_impl::reserveData();
        ParticleContainer_impl::resizeData();
    }

    //! \brief Same as the above, but accepts different refinement ratios in each direction.
    //!
    //! \param geom A Vector of Geometry objects, one for each level
    //! \param dmap A Vector of DistributionMappings, one for each level
    //! \param ba A Vector of BoxArrays, one for each level
    //! \param rr A Vector of IntVect refinement ratios, of size num_levels - 1. rr[n] gives the
    //! refinement ratio between levels n and n+1
    //!
    ParticleContainer_impl (const Vector<Geometry>            & geom,
                            const Vector<DistributionMapping> & dmap,
                            const Vector<BoxArray>            & ba,
                            const Vector<IntVect>             & rr)
        :
        ParticleContainerBase(geom, dmap, ba, rr)
    {
        Initialize ();
        ParticleContainer_impl::reserveData();
        ParticleContainer_impl::resizeData();
    }

    ~ParticleContainer_impl () override = default;

    ParticleContainer_impl ( const ParticleContainer_impl &) = delete;
    ParticleContainer_impl& operator= ( const ParticleContainer_impl & ) = delete;

    ParticleContainer_impl ( ParticleContainer_impl && ) noexcept = default;
    ParticleContainer_impl& operator= ( ParticleContainer_impl && ) noexcept = default;


    //! \brief Define a default-constructed ParticleContainer using a ParGDB object.
    //! The container will track changes in the grid structure of the ParGDB automatically.
    //!
    //! \param gdb A pointer to a ParGDBBase, which contains pointers to the Geometry,
    //! DistributionMapping, and BoxArray objects that define the AMR hierarchy. Usually,
    //! this is generated by an AmrCore or AmrLevel object.
    //!
    void Define (ParGDBBase* gdb)
    {
        this->ParticleContainerBase::Define(gdb);
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
        this->ParticleContainerBase::Define(geom, dmap, ba);
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
        this->ParticleContainerBase::Define(geom, dmap, ba, rr);
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
                 const Vector<IntVect>             & rr)
    {
        this->ParticleContainerBase::Define(geom, dmap, ba, rr);
        reserveData();
        resizeData();
    }

    //! \brief The total number of tiles on this rank on this level
    int numLocalTilesAtLevel (int lev) const {
        return (lev < m_particles.size()) ? m_particles[lev].size() : 0;
    }

    /**
     * \brief This reserves data in the vector of dummy MultiFabs used by
     * the ParticleContainer for the maximum number of levels possible.
     */
    void reserveData () override;

    /**
     * \brief This resizes the vector of dummy MultiFabs used by the
     * ParticleContainer for the current number of levels and calls
     * RedefineDummyMF on each level. Note that this must be done prior
     * to using ParticleIterator.
     */
    void resizeData () override;

    void InitFromAsciiFile (const std::string& file, int extradata,
                            const IntVect* Nrep = nullptr);

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
    void InitRandom (Long icount, ULong iseed,
                     const ParticleInitData& pdata,
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
    void InitRandomPerBox (Long icount, ULong iseed, const ParticleInitData& pdata);


    /**
    * \brief
    * This initializes the particle container with one particle per cell,
    * where the other particle data and attributes are all constant. The
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

    void Increment (MultiFab& mf, int level);

    Long IncrementWithTotal (MultiFab& mf, int level, bool local = false);

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
    * This usually happens after initialization or when doing dynamic load balancing.
    *
    * \param lev_min The minimum level consider. Particles on levels less than this will not be
    *                touched, and particles on finer levels will not be assigned to levels less
    *                than this, either.
    *                Default: 0.
    * \param lev_max The maximum level consider. Particles on levels greater than this will not be
    *                touched, and particles on coarser levels will not be assigned to levels greater
    *                than this, either. If negative, will use the finest level in the hierarchy.
    *                Default: -1.
    * \param nGrow If particles are within nGrow cells of their current box, they will not moved.
    *              This is useful for subcycling methods, when fine level particles need to be
    *              redistributed but are not necessarily at the same time as those on the coarse
    *              level.
    *              Default: 0
    * \param local If 0, this will be a non-local redistribute, meaning that particle can potentially
    *              go to any other box in the simulation. If > 0, this is the maximum number of cells
    *              a particle can have moved since the last Redistribute() call. Knowing this number
    *              allows an optimized MPI communication pattern to be used.
    */
    void Redistribute (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0,
                       bool remove_negative=true);


    /**
     * \brief Reorder particles on the tile given by lev and mfi using a the permutations array.
     *
     * permutations is a pointer to an array on the GPU of size numParticles() with
     * permutations[new index] = old index.
     *
     * \param lev
     * \param mfi
     * \param permutations
     *
     */
    template <class index_type>
    void ReorderParticles (int lev, const MFIter& mfi, const index_type* permutations);

    /**
     * \brief Sort particles on each tile such that particles adjacent in memory
     * are likely to map to adjacent cells. This ordering can be beneficial for performance
     * on GPU when deposition quantities onto a grid.
     *
     * idx_type = {0, 0, 0}: Sort particles to a cell centered grid
     * idx_type = {1, 1, 1}: Sort particles to a node centered grid
     * idx_type = {2, 2, 2}: Compromise between a cell and node centered grid.
     * This last option uses more memory than the fist two.
     * Mixed versions are also possible.
     *
     * \param idx_type
     *
     */
    void SortParticlesForDeposition (IntVect idx_type);

    /**
     * \brief Sort the particles on each tile by cell, using Fortran ordering.
     */
    void SortParticlesByCell ();

    /**
     * \brief Sort the particles on each tile by groups of cells, given an IntVect bin_size
     *
     * If bin_size is the zero vector, this operation is a no-op.
     *
     */
    void SortParticlesByBin (IntVect bin_size);

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

    std::array<Long, 3> ByteSpread () const;

    std::array<Long, 3> PrintCapacity () const;

    void ShrinkToFit ();

    /**
    * \brief Returns # of particles at specified the level.
    *
    * If "only_valid" is true it only counts valid particles.
    *
    * \param level
    * \param only_valid
    * \param only_local
    */

    Long NumberOfParticlesAtLevel (int level, bool only_valid = true, bool only_local = false) const;

    Vector<Long> NumberOfParticlesInGrid  (int level, bool only_valid = true, bool only_local = false) const;

    /**
    * \brief Returns # of particles at all levels
    *
    * If "only_valid" is true it only counts valid particles.
    *
    * \param only_valid
    * \param only_local
    */
    Long TotalNumberOfParticles (bool only_valid=true, bool only_local=false) const;


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
    * \brief Creates virtual particles for a given level that represent
    * in some capacity all particles at finer levels
    *
    * \param level
    * \param virts
    */
    void CreateVirtualParticles (int level, ParticleTileType& virts) const;

    /**
    * \brief Create ghost particles for a given level that are copies of particles
    * near coarse-\>fine boundaries in level-1
    *
    * \param level
    * \param ngrow
    * \param ghosts
    */
    void CreateGhostParticles (int level, int ngrow, ParticleTileType& ghosts) const;

    /**
    * \brief Add particles from a pbox to the grid at this level
    *
    * \param particles
    * \param level
    * \param nGrow
    */
    void AddParticlesAtLevel (ParticleTileType& particles, int level, int nGrow=0);


    /**
    * \brief Clear all the particles in this container. This does not free memory.
    */
    void clearParticles ();


    /**
    * \brief Copy particles from other to this ParticleContainer. Will clear all the
    * particles from this container first. local controls whether or not to call
    * Redistribute() after adding the particles.
    *
    * \param other the other pc to copy from
    * \param local whether to call redistribute after
    */
    template <class PCType,
              std::enable_if_t<IsParticleContainer<PCType>::value, int> foo = 0>
    void copyParticles (const PCType& other, bool local=false);

    /**
    * \brief Add particles from other to this ParticleContainer. local controls
    * whether or not to call Redistribute after adding the particles.
    *
    * \param other the other pc to copy from
    * \param local whether to call redistribute after
    */
    template <class PCType,
              std::enable_if_t<IsParticleContainer<PCType>::value, int> foo = 0>
    void addParticles (const PCType& other, bool local=false);

    /**
    * \brief Copy particles from other to this ParticleContainer. Will clear all the
    * particles from this container first. local controls whether or not to call
    * Redistribute() after adding the particles.
    *
    * This version conditionally copies based on a predicate function applied to
    * each particle.
    *
    * \tparam callable that takes a SuperParticle and returns a bool
    *
    * \param other the other pc to copy from
    * \param f function to apply to each particle as a filter
    * \param local whether to call redistribute after
    */
    template <class F, class PCType,
              std::enable_if_t<IsParticleContainer<PCType>::value, int> foo = 0,
              std::enable_if_t<! std::is_integral_v<F>, int> bar = 0>
    void copyParticles (const PCType& other, F&&f, bool local=false);

    /**
    * \brief Add particles from other to this ParticleContainer. local controls
    * whether or not to call Redistribute after adding the particles.
    *
    * This version conditionally copies based on a predicate function applied to
    * each particle.
    *
    * \tparam callable that takes a SuperParticle and returns a bool
    *
    * \param other the other pc to copy from
    * \param f function to apply to each particle as a filter
    * \param local whether to call redistribute after
    */
    template <class F, class PCType,
              std::enable_if_t<IsParticleContainer<PCType>::value, int> foo = 0,
              std::enable_if_t<! std::is_integral_v<F>, int> bar = 0>
    void addParticles (const PCType& other, F const& f, bool local=false);

    /**
    * \brief Write a contiguous chunk of real particle data to an ostream.
    *
    * \param data A pointer to the start of the buffer to write
    * \param size The number of elements to write
    * \param os The ostream into which to write the data
    */
    void WriteParticleRealData (void* data, size_t size, std::ostream& os) const;

    /**
    * \brief Read a contiguous chunk of real particle data from an istream.
    *
    * \param data A pointer to the start of the buffer into which to read
    * \param size The number of elements to read
    * \param os The istream from which to read the data
    * \param rd A RealDescriptor describing the type of the floating point data
    */
    void ReadParticleRealData (void* data, size_t size, std::istream& is);

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names vector of real component names, optional
     * \param int_comp_names vector of int component names, optional
     */
    void Checkpoint (const std::string& dir, const std::string& name,
                     const Vector<std::string>& real_comp_names = Vector<std::string>(),
                     const Vector<std::string>& int_comp_names = Vector<std::string>()) const
    {
        Checkpoint(dir, name, true, real_comp_names, int_comp_names);
    }

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *        This version allows the particle component names to be passed in.
     *        This overload exists for backwards compatibility. The is_checkpoint parameter is ignored.
     */
    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint,
                     const Vector<std::string>& real_comp_names = Vector<std::string>(),
                     const Vector<std::string>& int_comp_names = Vector<std::string>()) const;

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *        This version allows some components to be toggled off, if they
     *        don't need to be stored in the chk file.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
      * \param write_real_comp for each real component, whether or not we include that component in the file
      * \param write_int_comp for each integer component, whether or not we include that component in the file
     * \param real_comp_names vector of real component names
     * \param int_comp_names vector of int component names
     */
    void Checkpoint (const std::string& dir, const std::string& name,
                     const Vector<int>& write_real_comp,
                     const Vector<int>& write_int_comp,
                     const Vector<std::string>& real_comp_names,
                     const Vector<std::string>& int_comp_names) const;

     /**
      * \brief Writes particle data to disk in the AMReX native format.
      *
      * \tparam F function type
      *
      * \param dir The base directory into which to write (i.e. "plt00000")
      * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
      * \param write_real_comp for each real component, whether or not we include that component in the file
      * \param write_int_comp for each integer component, whether or not we include that component in the file
      * \param real_comp_names for each real component, a name to label the data with
      * \param int_comp_names for each integer component, a name to label the data with
      * \param f callable that returns whether a given particle should be written or not
      * \param is_checkpoint whether the data is written to a checkpoint or plotfile
      */
    template <class F>
    void WriteBinaryParticleData (const std::string& dir,
                                  const std::string& name,
                                  const Vector<int>& write_real_comp,
                                  const Vector<int>& write_int_comp,
                                  const Vector<std::string>& real_comp_names,
                                  const Vector<std::string>&  int_comp_names,
                                  F&& f, bool is_checkpoint=false) const;

    void CheckpointPre ();

    void CheckpointPost ();

    /**
     * \brief Restart from checkpoint
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void Restart (const std::string& dir, const std::string& file);

    /**
     * \brief Older version, for backwards compatibility
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param is_checkpoint Whether the particle id and cpu are included in the file.
     */
    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint);

    /**
     * \brief This version of WritePlotFile writes all components and assigns component names
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void WritePlotFile (const std::string& dir, const std::string& name) const;

    /**
     * \brief This version of WritePlotFile writes all components and assigns component names
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, std::enable_if_t<!std::is_same_v<std::decay_t<F>, Vector<std::string>>>* = nullptr>
    void WritePlotFile (const std::string& dir, const std::string& name, F&& f) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components.
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
     * \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFile (const std::string& dir, const std::string& name,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>&  int_comp_names, F&& f) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify
     * the names of the components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     */
    void WritePlotFile (const std::string& dir, const std::string& name,
                        const Vector<std::string>& real_comp_names) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify
     * the names of the components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, std::enable_if_t<!std::is_same_v<std::decay_t<F>, Vector<std::string>>>* = nullptr>
    void WritePlotFile (const std::string& dir, const std::string& name,
                        const Vector<std::string>& real_comp_names, F&& f) const;

    /**
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specific components.
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
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specific components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFile (const std::string& dir,
                        const std::string& name,
                        const Vector<int>& write_real_comp,
                        const Vector<int>& write_int_comp, F&& f) const;

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

    /**
     * \brief This is the most general version of WritePlotFile, which takes component
     * names and flags for whether to write each variable as components. Note that
     * the user should pass in vectors containing names of all the components, whether
     * they are written or not.
     *
     *  This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFile (const std::string& dir,
                        const std::string& name,
                        const Vector<int>& write_real_comp,
                        const Vector<int>& write_int_comp,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>&  int_comp_names,
                        F&& f) const;

    void WritePlotFilePre ();

    void WritePlotFilePost ();

    void WriteAsciiFile (const std::string& file);

    /**
     * \brief Return the underlying Vector (over AMR levels) of ParticleLevels.
     *        Const version.
     */
    const Vector<ParticleLevel>& GetParticles () const { return m_particles; }

    /**
     * \brief Return the underlying Vector (over AMR levels) of ParticleLevels.
     *        Non-const version.
     */
    Vector      <ParticleLevel>& GetParticles ()       { return m_particles; }

    /**
     * \brief Return the ParticleLevel for level "lev". Const version.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev the level at which to get the particles
     */
    const ParticleLevel& GetParticles (int lev) const { return m_particles[lev]; }

    /**
     * \brief Return the ParticleLevel for level "lev". Non-const version.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev the level at which to get the particles
     */
    ParticleLevel      & GetParticles (int lev)       { return m_particles[lev]; }

    /**
     * \brief Return the ParticleTile for level "lev", grid "grid" and tile "tile."
     *        Const version.
     *
     *        Here, grid and tile are integers that give the index and LocalTileIndex
     *        of the tile you want.
     *
     *        This is a runtime error if a ParticleTile at "grid" and "tile" has not been
     *        created yet.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param grid the index of the grid at which to get the particles
     *        \param tile the LocalTileIndex of the tile at which to get the particles
     */
    const ParticleTileType& ParticlesAt (int lev, int grid, int tile) const
    { return m_particles[lev].at(std::make_pair(grid, tile)); }

    /**
     * \brief Return the ParticleTile for level "lev", grid "grid" and tile "tile."
     *        Non-const version.
     *
     *        Here, grid and tile are integers that give the index and LocalTileIndex
     *        of the tile you want.
     *
     *        This is a runtime error if a ParticleTile at "grid" and "tile" has not been
     *        created yet.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param grid the index of the grid at which to get the particles
     *        \param tile the LocalTileIndex of the tile at which to get the particles
     */
    ParticleTileType&       ParticlesAt (int lev, int grid, int tile)
    { return m_particles[lev].at(std::make_pair(grid, tile)); }

    /**
     * \brief Return the ParticleTile for level "lev" and Iterator "iter".
     *        Const version.
     *
     *        Here, iter is either an MFIter or ParIter object pointing to the
     *        tile you want.
     *
     *        This is a runtime error if a ParticleTile at "iter" has not been
     *        created yet.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param iter MFIter or ParIter pointing to the tile to return
     */
    template <class Iterator>
    const ParticleTileType& ParticlesAt (int lev, const Iterator& iter) const
        { return ParticlesAt(lev, iter.index(), iter.LocalTileIndex()); }

    /**
     * \brief Return the ParticleTile for level "lev" and Iterator "iter".
     *        Non-const version.
     *
     *        Here, iter is either an MFIter or ParIter object pointing to the
     *        tile you want.
     *
     *        This is a runtime error if a ParticleTile at "iter" has not been
     *        created yet.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param iter MFIter or ParIter pointing to the tile to return
     */
    template <class Iterator>
    ParticleTileType&       ParticlesAt (int lev, const Iterator& iter)
        { return ParticlesAt(lev, iter.index(), iter.LocalTileIndex()); }

    /**
     * \brief Define and return the ParticleTile for level "lev", grid "grid" and tile "tile."
     *
     *        Here, grid and tile are integers that give the index and LocalTileIndex
     *        of the tile you want.
     *
     *        If a ParticleTile at "grid" and "tile" has not been created yet,
     *        this function call will add it. This call will also allocate space
     *        for any runtime-added particle components.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param grid the index of the grid at which to get the particles
     *        \param tile the LocalTileIndex of the tile at which to get the particles
     */
    ParticleTileType& DefineAndReturnParticleTile (int lev, int grid, int tile)
    {
        m_particles[lev][std::make_pair(grid, tile)].define(NumRuntimeRealComps(), NumRuntimeIntComps(), &m_soa_rdata_names, &m_soa_idata_names);

        return ParticlesAt(lev, grid, tile);
    }

    /**
     * \brief Define and return the ParticleTile for level "lev", and Iterator "iter".
     *
     *        Here, iter is either an MFIter or ParIter object pointing to the
     *        tile you want.
     *
     *        If a ParticleTile at "grid" and "tile" has not been created yet,
     *        this function call will add it. This call will also allocate space
     *        for any runtime-added particle components.
     *
     *        The ParticleLevel must already exist, meaning that the "resizeData()"
     *        method of this ParticleContainer has been called.
     *
     *        Note that, when using a ParticleContainer that has been constructed
     *        with an AmrCore*, "resizeData()" must be called *after* the grids
     *        have been created, meaning after the call to AmrCore::InitFromScratch
     *        or AmrCore::InitFromCheckpoint has been made.
     *
     *        \param lev  the level at which to get the particles
     *        \param iter MFIter or ParIter pointing to the tile to return
     */
    template <class Iterator>
    ParticleTileType& DefineAndReturnParticleTile (int lev, const Iterator& iter)
    {
        auto index = std::make_pair(iter.index(), iter.LocalTileIndex());
        m_particles[lev][index].define(NumRuntimeRealComps(), NumRuntimeIntComps(), &m_soa_rdata_names, &m_soa_idata_names);
        return ParticlesAt(lev, iter);
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

    void AssignCellDensitySingleLevel (int rho_index, MultiFab& mf, int level,
                                       int ncomp=1, int particle_lvl_offset = 0) const;

    template <typename P, typename Assignor=CellAssignor>
    IntVect Index (const P& p, int lev) const;

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
    template <typename P>
    bool PeriodicShift (P& p) const;

    void SetLevelDirectoriesCreated (bool tf) { levelDirectoriesCreated = tf; }

    bool GetLevelDirectoriesCreated () const { return levelDirectoriesCreated; }

    void SetUsePrePost (bool tf) const {
      usePrePost = tf;
    }
    bool GetUsePrePost () const {
      return usePrePost;
    }

    int GetMaxNextIDPrePost () const { return maxnextidPrePost; }
    Long GetNParticlesPrePost () const { return nparticlesPrePost; }

    void SetUseUnlink (bool tf) const {
      doUnlink = tf;
    }

    bool GetUseUnlink () const {
      return doUnlink;
    }

    void RedistributeCPU (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0,
                          bool remove_negative=true);

    void RedistributeGPU (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0,
                          bool remove_negative=true);

    Long superParticleSize() const { return superparticle_size; }

    void AddRealComp (std::string const & name, int communicate=1)
    {
        // names must be unique
        auto const it = std::find(m_soa_rdata_names.begin(), m_soa_rdata_names.end(), name);
        if (it != m_soa_rdata_names.end()) {
            throw std::runtime_error("AddRealComp: name '" + name + "' is already present in the SoA.");
        }
        m_soa_rdata_names.push_back(name);

        m_runtime_comps_defined = true;
        m_num_runtime_real++;
        h_redistribute_real_comp.push_back(communicate);
        SetParticleSize();
        this->resizeData();

        // resize runtime SoA
        for (int lev = 0; lev < numLevels(); ++lev) {
            for (ParIterType pti(*this,lev); pti.isValid(); ++pti) {
                auto& tile = DefineAndReturnParticleTile(lev, pti);
                auto np = tile.numParticles();
                if (np > 0) {
                    auto& soa = tile.GetStructOfArrays();
                    soa.resize(np);
                }
            }
        }
    }

    void AddRealComp (int communicate=1)
    {
        AddRealComp(getDefaultCompNameReal<ParticleType>(NArrayReal+m_num_runtime_real), communicate);
    }

    void AddIntComp (std::string const & name, int communicate=1)
    {
        // names must be unique
        auto const it = std::find(m_soa_idata_names.begin(), m_soa_idata_names.end(), name);
        if (it != m_soa_idata_names.end()) {
            throw std::runtime_error("AddIntComp: name '" + name + "' is already present in the SoA.");
        }
        m_soa_idata_names.push_back(name);

        m_runtime_comps_defined = true;
        m_num_runtime_int++;
        h_redistribute_int_comp.push_back(communicate);
        SetParticleSize();
        this->resizeData();

        // resize runtime SoA
        for (int lev = 0; lev < numLevels(); ++lev) {
            for (ParIterType pti(*this,lev); pti.isValid(); ++pti) {
                auto& tile = DefineAndReturnParticleTile(lev, pti);
                auto np = tile.numParticles();
                if (np > 0) {
                    auto& soa = tile.GetStructOfArrays();
                    soa.resize(np);
                }
            }
        }
    }

    void AddIntComp (int communicate=1)
    {
        AddIntComp(getDefaultCompNameInt<ParticleType>(NArrayInt+m_num_runtime_int), communicate);
    }

    int NumRuntimeRealComps () const { return m_num_runtime_real; }
    int NumRuntimeIntComps  () const { return m_num_runtime_int;  }

    int NumRealComps () const { return NArrayReal + NumRuntimeRealComps(); }
    int NumIntComps  () const { return NArrayInt  + NumRuntimeIntComps() ; }

    /** Resize the Real runtime components (SoA)
     *
     * @param new_size new number of Real runtime components
     * @param communicate participate this component in redistribute
     */
    void ResizeRuntimeRealComp (int new_size, bool communicate);

    /** Resize the Int runtime components (SoA)
     *
     * @param new_size new number of integer runtime components
     * @param communicate participate this component in redistribute
     */
    void ResizeRuntimeIntComp (int new_size, bool communicate);

    /** type trait to translate one particle container to another, with changed allocator */
    template <template<class> class NewAllocator=amrex::DefaultAllocator>
    using ContainerLike = amrex::ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, NewAllocator>;

    /** Create an empty particle container
     *
     * This creates a new AMReX particle container type with same compile-time
     * and run-time attributes. But, it can change its allocator. This is
     * helpful when creating temporary particle buffers for filter operations
     * and device-to-host copies.
     *
     * @tparam Allocator AMReX allocator, e.g., amrex::PinnedArenaAllocator
     * @return an empty particle container
     * */
    template <template<class> class NewAllocator=amrex::DefaultAllocator>
    ContainerLike<NewAllocator>
    make_alike () const
    {
        ContainerLike<NewAllocator> tmp(m_gdb);

        std::vector<std::string> const real_names = this->GetRealSoANames();
        std::vector<std::string> real_ct_names(NArrayReal);
        for (int ic = 0; ic < NArrayReal; ++ic) { real_ct_names.at(ic) = real_names[ic]; }

        std::vector<std::string> const int_names = this->GetIntSoANames();
        std::vector<std::string> int_ct_names(NArrayInt);
        for (int ic = 0; ic < NArrayInt; ++ic) { int_ct_names.at(ic) = int_names[ic]; }

        tmp.SetSoACompileTimeNames(real_ct_names, int_ct_names);

        // add runtime real comps to tmp
        for (int ic = 0; ic < this->NumRuntimeRealComps(); ++ic) {
            tmp.AddRealComp(real_names.at(ic + NArrayReal));
        }

        // add runtime int comps to tmp
        for (int ic = 0; ic < this->NumRuntimeIntComps(); ++ic) {
            tmp.AddIntComp(int_names.at(ic + NArrayInt));
        }

        tmp.h_redistribute_real_comp = h_redistribute_real_comp;
        tmp.h_redistribute_int_comp = h_redistribute_int_comp;

        return tmp;
    }

    Vector<int> h_redistribute_real_comp;
    Vector<int> h_redistribute_int_comp;

    //! Variables for i/o optimization saved for pre and post checkpoint
    mutable bool levelDirectoriesCreated;
    mutable bool usePrePost;
    mutable bool doUnlink;
    int maxnextidPrePost;
    mutable int nOutFilesPrePost;
    Long nparticlesPrePost;
    Vector<Long> nParticlesAtLevelPrePost;
    mutable Vector<Vector<int>>  whichPrePost;
    mutable Vector<Vector<int>>  countPrePost;
    mutable Vector<Vector<Long>> wherePrePost;
    mutable std::string HdrFileNamePrePost;
    mutable Vector<std::string> filePrefixPrePost;

protected:

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
    template <typename P>
    bool Where (const P& prt, ParticleLocData& pld,
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
    template <typename P>
    bool EnforcePeriodicWhere (P& prt, ParticleLocData& pld,
                               int lev_min = 0, int lev_max = -1, int local_grid=-1) const;

public:
    void
    WriteParticles (int level, std::ofstream& ofs, int fnum,
                    Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                    const Vector<int>& write_real_comp, const Vector<int>& write_int_comp,
                    const Vector<std::map<std::pair<int, int>,IntVector>>& particle_io_flags, bool is_checkpoint) const;
#ifdef AMREX_USE_HDF5
#include "AMReX_ParticlesHDF5.H"
#endif

    /** Overwrite the default names for the compile-time SoA components */
    void SetSoACompileTimeNames (std::vector<std::string> const & rdata_name, std::vector<std::string> const & idata_name);

    /** Get the names for the real SoA components **/
    std::vector<std::string> GetRealSoANames () const {return m_soa_rdata_names;}

    /** Get the names for the int SoA components **/
    std::vector<std::string> GetIntSoANames () const {return m_soa_idata_names;}

    /** Check if a container has a ParticleReal component
     *
     * @param name component name to check
     * @return true if found, else false
     */
    bool HasRealComp (std::string const & name);

    /** Check if a container has an Integer component
     *
     * @param name component name to check
     * @return true if found, else false
     */
    bool HasIntComp (std::string const & name);

    /** Get the ParticleReal SoA index of a component
     *
     * This throws a runtime exception if the component does not exist.
     *
     * @param name component name to query index for
     * @return zero-based index
     */
    int GetRealCompIndex (std::string const & name);

    /** Get the Integer SoA index of a component
     *
     * This throws a runtime exception if the component does not exist.
     *
     * @param name component name to query index for
     * @return zero-based index
     */
    int GetIntCompIndex (std::string const & name);

protected:

    template <class RTYPE>
    void ReadParticles (int cnt, int grd, int lev, std::ifstream& ifs, int finest_level_in_file, bool convert_ids);

    void SetParticleSize ();

    DenseBins<typename ParticleTileType::ParticleTileDataType> m_bins;

private:
    virtual void particlePostLocate (ParticleType& /*p*/, const ParticleLocData& /*pld*/,
                                     const int /*lev*/) {}

    virtual void correctCellVectors (int /*old_index*/, int /*new_index*/,
                                     int /*grid*/, const ParticleType& /*p*/) {}

    void RedistributeMPI (std::map<int, Vector<char> >& not_ours,
                          int lev_min = 0, int lev_max = 0, int nGrow = 0, int local=0);

    template <typename P>
    void locateParticle (P& p, ParticleLocData& pld,
                         int lev_min, int lev_max, int nGrow, int local_grid=-1) const;

    void Initialize ();

    bool m_runtime_comps_defined{false};
    int m_num_runtime_real{0};
    int m_num_runtime_int{0};

    size_t particle_size, superparticle_size;
    int num_real_comm_comps, num_int_comm_comps;
    Vector<ParticleLevel> m_particles;

    // names of both compile-time and runtime Real and Int SoA data
    std::vector<std::string> m_soa_rdata_names;
    std::vector<std::string> m_soa_idata_names;
};

template <int T_NStructReal, int T_NStructInt, int T_NArrayReal, int T_NArrayInt, template<class> class Allocator, class CellAssignor>
using ParticleContainer = ParticleContainer_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <int T_NArrayReal, int T_NArrayInt, template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParticleContainerPureSoA = ParticleContainer_impl<SoAParticle<T_NArrayReal, T_NArrayInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

}

#include "AMReX_ParticleInit.H"
#include "AMReX_ParticleContainerI.H"
#include "AMReX_ParticleIO.H"

#ifdef AMREX_USE_HDF5
#include "AMReX_ParticleHDF5.H"
#endif

#endif /*_PARTICLES_H_*/
