#ifndef BL_MFITER_H_
#define BL_MFITER_H_

#include <memory>

#include <AMReX_Arena.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_IntVect.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_RealBox.H>

#include <AMReX_Gpu.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

#ifdef AMREX_USE_GPU
    inline bool TilingIfNotGPU () { return Gpu::notInLaunchRegion(); }
#else
    inline constexpr bool TilingIfNotGPU () { return true; }
#endif

template<class T> class FabArray;

struct MFItInfo
{
    bool do_tiling;
    bool dynamic;
    IntVect tilesize;
    MFItInfo () 
        : do_tiling(false), dynamic(false), tilesize(IntVect::TheZeroVector()) {}
    MFItInfo& EnableTiling (const IntVect& ts = FabArrayBase::mfiter_tile_size) {
        do_tiling = true;
        tilesize = ts;
        return *this;
    }
    MFItInfo& SetDynamic (bool f) {
        dynamic = f;
        return *this;
    }
};

class MFIter
{
public:
    //! Flags that specify tiling modes. All these flags are off by default.
    enum Flags {
        //!Tiling: Enabling the tiling mode 
        Tiling        = 0x01,
        /** 
        * \brief AllBoxes: If on, all threads/workers loop over all boxes without tiling.
        * This essentially loops over indexMap.
        * Note that many functions won't work with this.
        */
        AllBoxes      = 0x02, 
        //! NoTeamBarrier: This option is for Team only. If on, there is no barrier in MFIter dtor.
        NoTeamBarrier = 0x04, 
        //! SkipInit: Used by MFGhostIter
	SkipInit      = 0x08  
    };

#ifdef AMREX_USE_GPU
    enum MFReducer { SUM = 0, MAX, MIN };
#endif

    /** 
    * \brief The default constructor does not enable tiling (flags_ is set to 0 by default). 
    * However, the tiling mode can be enabled by explicitly turning the tiling flag on (i.e. flags_ = Tiling). 
    * The tile size in this case is defined by FabArrayBase::mfiter_tile_size
    */
    explicit MFIter (const FabArrayBase& fabarray,
		     unsigned char       flags_=0);

    //! Enable tiling with the default tile size, which is defined by FabArrayBase::mfiter_tile_size
    MFIter (const FabArrayBase& fabarray, 
	    bool                do_tiling); 

    //! Enable tiling with explicit tile size and flags (See type Flags for more information)
    MFIter (const FabArrayBase& fabarray, 
	    const IntVect&      tilesize,
	    unsigned char       flags_=0);

    MFIter (const BoxArray& ba, const DistributionMapping& dm, unsigned char flags_=0);

    MFIter (const BoxArray& ba, const DistributionMapping& dm, bool do_tiling); 

    MFIter (const BoxArray& ba, const DistributionMapping& dm,
	    const IntVect& tilesize, unsigned char flags_=0);

    MFIter (const FabArrayBase& fabarray, const MFItInfo& info);

    MFIter (const BoxArray& ba, const DistributionMapping& dm, const MFItInfo& info);

    MFIter (MFIter&& rhs) = default;

    // dtor
    ~MFIter ();

    //! Return the tile Box at the current index.
    Box tilebox () const;

    //! Return the tilebox with provided nodal flag
    Box tilebox (const IntVect& nodal) const;

    //! Return the tilebox with provided nodal flag and grown cells
    Box tilebox (const IntVect& nodal, const IntVect& ngrow) const;

    //! Return the dir-nodal (or all nodal if dir<0) Box at the current index.
    Box nodaltilebox (int dir=-1) const;

    //! Return the tile box at the current index grown to include ghost cells.
    Box growntilebox (int ng=-1000000) const;

    Box growntilebox (const IntVect& ng) const;

    //! Return the dir-nodal (or all nodal if dir<0) box grown to include ghost cells.
    Box grownnodaltilebox (int dir=-1, int ng=-1000000) const;

    Box grownnodaltilebox (int dir, const IntVect& ng) const;

    //! Return the valid Box in which the current tile resides.
    Box validbox () const { return fabArray.box((*index_map)[currentIndex]); }

    //! Return the Box of the FAB at which we currently point.
    Box fabbox () const { return fabArray.fabbox((*index_map)[currentIndex]); }

    //! Increment iterator to the next tile we own.
#if defined(_OPENMP)
    void operator++ () {
        if (dynamic) {
#pragma omp atomic capture
            currentIndex = nextDynamicIndex++;
        } else {
            ++currentIndex;
        }
    }
#elif !defined(AMREX_USE_GPU)
    void operator++ () { ++currentIndex; }
#else
    void operator++ ();
#endif

    //! Is the iterator valid i.e. is it associated with a FAB?
    bool isValid () const { return currentIndex < endIndex; }

    //! The index into the underlying BoxArray of the current FAB.
    int index () const { return (*index_map)[currentIndex]; }

    //! The number of indices.
    int length () const { return (endIndex - beginIndex); }

    //! The current local tile index in the current grid;
    int LocalTileIndex () const {return local_tile_index_map ? (*local_tile_index_map)[currentIndex] : 0;}
  
    //! The the number of tiles in the current grid;
    int numLocalTiles() const {return num_local_tiles ? (*num_local_tiles)[currentIndex] : 1;}

#ifdef AMREX_USE_GPU
    //! Maintain a list of values to reduce.
    template<typename T>
    T* add_reduce_value(T* val, MFReducer r) { return val; }

    Real* add_reduce_value(Real* val, MFReducer r);

    //! Reduce over the values in the list.
    void reduce();
#endif

    /**
    * \brief Return local index into the vector of fab pointers, m_fabs_v
    * When AllBoxes is on, local_index_map is a nullptr and local index is current index.
    */
    int LocalIndex () const { return local_index_map ? (*local_index_map)[currentIndex] : currentIndex; }

    //! Constant reference to FabArray over which we're iterating.
    const FabArrayBase& theFabArrayBase () const { return fabArray; }

    int tileIndex () const {return currentIndex;}

    const DistributionMapping& DistributionMap () const { return fabArray.DistributionMap(); }

protected:

    std::unique_ptr<FabArray<FArrayBox> > m_fa;  // This must be the first memeber!

    const FabArrayBase& fabArray;

    IntVect tile_size;

    unsigned char flags;
    int           currentIndex;
    int           beginIndex;
    int           endIndex;
    IndexType     typ;

    bool          dynamic;

    const Vector<int>* index_map;
    const Vector<int>* local_index_map;
    const Vector<Box>* tile_array;
    const Vector<int>* local_tile_index_map;
    const Vector<int>* num_local_tiles;

#ifdef AMREX_USE_GPU
    mutable Real* real_reduce_val;

    mutable MFReducer reducer;

    mutable Vector<Real> real_reduce_list;

    mutable Vector<Real*> real_device_reduce_list;
#endif

    static int nextDynamicIndex;
  
    void Initialize ();
};

//! Iterate over ghost cells.  Lots of MFIter functions do not work.
class MFGhostIter
    :
    public MFIter
{
public:
    explicit MFGhostIter (const FabArrayBase& fabarray);
private:
    void Initialize ();
    FabArrayBase::TileArray lta;
};

inline Arena* The_MFIter_Arena () { return The_Device_Arena(); }

}

#endif
