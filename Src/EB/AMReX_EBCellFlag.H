#ifndef AMREX_EBCELLFLAG_H_
#define AMREX_EBCELLFLAG_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabFactory.H>

#include <cstdint>
#include <map>
#include <iosfwd>

/**
 * \file AMReX_EBCellFlag.H
 *
 * Bit-packed flags describing whether cells are regular, covered, or cut by EB.
 */

namespace amrex {

//! Bit-coded helper that classifies EB cells and neighbor connectivity.
class EBCellFlag
{
public:

    //! Default-construct a flag with all bits cleared (treated as regular).
    constexpr EBCellFlag () noexcept = default;
    ~EBCellFlag () noexcept = default;
    constexpr EBCellFlag (const EBCellFlag& rhs) noexcept = default;
    constexpr EBCellFlag (EBCellFlag&& rhs) noexcept = default;
    constexpr EBCellFlag& operator= (const EBCellFlag& rhs) noexcept = default;
    constexpr EBCellFlag& operator= (EBCellFlag&& rhs) noexcept = default;

    //! Construct from a raw bit pattern.
    explicit constexpr EBCellFlag (uint32_t i) noexcept : flag(i) {}

    //! Assign a raw bit pattern to this flag.
    constexpr EBCellFlag& operator= (uint32_t i) noexcept { flag = i; return *this; }

    //! Disabled on purpose—EBCellFlag values are not additive.
    AMREX_GPU_HOST_DEVICE
    EBCellFlag& operator+= (const EBCellFlag& /* rhs */) {
        amrex::Abort("EBCellFlag::operator+= not supported");
        return *this;
    }

    //! Mark this cell as purely regular (single volume fraction, no cuts).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setRegular () noexcept {
        flag &= zero_lower_mask; //!< clean lower bits
        flag |= regular_bits;
        flag |= single_vof_bits;
    }

    //! Mark the cell as completely covered (solid).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setCovered () noexcept {
        flag = covered_value;
    }

    //! Mark the cell as single-valued (one cut)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setSingleValued () noexcept {
        flag &= zero_lower_mask;
        flag |= single_valued_bits;
        flag |= single_vof_bits;
    }

    //! Mark the cell as multi-valued and store the number of VoFs (\p n).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setMultiValued (int n) noexcept {
        flag &= zero_lower_mask;
        flag |= multi_valued_bits;
        AMREX_ASSERT(n >= 2 && n <= 7);
        flag |= static_cast<uint32_t>(n) << pos_numvofs;
    }

    //! Compare two flags for bitwise equality.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator==(const EBCellFlag& a_input) const noexcept
    {
        return (flag == (a_input.flag));
    }

    //! Compare two flags for inequality.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!=(const EBCellFlag& a_input) const noexcept
    {
        return (flag != (a_input.flag));
    }

    //! Number of volume-of-fluid entries for this cell (1 for regular or single-valued, 0 for covered).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int getNumVoFs () const noexcept
    {
        return static_cast<int>((flag & one_numvofs_mask) >> pos_numvofs);
    }

    //! True if the stored type bits correspond to a regular cell.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isRegular () const noexcept {
        return (flag & one_type_mask) == regular_bits;
    }

    //! True if the cell stores a single cut (one VoF).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isSingleValued () const noexcept {
        return (flag & one_type_mask) == single_valued_bits;
    }

    //! True if the cell has multiple VoFs.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isMultiValued () const noexcept {
        return (flag & one_type_mask) == multi_valued_bits;
    }

    //! True if the cell is fully covered (solid).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isCovered () const noexcept {
        return (flag & one_type_mask) == covered_bits;
    }

    //! Return whether the EB cut is connected to the neighbor offset by \p iv.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isConnected (const IntVect& iv) const noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return flag & (1 << n);
    }

    //! Return whether the EB cut connects to neighbor at offset (i,j,k).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isConnected (int i, int j, int k) const  noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return flag & (1 << n);
    }

    //! Inverse of isConnected(IntVect).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isDisconnected (const IntVect& iv) const noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return !(flag & (1 << n));
    }

    //! Inverse of isConnected(i,j,k).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isDisconnected (int i, int j, int k) const  noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return !(flag & (1 << n));
    }

    //! Clear all neighbor-connection bits.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected () noexcept {
        flag &= one_lower_mask;
    }

    //! Clear connectivity to neighbor at offset \p iv.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected (const IntVect& iv) noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag &= ~(1 << n);
    }

    //! Clear connectivity to neighbor at offset (i,j,k).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected (int i, int j, int k) noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag &= ~(1 << n);
    }

    //! Mark all neighbors as connected.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected () noexcept {
        flag |= zero_lower_mask;
    }

    //! Mark connectivity to neighbor offset by \p iv.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected (const IntVect& iv) noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag |= 1 << n;
    }

    //! Mark connectivity to neighbor at offset (i,j,k).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected (int i, int j, int k) noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag |= 1 << n;
    }

    //! Count how many neighbors are connected (including itself).
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int numNeighbors () const noexcept {
        int n = 1;
        n += isConnected(-1,-1, 0);
        n += isConnected( 0,-1, 0);
        n += isConnected( 1,-1, 0);
        //
        n += isConnected(-1, 0, 0);
        n += isConnected( 1, 0, 0);
        //
        n += isConnected(-1, 1, 0);
        n += isConnected( 0, 1, 0);
        n += isConnected( 1, 1, 0);
#if (AMREX_SPACEDIM == 3)
        n += isConnected(-1,-1,-1);
        n += isConnected( 0,-1,-1);
        n += isConnected( 1,-1,-1);
        //
        n += isConnected(-1, 0,-1);
        n += isConnected( 0, 0,-1);
        n += isConnected( 1, 0,-1);
        //
        n += isConnected(-1, 1,-1);
        n += isConnected( 0, 1,-1);
        n += isConnected( 1, 1,-1);
        //
        n += isConnected(-1,-1, 1);
        n += isConnected( 0,-1, 1);
        n += isConnected( 1,-1, 1);
        //
        n += isConnected(-1, 0, 1);
        n += isConnected( 0, 0, 1);
        n += isConnected( 1, 0, 1);
        //
        n += isConnected(-1, 1, 1);
        n += isConnected( 0, 1, 1);
        n += isConnected( 1, 1, 1);
#endif
        return n;
    }

    //! Raw 32-bit encoding stored in this flag.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint32_t getValue() const noexcept
    {
        return flag;
    }

    //! Convenience flag representing a completely regular cell.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr EBCellFlag TheDefaultCell () { return EBCellFlag{default_value}; }

    // DO NOT test covered cells using cell == TheCoveredCell()
    // Use cell.isCovered() instead.
    // Each covered cell is covered in its own way.
    //! Convenience flag representing a fully covered cell.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr EBCellFlag TheCoveredCell () { return EBCellFlag{covered_value}; }

private:

    //! masks for lowest 5 bits (i.e., bit number 0-4)
    static constexpr uint32_t one_lower_mask      =  0x1f;
    static constexpr uint32_t zero_lower_mask     = ~0x1f;

    //! masks lowest 2 bits (i.e., bit number 0-1)
    static constexpr uint32_t one_type_mask       =  0x3;
    static constexpr uint32_t zero_type_mask      = ~0x3;

    //! masks for bit number 2-4
    static constexpr uint32_t one_numvofs_mask    =  0x1c;
    static constexpr uint32_t zero_numvofs_mask   = ~0x1c;

    //! these represent cell types
    //! note that single-valued means single-value cut cell
    static constexpr uint32_t regular_bits        =  0x0;
    static constexpr uint32_t single_valued_bits  =  0x1;
    static constexpr uint32_t multi_valued_bits   =  0x2;
    static constexpr uint32_t covered_bits        =  0x3;

    //! this represent single vof (regular is considered as single vof too)
    static constexpr uint32_t single_vof_bits     =  0x4;

    //! There are 32 bits.  The lowest 2 bits are used for cell type:
    //! regular, single-valued, multi-valued, and covered.  The next 3
    //! bits are for the number of vofs.  The rest 27 bits are used for
    //! connectivity with neighbors.

    static constexpr int w_lower_mask = 5;
    static constexpr int w_type       = 2;
    static constexpr int w_numvofs    = 3;
    static constexpr int pos_numvofs  = 2;
#if 0
    static constexpr std::array<std::array<std::array<int,3>,3>,3> pos_ngbr
        {{ std::array<std::array<int,3>,3>{{{ w_lower_mask   , w_lower_mask+ 1, w_lower_mask+ 2 },
                                            { w_lower_mask+ 3, w_lower_mask+ 4, w_lower_mask+ 5 },
                                            { w_lower_mask+ 6, w_lower_mask+ 7, w_lower_mask+ 8 }}},
           std::array<std::array<int,3>,3>{{{ w_lower_mask+ 9, w_lower_mask+10, w_lower_mask+11 },
                                            { w_lower_mask+12, w_lower_mask+13, w_lower_mask+14 },
                                            { w_lower_mask+15, w_lower_mask+16, w_lower_mask+17 }}},
           std::array<std::array<int,3>,3>{{{ w_lower_mask+18, w_lower_mask+19, w_lower_mask+20 },
                                            { w_lower_mask+21, w_lower_mask+22, w_lower_mask+23 },
                                            { w_lower_mask+24, w_lower_mask+25, w_lower_mask+26 }}} }};
#endif

    //! regular connected with all neighbors
#if AMREX_SPACEDIM == 3
    static constexpr uint32_t default_value = zero_lower_mask | regular_bits | single_vof_bits;
#else
    static constexpr uint32_t default_value = 0x7fc004;  //!< zero out "3d" neighbors
#endif
//    static constexpr uint32_t covered_value = zero_lower_mask | covered_bits;
    static constexpr uint32_t covered_value = 0x40003;  //!< zero out all neighbors

    uint32_t flag = default_value;
};

template <>
struct IsStoreAtomic<EBCellFlag> : std::true_type {};

//! FAB storing EBCellFlag data per cell.
class EBCellFlagFab final
    : public BaseFab<EBCellFlag>
{
public:

    //! Construct an empty FAB that allocates from arena \p ar when defined later.
    explicit EBCellFlagFab (Arena* ar) noexcept;

    //! Construct and allocate immediately for box \p b.
    explicit EBCellFlagFab (const Box& b, int ncomp, Arena* ar);

    //! General constructor controlling allocation and sharing.
    explicit EBCellFlagFab (const Box& b,
                            int        ncomp=1,
                            bool       alloc=true,
                            bool       shared=false,
                            Arena*     ar = nullptr);

    //! Restricted constructor used when creating alias FABs.
    EBCellFlagFab (const EBCellFlagFab& rhs, MakeType make_type, int scomp, int ncomp);

    EBCellFlagFab () = default;
    EBCellFlagFab (EBCellFlagFab&& rhs) = default;
    ~EBCellFlagFab () override = default;

    EBCellFlagFab (const EBCellFlagFab&) = delete;
    EBCellFlagFab& operator= (const EBCellFlagFab&) = delete;
    EBCellFlagFab& operator= (EBCellFlagFab&&) = delete;

    //! Returns FabType (regular, covered, etc.) for the entire FAB.
    FabType getType () const noexcept { return m_type; }

    //! Returns FabType for a specific sub-box \p bx.
    FabType getType (const Box& bx) const noexcept;

    //! Recompute FabType using the given number of grow cells.
    void resetType (int ng);

    //! Count regular cells within \p bx.
    int getNumRegularCells (const Box& bx) const noexcept;

    //! Count cut cells within \p bx.
    int getNumCutCells (const Box& bx) const noexcept;

    //! Count covered cells within \p bx.
    int getNumCoveredCells (const Box& bx) const noexcept;

    //! Manually override the cached FabType.
    void setType (FabType t) noexcept { m_type = t; }

    struct NumCells {
        int nregular = 0;
        int nsingle = 0;
        int nmulti = 0;
        int ncovered = 0;
        FabType type = FabType::undefined;
    };

private:
    FabType m_type = FabType::undefined;
    mutable std::map<Box,NumCells> m_typemap;
};

//! Stream helper that prints the underlying flag bits.
std::ostream& operator<< (std::ostream& os, const EBCellFlag& flag);

}

#endif
