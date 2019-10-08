#ifndef AMREX_EBCELLFLAG_H_
#define AMREX_EBCELLFLAG_H_

#include <cstdint>
#include <iostream>
#include <map>
#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabFactory.H>

namespace amrex {

class EBCellFlag
{
public:

    AMREX_GPU_HOST_DEVICE
    EBCellFlag () noexcept {}

    EBCellFlag (const EBCellFlag& rhs) noexcept = default;

    AMREX_GPU_HOST_DEVICE
    explicit constexpr EBCellFlag (uint32_t i) noexcept : flag(i) {}

    EBCellFlag& operator= (const EBCellFlag& rhs) noexcept = default;

    AMREX_GPU_HOST_DEVICE
    EBCellFlag& operator= (uint32_t i) noexcept { flag = i; return *this; }

    AMREX_GPU_HOST_DEVICE
    EBCellFlag& operator+= (const EBCellFlag& /* rhs */) {
        amrex::Abort("EBCellFlag::operator+= not supported");
        return *this;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setRegular () noexcept {
        flag &= zero_lower_mask; //!< clean lower bits
        flag |= regular_bits;
        flag |= single_vof_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setCovered () noexcept {
        flag &= zero_lower_mask;
        flag |= covered_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setSingleValued () noexcept {
        flag &= zero_lower_mask;
        flag |= single_vauled_bits;
        flag |= single_vof_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setMultiValued (int n) noexcept {
        flag &= zero_lower_mask;
        flag |= multi_valued_bits;
        BL_ASSERT(n >= 2 && n <= 7);
        flag |= static_cast<uint32_t>(n) << pos_numvofs;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator==(const EBCellFlag& a_input) const noexcept
    {
        return (flag == (a_input.flag));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!=(const EBCellFlag& a_input) const noexcept
    {
        return (flag != (a_input.flag));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int getNumVoFs () const noexcept
    {
        return (flag & one_numvofs_mask) >> pos_numvofs;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isRegular () const noexcept {
        return (flag & one_type_mask) == regular_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isSingleValued () const noexcept {
        return (flag & one_type_mask) == single_vauled_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isMultiValued () const noexcept {
        return (flag & one_type_mask) == multi_valued_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isCovered () const noexcept {
        return (flag & one_type_mask) == covered_bits;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isConnected (const IntVect& iv) const noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return flag & (1 << n);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isConnected (int i, int j, int k) const  noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return flag & (1 << n);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isDisconnected (const IntVect& iv) const noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return !(flag & (1 << n));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isDisconnected (int i, int j, int k) const  noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return !(flag & (1 << n));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected () noexcept {
        flag &= one_lower_mask;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected (const IntVect& iv) noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag &= ~(1 << n);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setDisconnected (int i, int j, int k) noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag &= ~(1 << n);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected () noexcept {
        flag |= zero_lower_mask;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected (const IntVect& iv) noexcept {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag |= 1 << n;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setConnected (int i, int j, int k) noexcept {
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag |= 1 << n;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint32_t getValue() const noexcept
    {
        return flag;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr EBCellFlag TheDefaultCell () { return EBCellFlag{default_value}; }

    // DO NOT test covered cells using cell == TheCoveredCell()
    // Use cell.isCovered() instead.
    // Each covered cell is covered in its own way.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
    static constexpr uint32_t single_vauled_bits  =  0x1;
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

class EBCellFlagFab
    : public BaseFab<EBCellFlag>
{
public:

    explicit EBCellFlagFab (Arena* ar) noexcept;

    explicit EBCellFlagFab (const Box& b, int ncomp, Arena* ar);

    explicit EBCellFlagFab (const Box& b,
                            int        ncomp=1,
                            bool       alloc=true,
                            bool       shared=false,
                            Arena*     ar = nullptr);

    EBCellFlagFab (const EBCellFlagFab& rhs, MakeType make_type, int scomp, int ncomp);

    EBCellFlagFab () = default;
    EBCellFlagFab (EBCellFlagFab&& rhs) = default;
    ~EBCellFlagFab () = default;

    EBCellFlagFab (const EBCellFlagFab&) = delete;
    EBCellFlagFab& operator= (const EBCellFlagFab&) = delete;
    EBCellFlagFab& operator= (EBCellFlagFab&&) = delete;

    FabType getType () const noexcept { return m_type; }

    FabType getType (const Box& bx) const noexcept;

    void setType (FabType t) noexcept { m_type = t; }

private:
    FabType m_type = FabType::undefined;
    mutable std::map<Box,FabType> m_typemap;
};

std::ostream& operator<< (std::ostream& os, const EBCellFlag& flag);

}

#endif
