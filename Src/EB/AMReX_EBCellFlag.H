#ifndef AMREX_EBCELLFLAG_H_
#define AMREX_EBCELLFLAG_H_

#include <array>
#include <cstdint>
#include <iostream>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabFactory.H>

namespace amrex {

class EBCellFlag
{
public:

    EBCellFlag () = default;
    EBCellFlag (const EBCellFlag& rhs) noexcept = default;
    constexpr EBCellFlag (uint32_t i) : flag(i) {}
    EBCellFlag& operator= (const EBCellFlag& rhs) noexcept = default;
    EBCellFlag& operator= (uint32_t i) { flag = i; return *this; }

    AMREX_GPU_HOST_DEVICE
    EBCellFlag& operator+= (const EBCellFlag& /* rhs */) {
        amrex::Abort("EBCellFlag::operator+= not supported");
        return *this;
    }

    AMREX_GPU_HOST_DEVICE
    void setRegular () {
        flag &= zero_lower_mask; //!< clean lower bits
        flag |= regular_bits;
        flag |= single_vof_bits;
    }

    AMREX_GPU_HOST_DEVICE
    void setCovered () {
        flag &= zero_lower_mask;
        flag |= covered_bits;
    }

    AMREX_GPU_HOST_DEVICE
    void setSingleValued () {
        flag &= zero_lower_mask;
        flag |= single_vauled_bits;
        flag |= single_vof_bits;
    }

    AMREX_GPU_HOST_DEVICE
    void setMultiValued (int n) {
        flag &= zero_lower_mask;
        flag |= multi_valued_bits;
        BL_ASSERT(n >= 2 && n <= 7);
        flag |= static_cast<uint32_t>(n) << pos_numvofs;
    }

    AMREX_GPU_HOST_DEVICE
    bool operator==(const EBCellFlag& a_input) const
    {
        return (flag == (a_input.flag));
    }

    AMREX_GPU_HOST_DEVICE
    bool operator!=(const EBCellFlag& a_input) const
    {
        return (flag != (a_input.flag));
    }

    AMREX_GPU_HOST_DEVICE
    int getNumVoFs () const
    {
        return (flag & one_numvofs_mask) >> pos_numvofs;
    }

    AMREX_GPU_HOST_DEVICE
    bool isRegular () const {
        return (flag & one_type_mask) == regular_bits;
    }

    AMREX_GPU_HOST_DEVICE
    bool isSingleValued () const {
        return (flag & one_type_mask) == single_vauled_bits;
    }

    AMREX_GPU_HOST_DEVICE
    bool isMultiValued () const {
        return (flag & one_type_mask) == multi_valued_bits;
    }

    AMREX_GPU_HOST_DEVICE
    bool isCovered () const {
        return (flag & one_type_mask) == covered_bits;
    }

    AMREX_GPU_HOST_DEVICE
    bool isConnected (const IntVect& iv) const {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        return flag & (1 << n);
    }

    AMREX_GPU_HOST_DEVICE
    void setDisconnected () {
        flag &= one_lower_mask;
    }

    AMREX_GPU_HOST_DEVICE
    void setDisconnected (const IntVect& iv) {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag &= ~(1 << n);
    }

    AMREX_GPU_HOST_DEVICE
    void setConnected () {
        flag |= zero_lower_mask;
    }

    AMREX_GPU_HOST_DEVICE
    void setConnected (const IntVect& iv) {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = w_lower_mask + 13 + i + 3*j + k*9;  // pos_ngbr[k+1,j+1,i+1]
        flag |= 1 << n;
    }

    AMREX_GPU_HOST_DEVICE
    uint32_t getValue() const
    {
        return flag;
    }

    AMREX_GPU_HOST_DEVICE
    static constexpr EBCellFlag TheDefaultCell () { return {default_value}; }

    AMREX_GPU_HOST_DEVICE
    static constexpr EBCellFlag TheCoveredCell () { return {covered_value}; }

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


class EBCellFlagFab
    : public BaseFab<EBCellFlag>
{
public:
    using BaseFab<EBCellFlag>::BaseFab;

    FabType getType () const { return m_type; }
    FabType getType (const Box& bx) const;

    const Box& getRegion () const { return m_region; }

    void setType (FabType t) { m_type = t; }
    void setRegion (const Box& bx) { m_region = bx; }
private:
    FabType m_type = FabType::undefined;
    Box m_region;
};

std::ostream& operator<< (std::ostream& os, const EBCellFlag& flag);

}

#endif
