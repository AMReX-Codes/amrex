#ifndef AMREX_EBFACEFLAG_H_
#define AMREX_EBFACEFLAG_H_

#include <array>
#include <cstdint>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabFactory.H>

namespace amrex {

class EBFaceFlag
{
public:

    EBFaceFlag () = default;
    EBFaceFlag (const EBFaceFlag& rhs) noexcept = default;
    EBFaceFlag& operator= (const EBFaceFlag& rhs) noexcept = default;
    EBFaceFlag& operator= (uint32_t i) { flag = i; return *this; }

    void setRegular () {
        flag &= zero_lower_mask; // clean lower bits
        flag |= regular_bits;
        flag |= single_vof_bits;
    }

    void setCovered () {
        flag &= zero_lower_mask;
        flag |= covered_bits;
    }

    void setSingleValued () {
        flag &= zero_lower_mask;
        flag |= single_vauled_bits;
        flag |= single_vof_bits;
    }

    void setMultiValued (int n) {
        flag &= zero_lower_mask;
        flag |= multi_valued_bits;
        BL_ASSERT(n >= 2 && n <= 7);
        flag |= static_cast<uint32_t>(n) << pos_numvofs;
    }

    int getNumVoFs () const {
        return (flag & one_numvofs_mask) >> pos_numvofs;
    }

    bool isRegular () const {
        return (flag & one_type_mask) == regular_bits;
    }

    bool isSingleValued () const {
        return (flag & one_type_mask) == single_vauled_bits;
    }

    bool isMultiValued () const {
        return (flag & one_type_mask) == multi_valued_bits;
    }

    bool isCovered () const {
        return (flag & one_type_mask) == covered_bits;
    }

    void setDisconnected () {
        flag &= one_lower_mask;
    }

    void setDisconnected (int dir, const IntVect& iv) {
        int i, j;
#if (AMREX_SPACEDIM == 2)
        if (dir == 0) {
            i = iv[1];
            j = 0;
        } else {
            i = iv[0];
            j = 0;
        }
#else
        if (dir == 0) {
            i = iv[1];
            j = iv[2];
        } else if (dir == 1) {
            i = iv[0];
            j = iv[2];                
        } else {
            i = iv[0];
            j = iv[1];
        }
#endif
        const int n = pos_ngbr[j+1][i+1];
        flag &= ~(1 << n);
    }

    void setConnected () {
        flag |= zero_lower_mask;
    }

    void setConnected (int dir, const IntVect& iv) {
        int i, j;
#if (AMREX_SPACEDIM == 2)
        if (dir == 0) {
            i = iv[1];
            j = 0;
        } else {
            i = iv[0];
            j = 0;
        }
#else
        if (dir == 0) {
            i = iv[1];
            j = iv[2];
        } else if (dir == 1) {
            i = iv[0];
            j = iv[2];                
        } else {
            i = iv[0];
            j = iv[1];
        }
#endif
        const int n = pos_ngbr[j+1][i+1];
        flag |= 1 << n;
    }

private:

    // masks for lowest 5 bits (i.e., bit number 0-4)
    static constexpr uint32_t one_lower_mask      =  0x1f;
    static constexpr uint32_t zero_lower_mask     = ~0x1f;

    // masks lowest 2 bits (i.e., bit number 0-1)
    static constexpr uint32_t one_type_mask       =  0x3;
    static constexpr uint32_t zero_type_mask      = ~0x3;

    // masks for bit number 2-4
    static constexpr uint32_t one_numvofs_mask    =  0x1c;
    static constexpr uint32_t zero_numvofs_mask   = ~0x1c;

    // these represent face types
    static constexpr uint32_t regular_bits        =  0x0;
    static constexpr uint32_t single_vauled_bits  =  0x1;
    static constexpr uint32_t multi_valued_bits   =  0x2;
    static constexpr uint32_t covered_bits        =  0x3;

    // this represent single vof (regular is considered as single vof too)
    static constexpr uint32_t single_vof_bits     =  0x4;

    // There are 32 bits.  The lowest 2 bits are used for face type:
    // regular, single-valued, multi-valued, and covered.  The next 3
    // bits are for the number of vofs.  Nine of the rest 27 bits are
    // used for connectivity with neighbors.

    static constexpr int w_lower_mask = 5;
    static constexpr int w_type       = 2;
    static constexpr int w_numvofs    = 3;
    static constexpr int pos_numvofs  = 2;
    static constexpr std::array<std::array<int,3>,3> pos_ngbr
        {{ std::array<int,3>{{ w_lower_mask  , w_lower_mask+1, w_lower_mask+2 }},
           std::array<int,3>{{ w_lower_mask+3, w_lower_mask+4, w_lower_mask+5 }},
           std::array<int,3>{{ w_lower_mask+6, w_lower_mask+7, w_lower_mask+8 }} }};

    static constexpr uint32_t default_value = zero_lower_mask; // regular connected with all neighbors 

    uint32_t flag = default_value;
};


class EBFaceFlagFab
{
public:
    
    EBFaceFlagFab (const Box& bx, int, bool, bool);
    ~EBFaceFlagFab ();

    EBFaceFlagFab (const EBFaceFlagFab& rhs) = delete;
    EBFaceFlagFab& operator= (const EBFaceFlagFab& rhs) = delete;
    EBFaceFlagFab& operator= (EBFaceFlagFab&& rhs) = delete;
    
    EBFaceFlagFab (EBFaceFlagFab&& rhs) noexcept = default;

    Box box () const { return m_box; }

    BaseFab      <EBFaceFlag>& getFaceFlagFab (int dir)       { return m_fabs[dir]; }
    const BaseFab<EBFaceFlag>& getFaceFlagFab (int dir) const { return m_fabs[dir]; }

    const EBFaceFlag& operator() (int dir, const IntVect& iv) const { return m_fabs[dir](iv); }
    EBFaceFlag&       operator() (int dir, const IntVect& iv)       { return m_fabs[dir](iv); }

private:
    Box m_box;
    std::array<BaseFab<EBFaceFlag>,AMREX_SPACEDIM> m_fabs;
};

}

#endif
