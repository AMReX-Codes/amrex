#ifndef AMREX_EBCELLFLAG_H_
#define AMREX_EBCELLFLAG_H_

#include <array>
#include <cstdint>
#include <AMReX_IntVect.H>

namespace amrex {

class EBCellFlag
{
public:

    EBCellFlag () = default;
    EBCellFlag (const EBCellFlag& rhs) noexcept = default;
    EBCellFlag& operator= (const EBCellFlag& rhs) noexcept = default;
    EBCellFlag& operator= (uint32_t i) { flag = i; return *this; }

    void setRegular () {
        flag &= zero_lower_bits; // clean lower bits
        flag |= regular_bits;
    }

    void setCovered () {
        flag &= zero_lower_bits;
        flag |= covered_bits;
    }

    void setSingleValued () {
        flag &= zero_lower_bits;
        flag |= single_vauled_bits;
    }

    void setMultiValued () {
        flag &= zero_lower_bits;
        flag |= multi_valued_bits;
    }

    bool isRegular () const {
        return (flag & one_lower_bits) == regular_bits;
    }

    bool isSingleValued () const {
        return (flag & one_lower_bits) == single_vauled_bits;
    }

    bool isMultiValued () const {
        return (flag & one_lower_bits) == multi_valued_bits;
    }

    bool isCovered () const {
        return (flag & one_lower_bits) == covered_bits;
    }

    void setDisconnected () {
        flag &= one_lower_bits;
    }

    void setDisconnected (const IntVect& iv) {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = pos[k][j][i];
        flag &= ~(1 << n);
    }

    void setConnected () {
        flag |= zero_lower_bits;
    }

    void setConnected (const IntVect& iv) {
        int i=0, j=0, k=0;
        AMREX_D_TERM(i=iv[0];, j=iv[1];, k=iv[2]);
        const int n = pos[k+1][j+1][i+1];
        flag |= 1 << n;
    }

private:

    static constexpr int width = 4;
    static constexpr std::array<std::array<std::array<int,3>,3>,3> pos
        {{ std::array<std::array<int,3>,3>{{{ width   , width+ 1, width+ 2 },
                                            { width+ 3, width+ 4, width+ 5 },
                                            { width+ 6, width+ 7, width+ 8 }}},
           std::array<std::array<int,3>,3>{{{ width+ 9, width+10, width+11 },
                                            { width+12, width+13, width+14 },
                                            { width+15, width+16, width+17 }}},
           std::array<std::array<int,3>,3>{{{ width+18, width+19, width+20 },
                                            { width+21, width+22, width+23 },
                                            { width+24, width+25, width+26 }}} }};

    static constexpr uint32_t one_lower_bits      =  0xf;
    static constexpr uint32_t zero_lower_bits     = ~0xf;
    static constexpr uint32_t regular_bits        =  0x0;
    static constexpr uint32_t single_vauled_bits  =  0x1;
    static constexpr uint32_t multi_valued_bits   =  0x2;
    static constexpr uint32_t covered_bits        =  0x3;

    static constexpr uint32_t default_value = zero_lower_bits; // regular connected with all neighbors    
    
    uint32_t flag = default_value;
};

}

#endif
