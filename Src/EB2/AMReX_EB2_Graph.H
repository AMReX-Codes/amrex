#ifndef AMREX_EB2_GRAPH_H_
#define AMREX_EB2_GRAPH_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>

#include <unordered_map>
#include <utility>
#include <cstdint>

namespace amrex { namespace EB2 {

struct Cell
{
    RealArray m_volfrac;
    // and more such as boundary centroid and normal
};

struct Face
{
    RealArray m_facecent;
    Real m_areafrac;
};

using VertexID = std::pair<IntVect,int>;

struct Edge
{
    Face m_face;
    VertexID m_vid;
};

using Type_t = uint32_t;

namespace Type {
    static constexpr Type_t regular=0x0;
    static constexpr Type_t covered=0x1;
    static constexpr Type_t irregular=0x2;
    static constexpr Type_t unknown=0x3;
}

struct Neighbor
{
    void setXLo (Type_t t) { flag &= zero_xlo_mask; flag |= t; }
    void setXHi (Type_t t) { flag &= zero_xhi_mask; flag |= (t << 2); }
    void setYLo (Type_t t) { flag &= zero_ylo_mask; flag |= (t << 4); }
    void setYHi (Type_t t) { flag &= zero_yhi_mask; flag |= (t << 6); }
    void setZLo (Type_t t) { flag &= zero_zlo_mask; flag |= (t << 8); }
    void setZHi (Type_t t) { flag &= zero_zhi_mask; flag |= (t << 10); }

    Type_t getXLo () const { return (flag & one_xlo_mask); }
    Type_t getXHi () const { return (flag & one_xhi_mask) >> 2; }
    Type_t getYLo () const { return (flag & one_ylo_mask) >> 4; }
    Type_t getYHi () const { return (flag & one_yhi_mask) >> 6; }
    Type_t getZLo () const { return (flag & one_zlo_mask) >> 8; }
    Type_t getZHi () const { return (flag & one_zhi_mask) >> 10; }

private:

    Type_t flag = 0;

    static constexpr Type_t  one_xlo_mask =  0x3;
    static constexpr Type_t zero_xlo_mask = ~0x3;
    static constexpr Type_t  one_xhi_mask =  0xC;
    static constexpr Type_t zero_xhi_mask = ~0xC;
    static constexpr Type_t  one_ylo_mask =  0x30;
    static constexpr Type_t zero_ylo_mask = ~0x30;
    static constexpr Type_t  one_yhi_mask =  0xC0;
    static constexpr Type_t zero_yhi_mask = ~0xC0;
    static constexpr Type_t  one_zlo_mask =  0x300;
    static constexpr Type_t zero_zlo_mask = ~0x300;
    static constexpr Type_t  one_zhi_mask =  0xC00;
    static constexpr Type_t zero_zhi_mask = ~0xC00;
};

struct Vertex
{
    Cell m_cell;
    Neighbor m_nghb;
    Vector<Edge> m_edges;
};

using Graph = std::unordered_map<IntVect, Vector<Vertex>, IntVect::hasher>;

}}

#endif
