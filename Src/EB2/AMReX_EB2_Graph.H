#ifndef AMREX_EB2_GRAPH_H_
#define AMREX_EB2_GRAPH_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_FabFactory.H>

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

using CellType_t = uint32_t;

namespace CellType {
    static constexpr CellType_t regular=0x0;
    static constexpr CellType_t covered=0x1;
    static constexpr CellType_t irregular=0x10;
}

struct Neighbor
{
    void setXLo (uint32_t t) { flag &= zero_xlo_mask; flag |= t; }
    void setXHi (uint32_t t) { flag &= zero_xhi_mask; flag |= (t << 2); }
    void setYLo (uint32_t t) { flag &= zero_ylo_mask; flag |= (t << 4); }
    void setYHi (uint32_t t) { flag &= zero_yhi_mask; flag |= (t << 6); }
    void setZLo (uint32_t t) { flag &= zero_zlo_mask; flag |= (t << 8); }
    void setZHi (uint32_t t) { flag &= zero_zhi_mask; flag |= (t << 10); }

    uint32_t getXLo () const { return (flag & one_xlo_mask); }
    uint32_t getXHi () const { return (flag & one_xhi_mask) >> 2; }
    uint32_t getYLo () const { return (flag & one_ylo_mask) >> 4; }
    uint32_t getYHi () const { return (flag & one_yhi_mask) >> 6; }
    uint32_t getZLo () const { return (flag & one_zlo_mask) >> 8; }
    uint32_t getZHi () const { return (flag & one_zhi_mask) >> 10; }

private:

    uint32_t flag = 0;

    static constexpr uint32_t  one_xlo_mask =  0x3;
    static constexpr uint32_t zero_xlo_mask = ~0x3;
    static constexpr uint32_t  one_xhi_mask =  0xC;
    static constexpr uint32_t zero_xhi_mask = ~0xC;
    static constexpr uint32_t  one_ylo_mask =  0x30;
    static constexpr uint32_t zero_ylo_mask = ~0x30;
    static constexpr uint32_t  one_yhi_mask =  0xC0;
    static constexpr uint32_t zero_yhi_mask = ~0xC0;
    static constexpr uint32_t  one_zlo_mask =  0x300;
    static constexpr uint32_t zero_zlo_mask = ~0x300;
    static constexpr uint32_t  one_zhi_mask =  0xC00;
    static constexpr uint32_t zero_zhi_mask = ~0xC00;
};

struct Vertex
{
    Cell m_cell;
    Neighbor m_nghb;
    Vector<Edge> m_edges;
};


class Graph
{
public:

    Graph () = default;
    Graph (const Box& box, int, bool, bool)
        : m_box(box)
        {}

private:

    Box m_box;
    std::unordered_map<IntVect, Vector<Vertex>, IntVect::shift_hasher> m_data;
};

//template

}}

#endif
