#ifndef AMREX_EB2_GRAPH_H_
#define AMREX_EB2_GRAPH_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>

#include <unordered_map>
#include <utility>
#include <cstdint>

/**
 * \file AMReX_EB2_Graph.H
 *
 * Lightweight types that store EB irregular-cell geometry (volume fractions,
 * faces, neighbors, etc.) inside EB2 graphs.
 */

namespace amrex::EB2 {

//! Per-cell geometric information stored in EB graphs.
struct Cell
{
    RealArray m_volfrac;
    // and more such as boundary centroid and normal
};

//! Metadata associated with an EB face.
struct Face
{
    RealArray m_facecent;
    Real m_areafrac;
};

using VertexID = std::pair<IntVect,int>;

//! Edge connecting a graph vertex to a face and downstream vertex id.
struct Edge
{
    Face m_face;
    VertexID m_vid;
};

using Type_t = uint32_t;

//! Type flags representing regular/covered/etc.
namespace Type {
    static constexpr Type_t regular=0x0;
    static constexpr Type_t covered=0x1;
    static constexpr Type_t irregular=0x2;
    static constexpr Type_t unknown=0x3;
}

//! Packed neighbor structure that tracks the neighbor type in every direction.
struct Neighbor
{
    //! Mark the low-x neighbor type.
    void setXLo (Type_t t) { flag &= zero_xlo_mask; flag |= t; }
    //! Mark the high-x neighbor type.
    void setXHi (Type_t t) { flag &= zero_xhi_mask; flag |= (t << 2); }
    //! Mark the low-y neighbor type.
    void setYLo (Type_t t) { flag &= zero_ylo_mask; flag |= (t << 4); }
    //! Mark the high-y neighbor type.
    void setYHi (Type_t t) { flag &= zero_yhi_mask; flag |= (t << 6); }
    //! Mark the low-z neighbor type.
    void setZLo (Type_t t) { flag &= zero_zlo_mask; flag |= (t << 8); }
    //! Mark the high-z neighbor type.
    void setZHi (Type_t t) { flag &= zero_zhi_mask; flag |= (t << 10); }

    [[nodiscard]] Type_t getXLo () const { return (flag & one_xlo_mask); }
    [[nodiscard]] Type_t getXHi () const { return (flag & one_xhi_mask) >> 2; }
    [[nodiscard]] Type_t getYLo () const { return (flag & one_ylo_mask) >> 4; }
    [[nodiscard]] Type_t getYHi () const { return (flag & one_yhi_mask) >> 6; }
    [[nodiscard]] Type_t getZLo () const { return (flag & one_zlo_mask) >> 8; }
    [[nodiscard]] Type_t getZHi () const { return (flag & one_zhi_mask) >> 10; }

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

//! Vertex capturing one irregular cell plus its neighbor connectivity.
struct Vertex
{
    Cell m_cell;
    Neighbor m_nghb;
    Vector<Edge> m_edges;
};

//! EB graph type: map from cell index to all EB vertices stored there.
using Graph = std::unordered_map<IntVect, Vector<Vertex>, IntVect::hasher>;

}

#endif
