
#ifndef AMREX_COORDSYS_H_
#define AMREX_COORDSYS_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>

namespace amrex {

class FArrayBox;

/**
* \brief Coordinate System
*
* Routines for mapping between physical coordinate system and index space.
*/

class CoordSys
{
public:

    enum CoordType { undef = -1, cartesian = 0, RZ = 1, SPHERICAL = 2 };
    //! Nice ASCII output.
    friend std::ostream& operator<< (std::ostream&, const CoordSys& );

    //! Nice ASCII input.
    friend std::istream& operator>> (std::istream&, CoordSys& );

    //! Default constructor to undefined state.
    CoordSys ();

    //! Construct, specify a cell size, and low end of index range.
    explicit CoordSys (const Real* cell_dx);

    ~CoordSys () = default;
    CoordSys (const CoordSys& rhs) = default;
    CoordSys (CoordSys&& rhs) noexcept = default;
    CoordSys& operator= (const CoordSys& rhs) = default;
    CoordSys& operator= (CoordSys&& rhs) = default;

    //! Initialize after construction.
    void define (const Real* cell_dx);

    //! Set the CoordType.
    static void SetCoord (CoordType coord) { c_sys = coord; }

    //! Returns the CoordType.
    static CoordType Coord () noexcept { return c_sys; }

    //! Returns the CoordType as an int.
    static int CoordInt () noexcept;

    //! Is CoordType == SPHERICAL?
    static bool IsSPHERICAL () noexcept;

    //! Is CoordType == RZ?
    static bool IsRZ () noexcept;

    //! Is CoordType == cartesion?
    static bool IsCartesian () noexcept { BL_ASSERT(c_sys != undef); return (c_sys == cartesian); }

    //! Sets the offset for each coordinate direction.
    static void SetOffset (const Real* x_lo);

    //! Returns the offset.
    static const Real* Offset () noexcept { return offset; }

    //! Returns the offset for the specified coordinate direction.
    static Real Offset (int dir) noexcept { return offset[dir]; }

    //! Returns the cellsize for each coordinate direction.
    const Real* CellSize () const noexcept { BL_ASSERT(ok); return dx; }

    //! Returns the cellsize for the specified coordinate direction.
    Real CellSize (int dir) const noexcept { BL_ASSERT(ok); return dx[dir]; }

    GpuArray<Real,AMREX_SPACEDIM> CellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return { AMREX_D_DECL(dx[0],dx[1],dx[2]) };
    }

    //! Returns the inverse cellsize for each coordinate direction.
    const Real* InvCellSize () const noexcept { BL_ASSERT(ok); return inv_dx; }

    //! Returns the inverse cellsize for the specified coordinate direction.
    Real InvCellSize (int dir) const noexcept { BL_ASSERT(ok); return inv_dx[dir]; }

    GpuArray<Real,AMREX_SPACEDIM> InvCellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return { AMREX_D_DECL(inv_dx[0],inv_dx[1],inv_dx[2]) };
    }

    //! Returns location of cell center in specified direction.
    Real CellCenter (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*(0.5+ (Real)point);
    }

    //! Return location of cell center.
    void CellCenter (const IntVect& point, Vector<Real>& loc) const noexcept;

    //! Return location of cell center.
    void CellCenter (const IntVect& point, Real* loc) const noexcept;

    //! Returns location of lo edge in specified direction.
    Real LoEdge (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*point;
    }

    //! Equivalent to LoEdge(point[dir], dir).
    Real LoEdge (const IntVect& point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*point[dir];
    }

    //! Returns location of hi edge in specified direction.
    Real HiEdge (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*(point + 1);
    }

    //! Equivalent to HiEdge(point[dir], dir).
    Real HiEdge (const IntVect& point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*(point[dir] + 1);
    }

    //! Sets location of lo face into loc.
    void LoFace (const IntVect& point, int dir, Vector<Real>& loc) const noexcept;

    //! Sets location of lo face into loc.
    void LoFace (const IntVect& point, int dir, Real* loc) const noexcept;

    //! Sets location of hi face into loc.
    void HiFace (const IntVect& point, int dir, Vector<Real>& loc) const noexcept;

    //! Sets location of hi face into loc.
    void HiFace (const IntVect& point, int dir, Real* loc) const noexcept;

    //! Return location of lower left hand corner.
    void LoNode (const IntVect& point, Vector<Real>& loc) const noexcept;

    //! Return location of lower left hand corner.
    void LoNode (const IntVect& point, Real* loc) const noexcept;

    //! Return location of upper right hand corner.
    void HiNode (const IntVect& point, Vector<Real>& loc) const noexcept;

    //! Return location of upper right hand corner.
    void HiNode (const IntVect& point, Real* loc) const noexcept;
    /**
    * \brief Returns cell centered index of cell containing point.
    * This may return undesired results if point
    * is on a cell boundary.
    */
    IntVect CellIndex (const Real* point) const noexcept;
    /**
    * \brief Returns node centered index of lower left hand corner of
    * cell containing this point.
    */
    IntVect LowerIndex (const Real* point) const noexcept;
    /**
    * \brief Returns node centered index of upper right hand corner of
    * cell containing this point.
    */
    IntVect UpperIndex (const Real* point) const noexcept;
    /**
    * \brief Compute cell volumes in given region and place them into
    * input FAB.
    */
    void SetVolume (FArrayBox& vol, const Box& region) const;
    /**
    * \brief Compute cell volumes in given region and place them into
    * resize()d input FAB.
    */
    void GetVolume (FArrayBox& vol, const Box& region) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * place them into input FAB.
    */
    void SetDLogA (FArrayBox& dloga, const Box& region, int dir) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * return the results in the resize()d input FAB.
    */
    void GetDLogA (FArrayBox& dloga, const Box& region, int dir) const;

    //! Return the volume of the specified cell.
    Real Volume (const IntVect& point) const;

    //! Return the volume of the specified cell.
    Real Volume (const Real xlo[AMREX_SPACEDIM],
                 const Real xhi[AMREX_SPACEDIM]) const;

    /**
    * \brief Compute area of cell faces in given region and given
    * index direction and return the result in input FAB.
    */
    void SetFaceArea (FArrayBox& area, const Box& region, int dir) const;

    /**
    * \brief Compute area of cell faces in given region and given
    * index direction and return the result in resize()d input FAB.
    */
    void GetFaceArea (FArrayBox& area, const Box& region, int dir) const;

    //! Returns lo face area of given cell in direction dir.
    Real AreaLo (const IntVect& point, int dir) const noexcept;

    //! Returns hi face area of given cell in direction dir.
    Real AreaHi (const IntVect& point, int dir) const noexcept;

    /**
    * \brief Return array of physical locations of cell edges
    * in the resize()d input array.
    */
    void GetEdgeLoc (Vector<Real>& loc, const Box& region, int dir) const;

    /**
    * \brief Return array of physical locations of cell centers
    * in the resize()d input array.
    */
    void GetCellLoc (Vector<Real>& loc, const Box& region, int dir) const;

    /**
    * \brief Return array of volume coordinates at cell edges
    * for region in given direction.
    */
    void GetEdgeVolCoord (Vector<Real>& loc, const Box& region, int dir) const;

    /**
    * \brief Return array of volume coordinates at cell centers
    * for region in given direction.
    */
    void GetCellVolCoord (Vector<Real>& loc, const Box& region, int dir) const;

protected:
    //! Static data.
    static CoordType c_sys;
    static Real      offset[AMREX_SPACEDIM];

    //! Non-static data.
    Real dx[AMREX_SPACEDIM];
    Real inv_dx[AMREX_SPACEDIM];
    bool ok;
};

inline
bool
CoordSys::IsSPHERICAL () noexcept
{
    BL_ASSERT(c_sys != undef);
#if (AMREX_SPACEDIM <= 2)
    return (c_sys == SPHERICAL);
#endif
#if (AMREX_SPACEDIM == 3)
    return false;
#endif
}

inline
bool
CoordSys::IsRZ () noexcept
{
    BL_ASSERT(c_sys != undef);
#if (AMREX_SPACEDIM <= 2)
    return (c_sys == RZ);
#endif
#if (AMREX_SPACEDIM == 3)
    return false;
#endif
}

}

#endif /*_COORDSYS_H_*/
