
#ifndef AMREX_COORDSYS_H_
#define AMREX_COORDSYS_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>

#include <limits>

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

    //! Is ok?
    [[nodiscard]] bool Ok () const noexcept { return ok; }

    //! Set the CoordType.
    void SetCoord (CoordType coord) noexcept { c_sys = coord; }

    //! Returns the CoordType.
    [[nodiscard]] CoordType Coord () const noexcept { return c_sys; }

    //! Returns the CoordType as an int.
    [[nodiscard]] int CoordInt () const  noexcept { return static_cast<int>(c_sys); }

    //! Is CoordType == SPHERICAL?
    [[nodiscard]] bool IsSPHERICAL () const noexcept {
        BL_ASSERT(c_sys != undef); return (c_sys == SPHERICAL);
    }

    //! Is CoordType == RZ?
    [[nodiscard]] bool IsRZ () const noexcept {
        BL_ASSERT(c_sys != undef); return (c_sys == RZ);
    }

    //! Is CoordType == cartesion?
    [[nodiscard]] bool IsCartesian () const noexcept {
        BL_ASSERT(c_sys != undef); return (c_sys == cartesian);
    }

    //! Sets the offset for each coordinate direction.
    void SetOffset (const Real* x_lo) noexcept;

    //! Returns the offset.
    [[nodiscard]] const Real* Offset () const noexcept { return offset; }

    //! Returns the offset for the specified coordinate direction.
    [[nodiscard]] Real Offset (int dir) const noexcept { return offset[dir]; }

    //! Returns the cellsize for each coordinate direction.
    [[nodiscard]] const Real* CellSize () const noexcept { BL_ASSERT(ok); return dx; }

    //! Returns the cellsize for the specified coordinate direction.
    [[nodiscard]] Real CellSize (int dir) const noexcept { BL_ASSERT(ok); return dx[dir]; }

    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> CellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return {{ AMREX_D_DECL(dx[0],dx[1],dx[2]) }};
    }

    //! Returns the inverse cellsize for each coordinate direction.
    [[nodiscard]] const Real* InvCellSize () const noexcept { BL_ASSERT(ok); return inv_dx; }

    //! Returns the inverse cellsize for the specified coordinate direction.
    [[nodiscard]] Real InvCellSize (int dir) const noexcept { BL_ASSERT(ok); return inv_dx[dir]; }

    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> InvCellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return {{ AMREX_D_DECL(inv_dx[0],inv_dx[1],inv_dx[2]) }};
    }

    //! Returns location of cell center in specified direction.
    [[nodiscard]] Real CellCenter (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*((Real)0.5+ (Real)point);
    }

    //! Return location of cell center.
    void CellCenter (const IntVect& point, Vector<Real>& loc) const noexcept;

    //! Return location of cell center.
    void CellCenter (const IntVect& point, Real* loc) const noexcept;

    //! Returns location of lo edge in specified direction.
    [[nodiscard]] Real LoEdge (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*static_cast<Real>(point);
    }

    //! Equivalent to LoEdge(point[dir], dir).
    [[nodiscard]] Real LoEdge (const IntVect& point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*static_cast<Real>(point[dir]);
    }

    //! Returns location of hi edge in specified direction.
    [[nodiscard]] Real HiEdge (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*static_cast<Real>(point + 1);
    }

    //! Equivalent to HiEdge(point[dir], dir).
    [[nodiscard]] Real HiEdge (const IntVect& point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*static_cast<Real>(point[dir] + 1);
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
    void SetVolume (FArrayBox& a_volfab, const Box& region) const;
    /**
    * \brief Compute cell volumes in given region and place them into
    * resize()d input FAB.
    */
    void GetVolume (FArrayBox& vol, const Box& region) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * place them into input FAB.
    */
    void SetDLogA (FArrayBox& a_dlogafab, const Box& region, int dir) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * return the results in the resize()d input FAB.
    */
    void GetDLogA (FArrayBox& dloga, const Box& region, int dir) const;

    //! Return the volume of the specified cell.
    [[nodiscard]] Real Volume (const IntVect& point) const;

    //! Return the volume of the specified cell.
    [[nodiscard]] Real Volume (const Real xlo[AMREX_SPACEDIM],
                               const Real xhi[AMREX_SPACEDIM]) const;

    /**
    * \brief Compute area of cell faces in given region and given
    * index direction and return the result in input FAB.
    */
    void SetFaceArea (FArrayBox& a_areafab, const Box& region, int dir) const;

    /**
    * \brief Compute area of cell faces in given region and given
    * index direction and return the result in resize()d input FAB.
    */
    void GetFaceArea (FArrayBox& area, const Box& region, int dir) const;

    //! Returns lo face area of given cell in direction dir.
    [[nodiscard]] Real AreaLo (const IntVect& point, int dir) const noexcept;

    //! Returns hi face area of given cell in direction dir.
    [[nodiscard]] Real AreaHi (const IntVect& point, int dir) const noexcept;

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
    void GetEdgeVolCoord (Vector<Real>& vc, const Box& region, int dir) const;

    /**
    * \brief Return array of volume coordinates at cell centers
    * for region in given direction.
    */
    void GetCellVolCoord (Vector<Real>& vc, const Box& region, int dir) const;

protected:
    // c_sys and offset used to be static
    CoordType c_sys = undef;
    Real      offset[AMREX_SPACEDIM];

    Real dx[AMREX_SPACEDIM] = {AMREX_D_DECL(0.,0.,0.)};
    Real inv_dx[AMREX_SPACEDIM]
        = {AMREX_D_DECL(std::numeric_limits<Real>::infinity(),
                        std::numeric_limits<Real>::infinity(),
                        std::numeric_limits<Real>::infinity())};
    bool ok = false;
};

}

#endif /*_COORDSYS_H_*/
