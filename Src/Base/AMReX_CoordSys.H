
#ifndef AMREX_COORDSYS_H_
#define AMREX_COORDSYS_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_CoordSys.H
 * \brief Coordinate system that maps between physical space and index space.
 */

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>

#include <limits>

namespace amrex {

class FArrayBox;

/**
 * \class amrex::CoordSys
 * \ingroup amrex_real_space
 * \brief Coordinate system; provides mapping between physical space and index space.
 *
 * Supports Cartesian, cylindrical (RZ), and spherical coordinate systems.
 */
class CoordSys
{
public:

    //! Coordinate system type: undefined, Cartesian, cylindrical (RZ), or spherical.
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

    //! Is CoordType == cartesian?
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

    //! Returns the cell sizes as a GpuArray for use on host or device.
    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> CellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return {{ AMREX_D_DECL(dx[0],dx[1],dx[2]) }};
    }

    //! Returns the inverse cellsize for each coordinate direction.
    [[nodiscard]] const Real* InvCellSize () const noexcept { BL_ASSERT(ok); return inv_dx; }

    //! Returns the inverse cellsize for the specified coordinate direction.
    [[nodiscard]] Real InvCellSize (int dir) const noexcept { BL_ASSERT(ok); return inv_dx[dir]; }

    //! Returns the inverse cell sizes as a GpuArray for use on host or device.
    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> InvCellSizeArray () const noexcept {
        BL_ASSERT(ok);
        return {{ AMREX_D_DECL(inv_dx[0],inv_dx[1],inv_dx[2]) }};
    }

    //! Returns location of cell center in specified direction.
    [[nodiscard]] Real CellCenter (int point, int dir) const noexcept
    {
        BL_ASSERT(ok); return offset[dir] + dx[dir]*((Real)0.5+ (Real)point);
    }

    /**
     * \brief Return location of the specified cell center.
     * \param point Cell-centered index (i,j,k).
     * \param loc Output vector resized to \c AMREX_SPACEDIM holding coordinates.
     */
    void CellCenter (const IntVect& point, Vector<Real>& loc) const noexcept;

    /**
     * \brief Return location of the specified cell center.
     * \param point Cell-centered index (i,j,k).
     * \param loc Output pointer to an array of length \c AMREX_SPACEDIM.
     */
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

    /**
     * \brief Compute the physical location of the lower face in direction \p dir.
     * \param point Cell index of the owning cell.
     * \param dir Coordinate direction (0-based).
     * \param loc Output vector resized to \c AMREX_SPACEDIM holding coordinates.
     */
    void LoFace (const IntVect& point, int dir, Vector<Real>& loc) const noexcept;

    /**
     * \brief Compute the physical location of the lower face in direction \p dir.
     * \param point Cell index of the owning cell.
     * \param dir Coordinate direction (0-based).
     * \param loc Output pointer to an array of length \c AMREX_SPACEDIM.
     */
    void LoFace (const IntVect& point, int dir, Real* loc) const noexcept;

    /**
     * \brief Compute the physical location of the upper face in direction \p dir.
     * \param point Cell index of the owning cell.
     * \param dir Coordinate direction (0-based).
     * \param loc Output vector resized to \c AMREX_SPACEDIM holding coordinates.
     */
    void HiFace (const IntVect& point, int dir, Vector<Real>& loc) const noexcept;

    /**
     * \brief Compute the physical location of the upper face in direction \p dir.
     * \param point Cell index of the owning cell.
     * \param dir Coordinate direction (0-based).
     * \param loc Output pointer to an array of length \c AMREX_SPACEDIM.
     */
    void HiFace (const IntVect& point, int dir, Real* loc) const noexcept;

    /**
     * \brief Return the location of the lower-left node of the cell.
     * \param point Cell index of the owning cell.
     * \param loc Output vector resized to \c AMREX_SPACEDIM holding coordinates.
     */
    void LoNode (const IntVect& point, Vector<Real>& loc) const noexcept;

    /**
     * \brief Return the location of the lower-left node of the cell.
     * \param point Cell index of the owning cell.
     * \param loc Output pointer to an array of length \c AMREX_SPACEDIM.
     */
    void LoNode (const IntVect& point, Real* loc) const noexcept;

    /**
     * \brief Return the location of the upper-right node of the cell.
     * \param point Cell index of the owning cell.
     * \param loc Output vector resized to \c AMREX_SPACEDIM holding coordinates.
     */
    void HiNode (const IntVect& point, Vector<Real>& loc) const noexcept;

    /**
     * \brief Return the location of the upper-right node of the cell.
     * \param point Cell index of the owning cell.
     * \param loc Output pointer to an array of length \c AMREX_SPACEDIM.
     */
    void HiNode (const IntVect& point, Real* loc) const noexcept;
    /**
    * \brief Returns cell-centered index of cell containing \p point.
    *
    * This may return undesired results if point
    * is on a cell boundary.
    *
    * \param point Physical-space coordinates.
    */
    IntVect CellIndex (const Real* point) const noexcept;
    /**
    * \brief Returns nodal index of lower-left corner of
    * cell containing this point.
    *
    * \param point Physical-space coordinates.
    */
    IntVect LowerIndex (const Real* point) const noexcept;
    /**
    * \brief Returns nodal index of upper right hand corner of
    * cell containing this point.
    *
    * \param point Physical-space coordinates.
    */
    IntVect UpperIndex (const Real* point) const noexcept;
    /**
    * \brief Compute cell volumes in given region and place them into
    * input FAB.
    *
    * \param a_volfab Output fab; must be sized to cover \p region (not resized by this call).
    * \param region Region of cells over which to compute volumes.
    */
    void SetVolume (FArrayBox& a_volfab, const Box& region) const;
    /**
    * \brief Compute cell volumes in given region and place them into
    * resize()d input FAB.
    *
    * \param vol Output fab; resized to \p region before being filled with volumes.
    * \param region Region of cells over which to compute volumes.
    */
    void GetVolume (FArrayBox& vol, const Box& region) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * place them into input FAB.
    *
    * \param a_dlogafab Output fab; must be sized to cover \p region (not resized by this call).
    * \param region Region of cells over which to compute the derivative.
    * \param dir Coordinate direction of the area change.
    */
    void SetDLogA (FArrayBox& a_dlogafab, const Box& region, int dir) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    * return the results in the resize()d input FAB.
    *
    * \param dloga Output fab; resized to \p region before being filled with derivatives.
    * \param region Region of cells over which to compute the derivative.
    * \param dir Coordinate direction of the area change.
    */
    void GetDLogA (FArrayBox& dloga, const Box& region, int dir) const;

    /**
     * \brief Return the physical volume of the specified cell.
     * \param point Cell-centered index (i,j,k).
     */
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
        = {AMREX_D_DECL(std::numeric_limits<Real>::max(),
                        std::numeric_limits<Real>::max(),
                        std::numeric_limits<Real>::max())};
    bool ok = false;
};

}

#endif /*_COORDSYS_H_*/
