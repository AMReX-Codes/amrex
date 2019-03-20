
#ifndef BL_ORIENTATION_H
#define BL_ORIENTATION_H

#include <iosfwd>

#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex {

enum class Direction : int { AMREX_D_DECL(x = 0, y = 1, z = 2) };

class OrientationIter;

/**
* \brief Encapsulation of the Orientation of the Faces of a Box
*
* This class encapsulates the orientation of the faces of a Box by
* providing an ordering of each of the faces of a Box in AMREX_SPACEDIM
* dimensions.  This allows iterating over all the faces of a Box.  The
* ordering first traverses the AMREX_SPACEDIM low sides from direction 0 ..
* AMREX_SPACEDIM-1 and then the AMREX_SPACEDIM high sides from direction 0 ..
* AMREX_SPACEDIM-1.
*/

class Orientation
{
public:

    friend class OrientationIter;
    //! In each dimension a face is either low or high.
    enum Side { low = 0, high = 1 };
    //! The default constructor.
    AMREX_GPU_HOST_DEVICE
    constexpr Orientation () noexcept
        :
        val(-1) {}
    //! Set the orientation of a side.
    AMREX_GPU_HOST_DEVICE
    Orientation (int dir, Side side) noexcept
        :
        val(AMREX_SPACEDIM*side + dir)
    {
        BL_ASSERT(0 <= dir && dir < AMREX_SPACEDIM);
    }
    //! Logical equality.
    AMREX_GPU_HOST_DEVICE
    bool operator== (const Orientation& o) const noexcept { return val == o.val; }
    //! Logical inequality.
    AMREX_GPU_HOST_DEVICE
    bool operator!= (const Orientation& o) const noexcept { return val != o.val; }
    //! Less-than.
    AMREX_GPU_HOST_DEVICE
    bool operator<  (const Orientation& o) const noexcept { return val < o.val; }
    //! Less-than or equal.
    AMREX_GPU_HOST_DEVICE
    bool operator<= (const Orientation& o) const noexcept { return val <= o.val; }
    //! Greater-than.
    AMREX_GPU_HOST_DEVICE
    bool operator>  (const Orientation& o) const noexcept { return val > o.val; }
    //! Greater-than or equal.
    AMREX_GPU_HOST_DEVICE
    bool operator>= (const Orientation& o) const noexcept { return val >= o.val; }
    /**
    * \brief This conversion operator maps an orientation into a
    * unique integer in the range [0 .. 2*AMREX_SPACEDIM-1]
    * according to the above ordering.
    */
    AMREX_GPU_HOST_DEVICE
    operator int () const noexcept { return val; }
    //! Return opposite orientation.
    AMREX_GPU_HOST_DEVICE
    Orientation flip () const noexcept
    {
        return Orientation(val < AMREX_SPACEDIM ? val+AMREX_SPACEDIM : val-AMREX_SPACEDIM);
    }
    //! Returns the coordinate direction.
    AMREX_GPU_HOST_DEVICE
    int coordDir () const noexcept { return val%AMREX_SPACEDIM; }
    //! Returns the orientation of the face -- low or high.
    AMREX_GPU_HOST_DEVICE
    Side faceDir () const noexcept { return Side(val/AMREX_SPACEDIM); }
    //! Returns true if Orientation is low.
    AMREX_GPU_HOST_DEVICE
    bool isLow () const noexcept { return val < AMREX_SPACEDIM; }
    //! Returns true if Orientation is high.
    AMREX_GPU_HOST_DEVICE
    bool isHigh () const noexcept { return val >= AMREX_SPACEDIM; }
    //! Read from an istream.
    friend std::istream& operator>> (std::istream& os, Orientation& o);

private:
    //! Used internally.
    AMREX_GPU_HOST_DEVICE
    Orientation (int v) noexcept : val(v) {}
    //
    //! The data.
    int val;
};

//! Write to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const Orientation& o);


//! An Iterator over the Orientation of Faces of a Box
class OrientationIter
{

public:
    //! The default constructor.
    AMREX_GPU_HOST_DEVICE
    constexpr OrientationIter () noexcept
        :
        face(0) {}
    //! Construct an iterator on the Orientation.
    AMREX_GPU_HOST_DEVICE
    OrientationIter (const Orientation& _face) noexcept
        :
        face(_face) {}
    //! Reset (rewind) the iterator.
    AMREX_GPU_HOST_DEVICE
    void rewind () noexcept { face = 0; }
    //! Return the orientation of the face.
    AMREX_GPU_HOST_DEVICE
    Orientation operator() () const noexcept { BL_ASSERT(ok()); return Orientation(face); }
    //! Cast to void*.  Used to test if iterator is valid.
    AMREX_GPU_HOST_DEVICE
    operator void* () noexcept { return 0 <= face && face < 2*AMREX_SPACEDIM ? this : 0; }
    //! Pre-decrement.
    AMREX_GPU_HOST_DEVICE
    OrientationIter& operator-- () noexcept { BL_ASSERT(ok()); --face; return *this; }
    //! Pre-increment.
    AMREX_GPU_HOST_DEVICE
    OrientationIter& operator++ () noexcept { BL_ASSERT(ok()); ++face; return *this; }
    //! Post-decrement.
    AMREX_GPU_HOST_DEVICE
    OrientationIter operator-- (int) noexcept
    {
        BL_ASSERT(ok()); OrientationIter it(face); --face; return it;
    }
    //! Post-increment.
    AMREX_GPU_HOST_DEVICE
    OrientationIter operator++ (int) noexcept
    {
        BL_ASSERT(ok()); OrientationIter it(face); ++face; return it;
    }
    //! The equality operator.
    AMREX_GPU_HOST_DEVICE
    bool operator== (const OrientationIter& oi) const noexcept
    {
        BL_ASSERT(ok() && oi.ok()); return face == oi.face;
    }
    //! The inequality operator.
    AMREX_GPU_HOST_DEVICE
    bool operator!= (const OrientationIter& oi) const noexcept
    {
        BL_ASSERT(ok() && oi.ok()); return face != oi.face;
    }

private:

    int face;
    //! Construct an iterator on the face.
    AMREX_GPU_HOST_DEVICE
    OrientationIter (int _face) noexcept : face(_face) {}
    //! Is the iterator valid?
    AMREX_GPU_HOST_DEVICE
    bool ok () const noexcept { return 0 <= face && face < 2*AMREX_SPACEDIM; }
};

}

#endif /*BL_ORIENTATION_H*/
