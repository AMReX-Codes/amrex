
#ifndef BL_ORIENTATION_H
#define BL_ORIENTATION_H
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_GpuQualifiers.H>

#include <iosfwd>

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
    constexpr Orientation () noexcept = default;
    //! Set the orientation of a side.
    AMREX_GPU_HOST_DEVICE
    Orientation (int dir, Side side) noexcept
        :
        val(AMREX_SPACEDIM*side + dir)
    {
        BL_ASSERT(0 <= dir && dir < AMREX_SPACEDIM);
    }
    AMREX_GPU_HOST_DEVICE
    constexpr Orientation (Direction dir, Side side) noexcept
        : val(AMREX_SPACEDIM*side + static_cast<int>(dir))
    {}

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
    constexpr operator int () const noexcept { return val; }
    //! Return opposite orientation.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Orientation flip () const noexcept
    {
        return Orientation(val < AMREX_SPACEDIM ? val+AMREX_SPACEDIM : val-AMREX_SPACEDIM);
    }
    //! Returns the coordinate direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int coordDir () const noexcept { return val%AMREX_SPACEDIM; }
    //! Returns the orientation of the face -- low or high.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Side faceDir () const noexcept { return Side(val/AMREX_SPACEDIM); }
    //! Returns true if Orientation is low.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool isLow () const noexcept { return val < AMREX_SPACEDIM; }
    //! Returns true if Orientation is high.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool isHigh () const noexcept { return val >= AMREX_SPACEDIM; }
    //! Read from an istream.
    friend std::istream& operator>> (std::istream& is, Orientation& o);

    //! Int value of the x-lo-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int xlo () noexcept { return 0; }

    //! Int value of the x-hi-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int xhi () noexcept { return AMREX_SPACEDIM; }

    //! Int value of the y-lo-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int ylo () noexcept { return 1; }

    //! Int value of the y-hi-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int yhi () noexcept { return 1+AMREX_SPACEDIM; }

    //! Int value of the z-lo-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int zlo () noexcept { return 2; }

    //! Int value of the z-hi-face
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int zhi () noexcept { return 2+AMREX_SPACEDIM; }

private:
    //! Used internally.
    AMREX_GPU_HOST_DEVICE
    Orientation (int v) noexcept : val(v) {}
    //
    //! The data.
    int val = -1;
};

//! Write to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const Orientation& o);


//! An Iterator over the Orientation of Faces of a Box
class OrientationIter
{

public:
    //! The default constructor.
    constexpr OrientationIter () noexcept = default;
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
    Orientation operator() () const noexcept { BL_ASSERT(isValid()); return Orientation(face); }
    //! Cast to void*.  Used to test if iterator is valid.
    AMREX_GPU_HOST_DEVICE
    operator void* () noexcept { return 0 <= face && face < 2*AMREX_SPACEDIM ? this : nullptr; }
    //! Is the iterator valid?
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool isValid () const noexcept { return 0 <= face && face < 2*AMREX_SPACEDIM; }
    //! Pre-decrement.
    AMREX_GPU_HOST_DEVICE
    OrientationIter& operator-- () noexcept { BL_ASSERT(isValid()); --face; return *this; }
    //! Pre-increment.
    AMREX_GPU_HOST_DEVICE
    OrientationIter& operator++ () noexcept { BL_ASSERT(isValid()); ++face; return *this; }
    //! Post-decrement.
    AMREX_GPU_HOST_DEVICE
    OrientationIter operator-- (int) noexcept
    {
        BL_ASSERT(isValid()); OrientationIter it(face); --face; return it;
    }
    //! Post-increment.
    AMREX_GPU_HOST_DEVICE
    OrientationIter operator++ (int) noexcept
    {
        BL_ASSERT(isValid()); OrientationIter it(face); ++face; return it;
    }
    //! The equality operator.
    AMREX_GPU_HOST_DEVICE
    bool operator== (const OrientationIter& oi) const noexcept
    {
        BL_ASSERT(isValid() && oi.isValid()); return face == oi.face;
    }
    //! The inequality operator.
    AMREX_GPU_HOST_DEVICE
    bool operator!= (const OrientationIter& oi) const noexcept
    {
        BL_ASSERT(isValid() && oi.isValid()); return face != oi.face;
    }

private:

    int face = 0;
    //! Construct an iterator on the face.
    AMREX_GPU_HOST_DEVICE
    OrientationIter (int _face) noexcept : face(_face) {}
};

}

#endif /*BL_ORIENTATION_H*/
