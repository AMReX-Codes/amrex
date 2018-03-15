
#ifndef BL_ORIENTATION_H
#define BL_ORIENTATION_H

#include <iosfwd>

#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>

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
    Orientation ()
        :
        val(-1) {}
    //! Set the orientation of a side.
    Orientation (int dir, Side side)
        :
        val(AMREX_SPACEDIM*side + dir)
    {
        BL_ASSERT(0 <= dir && dir < AMREX_SPACEDIM);
    }
    //! Logical equality.
    bool operator== (const Orientation& o) const { return val == o.val; }
    //! Logical inequality.
    bool operator!= (const Orientation& o) const { return val != o.val; }
    //! Less-than.
    bool operator<  (const Orientation& o) const { return val < o.val; }
    //! Less-than or equal.
    bool operator<= (const Orientation& o) const { return val <= o.val; }
    //! Greater-than.
    bool operator>  (const Orientation& o) const { return val > o.val; }
    //! Greater-than or equal.
    bool operator>= (const Orientation& o) const { return val >= o.val; }
    /**
    * \brief This conversion operator maps an orientation into a
    * unique integer in the range [0 .. 2*AMREX_SPACEDIM-1]
    * according to the above ordering.
    */
    operator int () const { return val; }
    //! Return opposite orientation.
    Orientation flip () const
    {
        return Orientation(val < AMREX_SPACEDIM ? val+AMREX_SPACEDIM : val-AMREX_SPACEDIM);
    }
    //! Returns the coordinate direction.
    int coordDir () const { return val%AMREX_SPACEDIM; }
    //! Returns the orientation of the face -- low or high.
    Side faceDir () const { return Side(val/AMREX_SPACEDIM); }
    //! Returns true if Orientation is low.
    bool isLow () const { return val < AMREX_SPACEDIM; }
    //! Returns true if Orientation is high.
    bool isHigh () const { return val >= AMREX_SPACEDIM; }
    //! Read from an istream.
    friend std::istream& operator>> (std::istream& os, Orientation& o);

private:
    //! Used internally.
    Orientation (int v) : val(v) {}
    //
    // The data.
    //
    int val;
};

//! Write to an ostream in ASCII format.
std::ostream& operator<< (std::ostream& os, const Orientation& o);


//! An Iterator over the Orientation of Faces of a Box
class OrientationIter
{

public:
    //! The default constructor.
    OrientationIter ()
        :
        face(0) {}
    //! Construct an iterator on the Orientation.
    OrientationIter (const Orientation& _face)
        :
        face(_face) {}
    //! Reset (rewind) the iterator.
    void rewind () { face = 0; }
    //! Return the orientation of the face.
    Orientation operator() () const { BL_ASSERT(ok()); return Orientation(face); }
    //! Cast to void*.  Used to test if iterator is valid.
    operator void* () { return 0 <= face && face < 2*AMREX_SPACEDIM ? this : 0; }
    //! Pre-decrement.
    OrientationIter& operator-- () { BL_ASSERT(ok()); --face; return *this; }
    //! Pre-increment.
    OrientationIter& operator++ () { BL_ASSERT(ok()); ++face; return *this; }
    //! Post-decrement.
    OrientationIter operator-- (int)
    {
        BL_ASSERT(ok()); OrientationIter it(face); --face; return it;
    }
    //! Post-increment.
    OrientationIter operator++ (int)
    {
        BL_ASSERT(ok()); OrientationIter it(face); ++face; return it;
    }
    //! The equality operator.
    bool operator== (const OrientationIter& oi) const
    {
        BL_ASSERT(ok() && oi.ok()); return face == oi.face;
    }
    //! The inequality operator.
    bool operator!= (const OrientationIter& oi) const
    {
        BL_ASSERT(ok() && oi.ok()); return face != oi.face;
    }

private:

    int face;
    //! Construct an iterator on the face.
    OrientationIter (int _face) : face(_face) {}
    //! Is the iterator valid?
    bool ok () const { return 0 <= face && face < 2*AMREX_SPACEDIM; }
};

}

#endif /*BL_ORIENTATION_H*/
