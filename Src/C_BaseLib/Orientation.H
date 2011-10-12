
#ifndef BL_ORIENTATION_H
#define BL_ORIENTATION_H

#include <iosfwd>

#include <BLassert.H>
#include <SPACE.H>

class OrientationIter;

//
// Encapsulation of the Orientation of the Faces of a Box
//
// This class encapsulates the orientation of the faces of a Box by
// providing an ordering of each of the faces of a Box in BL_SPACEDIM
// dimensions.  This allows iterating over all the faces of a Box.  The
// ordering first traverses the BL_SPACEDIM low sides from direction 0 ..
// BL_SPACEDIM-1 and then the BL_SPACEDIM high sides from direction 0 ..
// BL_SPACEDIM-1.
//

class Orientation
{
public:

    friend class OrientationIter;
    //
    // In each dimension a face is either low or high.
    //
    enum Side { low = 0, high = 1 };
    //
    // The default constructor.
    //
    Orientation ();
    //
    // Set the orientation of a side.
    //
    Orientation (int  dir,
                 Side side);
    //
    // Logical equality.
    //
    bool operator== (const Orientation& o) const;
    //
    // Logical inequality.
    //
    bool operator!= (const Orientation& o) const;
    //
    // Less-than.
    //
    bool operator<  (const Orientation& o) const;
    //
    // Less-than or equal.
    //
    bool operator<= (const Orientation& o) const;
    //
    // Greater-than.
    //
    bool operator>  (const Orientation& o) const;
    //
    // Greater-than or equal.
    //
    bool operator>= (const Orientation& o) const;
    //
    // This conversion operator maps an orientation into a
    // unique integer in the range [0 .. 2*BL_SPACEDIM-1]
    // according to the above ordering.
    //
    operator int () const;
    //
    // Return opposite orientation.
    //
    const Orientation flip () const;
    //
    // Returns the coordinate direction.
    //
    int coordDir () const;
    //
    // Returns the orientation of the face -- low or high.
    //
    Side faceDir () const;
    //
    // Returns true if Orientation is low.
    //
    bool isLow () const;
    //
    // Returns true if Orientation is high.
    //
    bool isHigh () const;
    //
    // Read from an istream.
    //
    friend std::istream& operator>> (std::istream& os, Orientation& o);

private:
    //
    // Used internally.
    //
    Orientation (int val);
    //
    // The data.
    //
    int val;
};

//
// Write to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream& os, const Orientation& o);

//
// An Iterator over the Orientation of Faces of a Box
//

class OrientationIter
{

public:
    //
    // The default constructor.
    //
    OrientationIter ();
    //
    // Construct an iterator on the Orientation.
    //
    OrientationIter (const Orientation& _face);
    //
    // Reset (rewind) the iterator.
    //
    void rewind ();
    //
    // Return the orientation of the face.
    //
    Orientation operator() () const;
    //
    // Cast to void*.  Used to test if iterator is valid.
    //
    operator void* ();
    //
    // Pre-decrement.
    //
    OrientationIter& operator-- ();
    //
    // Pre-increment.
    //
    OrientationIter& operator++ ();
    //
    // Post-decrement.
    //
    OrientationIter operator-- (int);
    //
    // Post-increment.
    //
    OrientationIter operator++ (int);
    //
    // The equality operator.
    //
    bool operator== (const OrientationIter& oi) const;
    //
    // The inequality operator.
    //
    bool operator!= (const OrientationIter& oi) const;

protected:

    int face;
    //
    // Construct an iterator on the face.
    //
    OrientationIter (int _face);
    //
    // Is the iterator valid?
    //
    bool ok () const;
};

inline
Orientation::Orientation (int _val)
    :
    val(_val)
{}

inline
Orientation::Orientation ()
    :
    val(-1)
{}

inline
Orientation::Orientation (int  _dir,
                          Side _side)
    :
    val(BL_SPACEDIM*_side + _dir)
{
    BL_ASSERT(0 <= _dir && _dir < BL_SPACEDIM);
}

inline
bool
Orientation::operator== (const Orientation& o) const
{
    return val == o.val;
}

inline
bool
Orientation::operator!= (const Orientation& o) const
{
    return val != o.val;
}

inline
bool
Orientation::operator<  (const Orientation& o) const
{
    return val < o.val;
}

inline
bool
Orientation::operator<= (const Orientation& o) const
{
    return val <= o.val;
}

inline
bool
Orientation::operator>  (const Orientation& o) const
{
    return val > o.val;
}

inline
bool
Orientation::operator>= (const Orientation& o) const
{
    return val >= o.val;
}

inline
Orientation::operator int () const
{
    return val;
}

inline
int
Orientation::coordDir () const
{
    return val%BL_SPACEDIM;
}

inline
Orientation::Side
Orientation::faceDir () const
{
    return Side(val/BL_SPACEDIM);
}

inline
bool
Orientation::isLow () const
{
    return val < BL_SPACEDIM;
}

inline
bool
Orientation::isHigh () const
{
    return val >= BL_SPACEDIM;
}

inline
const Orientation
Orientation::flip () const
{
    return Orientation(val < BL_SPACEDIM ? val+BL_SPACEDIM : val-BL_SPACEDIM);
}

inline
OrientationIter::OrientationIter (int _face)
    :
    face(_face)
{}

inline
OrientationIter::OrientationIter ()
    :
    face(0)
{}

inline
OrientationIter::OrientationIter (const Orientation& _face)
    :
    face(_face)
{}

inline
bool
OrientationIter::ok () const
{
    return 0 <= face && face < 2*BL_SPACEDIM;
}

inline
void
OrientationIter::rewind ()
{
    face = 0;
}

inline
Orientation
OrientationIter::operator() () const
{
    BL_ASSERT(ok());
    return Orientation(face);
}

inline
OrientationIter::operator void* ()
{
    return 0 <= face && face < 2*BL_SPACEDIM ? this : 0;
}

inline
OrientationIter&
OrientationIter::operator-- ()
{
    BL_ASSERT(ok());
    --face;
    return *this;
}

inline
OrientationIter&
OrientationIter::operator++ ()
{
    BL_ASSERT(ok());
    ++face;
    return *this;
}

inline
OrientationIter
OrientationIter::operator-- (int)
{
    BL_ASSERT(ok());
    OrientationIter it(face);
    --face;
    return it;
}

inline
OrientationIter
OrientationIter::operator++ (int)
{
    BL_ASSERT(ok());
    OrientationIter it(face);
    ++face;
    return it;
}

inline
bool
OrientationIter::operator== (const OrientationIter& oi) const
{
    BL_ASSERT(ok() && oi.ok());
    return face == oi.face;
}

inline
bool
OrientationIter::operator!= (const OrientationIter& oi) const
{
    BL_ASSERT(ok() && oi.ok());
    return face != oi.face;
}

#endif /*BL_ORIENTATION_H*/
