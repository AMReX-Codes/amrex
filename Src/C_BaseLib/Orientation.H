//BL_COPYRIGHT_NOTICE

#ifndef _ORIENTATION_H_
#define _ORIENTATION_H_

//
// $Id: Orientation.H,v 1.1 1997-09-12 18:00:13 lijewski Exp $
//

#include <iostream.h>
#include <Assert.H>
#include <Boolean.H>
#include <SPACE.H>

class OrientationIter;

//
//@Man:
//@Memo: Encapsulation of the Orientation of the Faces of a Box
/*@Doc:

  This class encapsulates the orientation of the faces of a Box by
  providing an ordering of each of the faces of a Box in BL\_SPACEDIM
  dimensions.  This allows iterating over all the faces of a Box.  The
  ordering first traverses the BL\_SPACEDIM low sides from direction 0 ..
  BL\_SPACEDIM-1 and then the BL\_SPACEDIM high sides from direction 0 ..
  BL\_SPACEDIM-1.
*/

class Orientation
{
public:
    friend class OrientationIter;
    //
    //@ManDoc: In each dimension a face is either low or high.
    //
    enum Side { low = 0, high = 1 };
    //
    //@ManDoc: The default constructor.
    //
    Orientation ();
    //
    //@ManDoc: Set the orientation of a side.
    //
    Orientation (int  dir,
                 Side side);
    //
    //@ManDoc: The copy constructor.
    //
    Orientation (const Orientation& rhs);
    //
    //@ManDoc: The assignment operator.
    //
    Orientation& operator= (const Orientation& rhs);
    //
    //@ManDoc: Logical equality.
    //
    bool operator== (const Orientation& o) const;
    //
    //@ManDoc: Logical inequality.
    //
    bool operator!= (const Orientation& o) const;
    //
    //@ManDoc: Less-than.
    //
    bool operator<  (const Orientation& o) const;
    //
    //@ManDoc: Less-than or equal.
    //
    bool operator<= (const Orientation& o) const;
    //
    //@ManDoc: Greater-than.
    //
    bool operator>  (const Orientation& o) const;
    //
    //@ManDoc: Greater-than or equal.
    //
    bool operator>= (const Orientation& o) const;

    /*@ManDoc: This conversion operator maps an orientation into a
               unique integer in the range [0 .. 2*SpaceDim-1]
               according to the above ordering.
    */
    operator int () const;
    //
    //@ManDoc: Return opposite orientation.
    //
    Orientation flip () const;
    //
    //@ManDoc: Returns the coordinate direction.
    //
    int coordDir () const;
    //
    //@ManDoc: Returns the orientation of the face -- low or high.
    //
    Side faceDir () const;
    //
    //@ManDoc: Returns true if Orientation is low.
    //
    bool isLow () const;
    //
    //@ManDoc: Returns true if Orientation is high.
    //
    bool isHigh () const;
    //
    //@ManDoc: Write to an ostream in ASCII format.
    //
    friend ostream& operator<< (ostream&           os,
                                const Orientation& o);
    //
    //@ManDoc: Read from an istream.
    //
    friend istream& operator>> (istream&     os,
                                Orientation& o);
protected:
    //
    // Used internally.
    //
    Orientation (int val);

private:
    int val;
};

//
//@Man:
//@Memo: An Iterator over the Orientation of Faces of a Box
//@Doc:
//

class OrientationIter
{

public:
    //
    //@ManDoc: The default constructor.
    //
    OrientationIter ();
    //
    //@ManDoc: Construct an iterator on the Orientation.
    //
    OrientationIter (const Orientation& _face);
    //
    //@ManDoc: The copy constructor.
    //
    OrientationIter (const OrientationIter& it);
    //
    //@ManDoc: The assignment operator.
    //
    OrientationIter& operator= (const OrientationIter& it);
    //
    //@ManDoc: Reset (rewind) the iterator.
    //
    void rewind ();
    //
    //@ManDoc: Return the orientation of the face.
    //
    Orientation operator() () const;
    //
    //@ManDoc: Cast to void*.  Used to test if iterator is valid.
    //
    operator void* ();
    //
    //@ManDoc: Pre-decrement.
    //
    OrientationIter& operator-- ();
    //
    //@ManDoc: Pre-increment.
    //
    OrientationIter& operator++ ();
    //
    //@ManDoc: Post-decrement.
    //
    OrientationIter operator-- (int);
    //
    //@ManDoc: Post-increment.
    //
    OrientationIter operator++ (int);
    //
    //@ManDoc: The equality operator.
    //
    bool operator== (const OrientationIter& oi) const;
    //
    //@ManDoc: The inequality operator.
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

//
// Inlines.
//

inline
Orientation::Orientation (int _val)
    : val(_val)
{}

inline
Orientation::Orientation ()
    : val(-1)
{}

inline
Orientation::Orientation (int  _dir,
                          Side _side)
    : val(SpaceDim*_side + _dir)
{
    assert(0 <= _dir && _dir < SpaceDim);
}

inline
Orientation::Orientation (const Orientation& o)
    : val(o.val)
{}

inline
Orientation&
Orientation::operator= (const Orientation& o)
{
    val = o.val;
    return *this;
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
Orientation
Orientation::flip () const
{
    return Orientation(val < SpaceDim ? val+SpaceDim : val-SpaceDim);
}

inline
int
Orientation::coordDir () const
{
    return val%SpaceDim;
}

inline
Orientation::Side
Orientation::faceDir () const
{
    return Side(val/SpaceDim);
}

inline
bool
Orientation::isLow () const
{
    return val < SpaceDim;
}

inline
bool
Orientation::isHigh () const
{
    return val >= SpaceDim;
}

inline
OrientationIter::OrientationIter (int _face)
    : face(_face)
{}

inline
OrientationIter::OrientationIter ()
    : face(0)
{}

inline
OrientationIter::OrientationIter (const Orientation& _face)
    : face(_face)
{}

inline
bool
OrientationIter::ok () const
{
    return 0 <= face && face < 2*SpaceDim;
}

inline
OrientationIter::OrientationIter (const OrientationIter& it)
{
    assert(it.ok());
    face = it.face;
}

inline
OrientationIter&
OrientationIter::operator= (const OrientationIter& it)
{
    assert(it.ok());
    face = it.face;
    return *this;
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
    assert(ok());
    return Orientation(face);
}

inline
OrientationIter::operator void* ()
{
    return 0 <= face && face < 2*SpaceDim ? this : 0;
}

inline
OrientationIter&
OrientationIter::operator-- ()
{
    assert(ok());
    --face;
    return *this;
}

inline
OrientationIter&
OrientationIter::operator++ ()
{
    assert(ok());
    ++face;
    return *this;
}

inline
OrientationIter
OrientationIter::operator-- (int)
{
    assert(ok());
    return OrientationIter(face--);
}

inline
OrientationIter
OrientationIter::operator++ (int)
{
    assert(ok());
    return OrientationIter(face++);
}

inline
bool
OrientationIter::operator== (const OrientationIter& oi) const
{
    assert(ok() && oi.ok());
    return face == oi.face;
}

inline
bool
OrientationIter::operator!= (const OrientationIter& oi) const
{
    assert(ok() && oi.ok());
    return face != oi.face;
}

#endif /*_ORIENTATION_H_*/
