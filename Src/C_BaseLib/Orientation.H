#ifndef BL_ORIENTATION_H
#define BL_ORIENTATION_H

//
// $Id: Orientation.H,v 1.13 2001-07-19 20:02:46 lijewski Exp $
//

#include <iostream>

#include <BLassert.H>
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
               unique integer in the range [0 .. 2*BL\_SPACEDIM-1]
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
    //@ManDoc: Read from an istream.
    //
    friend std::istream& operator>> (std::istream& os, Orientation& o);

protected:
    //
    // Used internally.
    //
    Orientation (int val);

private:
    int val;
};

//
//@ManDoc: Write to an ostream in ASCII format.
//
std::ostream& operator<< (std::ostream& os, const Orientation& o);

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

#endif /*BL_ORIENTATION_H*/
