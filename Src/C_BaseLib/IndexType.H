#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H

//
// $Id: IndexType.H,v 1.9 2001-07-17 23:02:22 lijewski Exp $
//

#include <iostream>

#include <IntVect.H>
#include <SPACE.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

//
//@Man:
//@Memo: Cell-Based or Node-Based Indices
/*@Doc:

  The class IndexType defines an index as being cell based or node (edge)
  based in each of the BL\_SPACEDIM directions.  This class defines an
  enumerated type CellIndex to be either CELL or NODE; i.e. each of the
  BL\_SPACEDIM dimensions must be either CELL or NODE.
*/

class IndexType
{
public:
    //
    //@ManDoc: The cell index type: one of CELL or NODE.
    //
    enum CellIndex { CELL = 0, NODE = 1 };
    //
    //@ManDoc: The default constructor
    //
    IndexType ();
    //
    //@ManDoc: The copy constructor.
    //
    IndexType (const IndexType& rhs);
    //
    //@ManDoc: Construct an IndexType identical to an IntVect.
    //
    explicit IndexType (const IntVect& iv);
    //
    //@ManDoc: The assignment operator.
    //
    IndexType& operator= (const IndexType& rhs);

    /*@ManDoc: Construct an IndexType given an explicit CellIndex for
               each direction.  D\_DECL is a macro that sets the constructor
               to take BL\_SPACEDIM arguments.
    */
    IndexType (D_DECL(CellIndex i, CellIndex j, CellIndex k));
    //
    //@ManDoc: Set IndexType to be NODE based in direction dir.
    //
    void set (int dir);
    //
    //@ManDoc: Set IndexType to be CELL based in direction dir.
    //
    void unset (int dir);
    //
    //@ManDoc: True if IndexType is NODE based in direction dir.
    //
    bool test (int dir) const;
    //
    //@ManDoc: Set NODE based in all directions.
    //
    void setall ();
    //
    //@ManDoc: Set CELL based in all directions.
    //
    void clear ();
    //
    //@ManDoc: True if this IndexType is NODE based in any direction.
    //
    bool any () const;
    //
    //@ManDoc: True if IndexType is valid.
    //
    bool ok () const;
    //
    //@ManDoc: Change from CELL to NODE or NODE to CELL in direction dir.
    //
    void flip (int i);
    //
    //@ManDoc: True if IndexTypes are identical.
    //
    bool operator== (const IndexType& t) const;
    //
    //@ManDoc: True if IndexTypes are not identical.
    //
    bool operator!= (const IndexType& t) const;
    //
    //@ManDoc: True if the IndexType is CELL based in all directions.
    //
    bool cellCentered () const;
    //
    //@ManDoc: True if the IndexType is NODE based in all directions.
    //
    bool nodeCentered () const;
    //
    //@ManDoc: Set IndexType to CellIndex type t in direction dir.
    //
    void setType (int       dir,
                  CellIndex t);
    //
    //@ManDoc: Returns  the CellIndex in direction dir.
    //
    CellIndex ixType (int dir) const;
    //
    //@ManDoc: Return an integer representing the IndexType in direction dir.
    //
    int operator[] (int dir) const;
    //
    //@ManDoc: Fill an IntVect of size BL\_SPACEDIM with IndexTypes.
    //
    IntVect ixType () const;

    /*@ManDoc: This static member function returns an IndexType object of value
               IndexType::CELL.  It is provided as a convenience to our users
               when defining a Box all of whose faces should be of type
               IndexType::CELL.
      */
    static IndexType TheCellType ();

    /*@ManDoc: This static member function returns an IndexType object of value
               IndexType::NODE.  It is provided as a convenience to our users
               when defining a Box all of whose faces should be of type
               IndexType::NODE.
      */
    static IndexType TheNodeType ();
    //
    //@ManDoc: Write an IndexType to an ostream in ASCII.
    //
    friend std::ostream& operator<< (std::ostream&    os,
                                     const IndexType& itype);
    //
    //@ManDoc: Read an IndexType from an istream.
    //
    friend std::istream& operator>> (std::istream& is,
                                     IndexType&    itype);
private:
    //
    // Returns 1<<k.
    //
    static int mask (int k);
    //
    // An integer holding the CellIndex in bits 0 - BL\_SPACEDIM-1.
    //
    unsigned int itype;
};

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_INDEXTYPE_H*/
