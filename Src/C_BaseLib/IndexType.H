//BL_COPYRIGHT_NOTICE

#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H

//
// $Id: IndexType.H,v 1.7 2000-04-24 17:52:35 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::flush;
#else
#include <iostream.h>
#endif

#include <Misc.H>
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
    friend ostream& operator<< (ostream&         os,
                                const IndexType& itype);
    //
    //@ManDoc: Read an IndexType from an istream.
    //
    friend istream& operator>> (istream&   is,
                                IndexType& itype);
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

//
// Inlines.
//

inline
int
IndexType::mask (int k)
{
    return 1<<k;
}

inline
IndexType::IndexType ()
    : itype(0)
{}

inline
IndexType::IndexType (const IndexType& bt)
    : itype(bt.itype)
{}

inline
IndexType& IndexType::operator= (const IndexType& bt)
{
    itype = bt.itype;
    return *this;
}

inline
IndexType::IndexType (const IntVect& iv)
{
    itype = D_TERM((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2));
}

inline
IndexType::IndexType (D_DECL(CellIndex i, CellIndex j, CellIndex k))
{
    itype = D_TERM(i, | (j<<1), | (k<<2));
}

inline
void
IndexType::set (int dir)
{
    itype |= mask(dir);
}

inline
void
IndexType::unset (int dir)
{
    itype &= ~mask(dir);
}

inline
bool
IndexType::test (int dir) const
{
    return (itype & mask(dir)) != 0;
}

inline
void
IndexType::setall ()
{
    itype = (1 << BL_SPACEDIM) - 1;
}

inline
void
IndexType::clear ()
{
    itype = 0;
}

inline
bool
IndexType::any () const
{
    return itype != 0;
}

inline
bool
IndexType::ok () const
{
    return itype < (1 << BL_SPACEDIM);
}

inline
void
IndexType::flip (int i)
{
    itype ^= mask(i);
}

inline
bool
IndexType::operator== (const IndexType& t) const
{
    return t.itype == itype;
}

inline
bool
IndexType::operator!= (const IndexType& t) const
{
    return t.itype != itype;
}

inline
bool
IndexType::cellCentered () const
{
    return itype == 0;
}

inline
bool
IndexType::nodeCentered () const
{
    return itype == (1<<BL_SPACEDIM)-1;
}

inline
void
IndexType::setType (int       dir,
                    CellIndex t)
{
    t == CELL ? unset(dir) : set(dir);
}

inline
IndexType::CellIndex
IndexType::ixType (int dir) const
{
    return (CellIndex) ((itype & (1<<dir)) >> dir);
}

inline
int
IndexType::operator[] (int dir) const
{
    return test(dir);
}

inline
IntVect
IndexType::ixType () const
{
    return IntVect(D_DECL(itype&1, (itype>>1)&1, (itype>>2)&1));
}

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_INDEXTYPE_H*/
