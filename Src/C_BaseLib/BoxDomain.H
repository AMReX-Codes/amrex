#ifndef BL_BOXDOMAIN_H
#define BL_BOXDOMAIN_H
//
// $Id: BoxDomain.H,v 1.27 2002-03-26 20:04:06 lijewski Exp $
//
#include <iosfwd>

#include <IndexType.H>
#include <Box.H>
#include <IntVect.H>
#include <BoxList.H>
//
// Forward declaration.
//
class BoxDomain;

namespace BoxLib
{
    //
    //@ManDoc: Returns the complement of BoxDomain bl in Box b.
    //
    BoxDomain complementIn (const Box& b, const BoxDomain& bl);

    /*@ManDoc: Coarsen all Boxes in the domain by the refinement ratio.
               The result is placed into a new BoxDomain.
    */
    void coarsen (BoxDomain& dest, const BoxDomain& fin, int ratio);

    /*@ManDoc: Refine all Boxes in the domain by the refinement ratio
               and return the result in dest.
    */
    void refine (BoxDomain&       dest,
                 const BoxDomain& fin,
                 int              ratio);

    /*@ManDoc: Compute the intersection of BoxDomain fin with Box b
               and place the result into BoxDomain dest. 
    */
    void intersect (BoxDomain&       dest,
                    const BoxDomain& fin,
                    const Box&       b);

    /*@ManDoc: Grow each Box in BoxDomain fin by size sz and place
               the result into BoxDomain dest.
    */
    void accrete (BoxDomain&       dest,
                  const BoxDomain& fin,
                  int              sz = 1);
}

//
//@ManDoc: Output a BoxDomain to an ostream is ASCII format.
//
std::ostream& operator<< (std::ostream& os, const BoxDomain& bd);

//
//@Man:
//@Memo: A List of Disjoint Boxes
/*@Doc:

  A BoxDomain is a BoxList with the restriction that Boxes in the list
  are disjoint.

  Note that a BoxDomain is NOT a BoxList due to the protected inheritance.

  This is a concrete class, not a polymorphic one.
*/

class BoxDomain
    :
    protected BoxList
{
public:

    typedef BoxList::const_iterator const_iterator;
    //
    //@ManDoc: Construct an empty BoxDomain of IndexType::TheCellType().
    //
    BoxDomain ();
    //
    /// A Box Domain out of single Box
    //
    BoxDomain (const Box& bx);
    //
    //@ManDoc: Construct an empty BoxDomain of IndexType itype.
    //
    explicit BoxDomain (IndexType itype);
    //
    //@ManDoc: Add a Box to the domain.
    //
    void add (const Box& b);
    //
    //@ManDoc: Add all Boxes in the BoxList to the domain.
    //
    void add (const BoxList& bl);
    //
    //@ManDoc: Remove a box from the domain.
    //
    BoxDomain& rmBox (const Box& b);

    /*@ManDoc: Returns True if this BoxDomain is valid.  Checks for
               disjointness. Returns true if the BoxDomain is empty.
    */
    bool ok () const;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::clear;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::size;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::isEmpty;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::isNotEmpty;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::minimalBox;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::ixType;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::minimize;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::simplify;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::begin;
    //
    //@ManDoc: Promoted from BoxList.
    //
    using BoxList::end;
    //
    //@ManDoc: Promoted from BoxList
    //
    using BoxList::contains;
    //
    //@ManDoc: Promoted from BoxList
    //
    using BoxList::shift;
    //
    //@ManDoc: Promoted from BoxList
    //
    using BoxList::shiftHalf;
    //
    //@ManDoc: Create a BoxList from this BoxDomain.
    //
    BoxList boxList () const;
    //
    //@ManDoc: Are the BoxDomains equal?
    //
    bool operator== (const BoxDomain& rhs) const;
    //
    //@ManDoc: Are the BoxDomains not equal?
    //
    bool operator!= (const BoxDomain& rhs) const;
    //
    //@ManDoc: Refine all Boxes in the domain by the refinement ratio.
    //
    BoxDomain& refine (int ratio);
    //
    //@ManDoc: Coarsen all Boxes in the domain by the refinement ratio.
    //
    BoxDomain& coarsen (int ratio);
    //
    //@ManDoc: Intersect this BoxDomain with Box b.
    //
    BoxDomain& intersect (const Box& b);
    //
    //@ManDoc: Grow each Box by size sz.
    //
    BoxDomain& accrete (int sz);
    //
    //@ManDoc: Creates the complement of BoxDomain bl in Box b.
    //
    BoxDomain& complementIn (const Box&       b,
                             const BoxDomain& bl);
};

#endif /*BL_BOXDOMAIN_H*/
