#ifndef BL_BOXDOMAIN_H
#define BL_BOXDOMAIN_H

//
// $Id: BoxDomain.H,v 1.17 2001-07-23 21:15:16 car Exp $
//

#include <cstdlib>
#include <iostream>

#include <IndexType.H>
#include <Box.H>
#include <IntVect.H>
#include <BoxList.H>

//
// Forward declaration of our BoxDomain iterator class.
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
    typedef BoxList::Iterator Iterator;
    typedef BoxList::ConstIterator ConstIterator;
    //
    //@ManDoc: Construct an empty BoxDomain of IndexType::TheCellType().
    //
    BoxDomain ();
    //
    //@ManDoc: Construct an empty BoxDomain of IndexType itype.
    //
    explicit BoxDomain (IndexType itype);
    //
    //@ManDoc: The copy constructor.
    //
    BoxDomain (const BoxDomain& rhs);
    //
    //@ManDoc: The assignment operator.
    //
    BoxDomain& operator= (const BoxDomain& rhs);
    //
    //@ManDoc: The destructor.
    //
    ~BoxDomain ();
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
    using BoxList::length;
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
    //@ManDoc: True if IntVect is within any of the Boxes in this BoxDomain.
    //
    bool contains (const IntVect& v) const;

    /*@ManDoc: True if the Box is within any of the Boxes in this BoxDomain.
               The Box must also have the same IndexType as those in this
               BoxDomain.
    */
    bool contains (const Box& b) const;
    //
    //@ManDoc: True if all Boxes in bl are contained in this BoxDomain.
    //
    bool contains (const BoxList& bl) const;
    //
    //@ManDoc: Applies Box::shift(int,int) to each Box in the BoxDomain.
    //
    BoxDomain& shift (int dir,
                      int nzones);
    //
    //@ManDoc: Applies Box::shiftHalf(int,int) to each Box in the BoxDomain.
    //
    BoxDomain& shiftHalf (int dir,
                          int num_halfs);

    //
    //@ManDoc: Applies Box::shiftHalf(IntVect) to each Box in BoxDomain.
    //
    BoxDomain& shiftHalf (const IntVect& iv);
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
