
#ifndef _BOUNDCOND_H_
#define _BOUNDCOND_H_

//
// $Id: BoundCond.H,v 1.8 2000-10-02 20:51:15 lijewski Exp $
//

#include <BLassert.H>

//@Man:
/*@Memo:
  Maintain an identifier for boundary condition types.
*/        
/*@Doc:
        This is a placeholder for more extensive boundary condition
        implementations, which might include stencils, etc.  Presently,
        boundary conditions are specified via an integer identifier.
        This class maintains that integer.
*/

class BoundCond
{
public:
    //
    //@ManDoc: Initializes boundary condition type to default.
    //
    BoundCond () : bctype(-1) {}
    //
    //@ManDoc: Initializes boundary condition type to specified value.
    //
    BoundCond (int _bctype) : bctype(_bctype) {}
    //
    //@ManDoc: Copy constructor.
    //
    BoundCond (const BoundCond& rhs) : bctype(rhs.bctype) {}
    //
    //@ManDoc: Copy assignment operator.
    //
    BoundCond& operator= (const BoundCond& rhs)
    {
        if (!(this == &rhs))
        {
            bctype = rhs.bctype;
        }
        return *this;
    }
    //
    //@ManDoc: Return index of boundary condition type.
    //
    operator int () const
    {
        BL_ASSERT(bctype != -1);
        return bctype;
    }
    //
    //@ManDoc: Set index of boundary condition type.
    //
    BoundCond& operator= (const int _bctype)
    {
        bctype = _bctype;
        return *this;
    }
private:
    //
    // Index of boundary condition type.
    //
    int bctype;
};

#endif /*_BOUNDCOND_H_*/
