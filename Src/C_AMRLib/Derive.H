//BL_COPYRIGHT_NOTICE

#ifndef _Derive_H_
#define _Derive_H_

//
// $Id: Derive.H,v 1.12 1998-05-11 17:06:24 lijewski Exp $
//

#include <REAL.H>
#include <List.H>
#include <aString.H>
#include <Box.H>
#include <Interpolater.H>
#include <StateDescriptor.H>

extern "C"
{
//
// Type of extern "C" function called by DeriveRec to compute derived quantity.
//
typedef void (*DeriveFunc)(REAL* data, ARLIM_P(dlo), ARLIM_P(dhi),
                           const int* nvar, const REAL* compdat,
                           ARLIM_P(compdat_lo), ARLIM_P(compdat_hi),
                           const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const REAL* delta, const REAL* xlo,
                           const REAL* time, const REAL* dt, 
                           const int* bcrec,
                           const int* level, const int* grid_no) ;
}

//
//@Man:
//@Memo: Derived Type Record
/*@Doc:

  Computes quantities derived from state data.

  DeriveRec is designed to compute quantities which can be derived
  from the state data contained in AmrLevel and its derivatives. Some
  examples might be kinetic energy, vorticity, concentration gradients ...
*/

class DeriveRec
{
   friend class DeriveList;

public:
    //
    //@ManDoc: A userful typedef -- pointer to function taking and returning a Box. 
    //
    typedef Box (*DeriveBoxMap)(const Box&);
    //
    //@ManDoc: The destructor.
    //
    ~DeriveRec ();
    //
    //@ManDoc: The name of the derived type.
    //
    const aString& name () const { return derive_name; }
    //
    //@ManDoc: The IndexType of the derived type.
    //
    IndexType deriveType () const { return der_type; }
    //
    //@ManDoc: The DeriveFunc used to calculate the derived type.
    //
    DeriveFunc derFunc () const { return func; }
    //
    //@ManDoc: Maps state data box to derived data box.
    //
    DeriveBoxMap boxMap () const { return bx_map; }
    //
    //@ManDoc: Type of interpolater to use in computing derived type.
    //
    Interpolater* interp () const { return mapper; }
    //
    //@ManDoc: Number of components in the derived type.
    //
    int numDerive () const { return n_derive; }
    //
    //@ManDoc: Number of different chunks of state data needed for derived type.
    //
    int numRange () const { return nsr; }
    //
    //@ManDoc: Total number of state variables needed for derived type.
    //
    int numState () const { return n_state; }
    //
    //@ManDoc: The boundary conditions.
    //
    const int* getBC () const { return bcr; }

    /*@ManDoc: Sets state\_indx, src\_comp and num\_comp for the kth
               range (or chunk) of state data components needed to calculate
               the derived quantity.
    */
    void getRange (int  k,
                   int& state_indx,
                   int& src_comp,
                   int& num_comp) const;
protected:
    //
    // Constructor.
    //
    DeriveRec (const aString& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc     der_func,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);

    void addRange (const DescriptorList& d_list,
                   int                   state_indx,
                   int                   src_comp,
                   int                   num_comp);

    void buildBC (const DescriptorList& d_list);
    //
    // An element of a linked list to point to state quantities in AmrLevels.
    //
    struct StateRange
    {
        int         typ;
        int         sc;
        int         nc;
        StateRange* next;
    };
    //
    // Name of derived quantity.
    //
    aString derive_name;
    //
    // Type of derived quantity.
    //
    IndexType der_type;
    //
    // Number of components in derived quantity.
    //
    int n_derive;
    //
    // Function that computes derived quantity from state variables.
    //
    DeriveFunc func;
    //
    // Interpolater for mapping crse grid derived data to finer levels.
    //
    Interpolater* mapper;
    //
    // Box mapper that specifies constituent region given derived region.
    //
    DeriveBoxMap bx_map;
    //
    // Total number of state variables.
    //
    int n_state;
    //
    // Number of state ranges.
    //
    int nsr;
    //
    // List of state data subranges required to derive given quantity.
    //
    StateRange* rng;
    //
    // Array of bndry types.
    //
    int* bcr;
};

//
//@Man:
//@Memo: A list of DeriveRecs.
/*@Doc:

  DeriveList manages and provides access to the list of DeriveRecs.
*/

class DeriveList
{
public:
    //
    //@ManDoc: The default constructor.
    //
    DeriveList () {}
    //
    //@ManDoc: The destructor.
    //
    ~DeriveList () {}
    //
    //@ManDoc: Determines whether quantity identified by <name> is in the registry.
    //
    bool canDerive (const aString& name) const;
    //
    //@ManDoc: Access the particular record in registry.
    //
    const DeriveRec* get (const aString& name) const;
    //
    //@ManDoc: Adds another entry to the registry.
    //
    void add (const aString&          name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp)
    {
        lst.add(DeriveRec(name,
                          result_type,
                          nvar_derive,
                          der_func,
                          box_map,
                          interp));
    }
    //
    //@ManDoc: Adds another StateRange to the DeriveRec identified by <name>.
    //
    void addComponent (const aString&        name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   start_comp,
                       int                   ncomp);

    List<DeriveRec>& dlist () { return lst; }

protected:

    List<DeriveRec> lst;

private:
    //
    // Disallowed.
    //
    DeriveList (const DeriveList&);
    DeriveList& operator= (const DeriveList&);
};


#endif /*_Derive_H_*/
