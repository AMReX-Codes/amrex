#ifndef _Derive_H_
#define _Derive_H_
//
// $Id: Derive.H,v 1.15 2001-08-01 21:50:45 lijewski Exp $
//
#include <list>
#include <string>

#include <ArrayLim.H>
#include <REAL.H>
#include <Box.H>
#include <Interpolater.H>

class DescriptorList;

extern "C"
{
    //
    // Type of extern "C" function called by DeriveRec to compute derived quantity.
    //
    typedef void (*DeriveFunc)(Real* data, ARLIM_P(dlo), ARLIM_P(dhi),
                               const int* nvar, const Real* compdat,
                               ARLIM_P(compdat_lo), ARLIM_P(compdat_hi),
                               const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* delta, const Real* xlo,
                               const Real* time, const Real* dt, 
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
    //@ManDoc: A pointer to function taking and returning a Box. 
    //
    typedef Box (*DeriveBoxMap)(const Box&);
    //
    //@ManDoc: The destructor.
    //
    ~DeriveRec ();
    //
    //@ManDoc: The name of the derived type.
    //
    const std::string& name () const;
    //
    //@ManDoc: The names of components
    //
    const std::string& variableName (int comp) const;
    //
    //@ManDoc: The IndexType of the derived type.
    //
    IndexType deriveType () const;
    //
    //@ManDoc: The DeriveFunc used to calculate the derived type.
    //
    DeriveFunc derFunc () const;
    //
    //@ManDoc: Maps state data box to derived data box.
    //
    DeriveBoxMap boxMap () const;
    //
    //@ManDoc: Type of interpolater to use in computing derived type.
    //
    Interpolater* interp () const;
    //
    //@ManDoc: Number of components in the derived type.
    //
    int numDerive () const;
    //
    //@ManDoc: Number of different chunks of state data needed for derived type.
    //
    int numRange () const;
    //
    //@ManDoc: Total number of state variables needed for derived type.
    //
    int numState () const;
    //
    //@ManDoc: The boundary conditions.
    //
    const int* getBC () const;

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
    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc     der_func,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);
    //
    // Constructor.
    //
    DeriveRec (const std::string&  name,
               IndexType           result_type,
               int                 nvar_derive,
	       Array<std::string>& var_names,
               DeriveFunc          der_func,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

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
    std::string derive_name;
    //    
    // Names of derived variables
    //
    Array<std::string> variable_names;  
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
    DeriveList ();
    //
    //@ManDoc: Determines whether quantity identified by <name> is in the registry.
    //
    bool canDerive (const std::string& name) const;
    //
    //@ManDoc: Access the particular record in registry.
    //
    const DeriveRec* get (const std::string& name) const;
    //
    //@ManDoc: Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);
    //
    //@ManDoc: Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
	      Array<std::string>&     var_names,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);
    //
    //@ManDoc: Adds another StateRange to the DeriveRec identified by <name>.
    //
    void addComponent (const std::string&    name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   start_comp,
                       int                   ncomp);
 
    std::list<DeriveRec>& dlist ();

private:
    //
    // Disallowed.
    //
    DeriveList (const DeriveList&);
    DeriveList& operator= (const DeriveList&);

    std::list<DeriveRec> lst;
};

#endif /*_Derive_H_*/
