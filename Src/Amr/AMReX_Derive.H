#ifndef AMREX_Derive_H_
#define AMREX_Derive_H_

#include <list>
#include <string>

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Interpolater.H>

extern "C"
{
    //
    // Type of extern "C" function called by DeriveRec to compute derived quantity.
    //
    // Note that ARLIM_P will be preprocessed into DIM const int&'s.
    typedef void (*DeriveFunc)(amrex::Real* data, ARLIM_P(dlo), ARLIM_P(dhi),
                               const int* nvar, const amrex::Real* compdat,
                               ARLIM_P(compdat_lo), ARLIM_P(compdat_hi),
                               const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const amrex::Real* delta, const amrex::Real* xlo,
                               const amrex::Real* time, const amrex::Real* dt, 
                               const int* bcrec,
                               const int* level, const int* grid_no) ;
    // This is dimension agnostic.  For example, dlo always has three elements.
    typedef void (*DeriveFunc3D)(amrex::Real* data, const int* dlo, const int* dhi, const int* nvar, 
				 const amrex::Real* compdat, const int* clo, const int* chi, const int* ncomp,
				 const int* lo, const int* hi,
				 const int* domain_lo, const int* domain_hi,
				 const amrex::Real* delta, const amrex::Real* xlo,
				 const amrex::Real* time, const amrex::Real* dt, 
				 const int* bcrec,
				 const int* level, const int* grid_no) ;
}

namespace amrex {

class DescriptorList;

//
// Derived Type Record
//
// Computes quantities derived from state data.
//
// DeriveRec is designed to compute quantities which can be derived
// from the state data contained in AmrLevel and its derivatives. Some
// examples might be kinetic energy, vorticity, concentration gradients ...
//

class DeriveRec
{
   friend class DeriveList;

public:
    //
    // A pointer to function taking and returning a Box. 
    //
    typedef Box (*DeriveBoxMap)(const Box&);

    static Box TheSameBox (const Box& box);
    static Box GrowBoxByOne (const Box& box);
    //
    // The destructor.
    //
    ~DeriveRec ();
    //
    // The name of the derived type.
    //
    const std::string& name () const;
    //
    // The names of components
    //
    const std::string& variableName (int comp) const;
    //
    // The IndexType of the derived type.
    //
    IndexType deriveType () const;
    //
    // The DeriveFunc used to calculate the derived type.
    //
    DeriveFunc   derFunc   () const;
    DeriveFunc3D derFunc3D () const;
    //
    // Maps state data box to derived data box.
    //
    DeriveBoxMap boxMap () const;
    //
    // Type of interpolater to use in computing derived type.
    //
    Interpolater* interp () const;
    //
    // Number of components in the derived type.
    //
    int numDerive () const;
    //
    // Number of different chunks of state data needed for derived type.
    //
    int numRange () const;
    //
    // Total number of state variables needed for derived type.
    //
    int numState () const;
    //
    // The boundary conditions.
    //
    const int* getBC () const;
    //
    // Sets state_indx, src_comp and num_comp for the kth
    // range (or chunk) of state data components needed to calculate
    // the derived quantity.
    //
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
    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc3D   der_func_3d,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);

    //
    // Constructor without a Fortran function
    //
    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveRec::DeriveBoxMap box_map);
    //
    // Constructor.
    //
    DeriveRec (const std::string&  name,
               IndexType           result_type,
               int                 nvar_derive,
	       Vector<std::string>& var_names,
               DeriveFunc          der_func,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);
    DeriveRec (const std::string&  name,
               IndexType           result_type,
               int                 nvar_derive,
	       Vector<std::string>& var_names,
               DeriveFunc3D        der_func_3d,
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

private:
    //
    // Name of derived quantity.
    //
    std::string derive_name;
    //    
    // Names of derived variables
    //
    Vector<std::string> variable_names;  
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
    DeriveFunc   func;
    DeriveFunc3D func_3d;
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
// A list of DeriveRecs.
//
// DeriveList manages and provides access to the list of DeriveRecs.
//

class DeriveList
{
public:
    //
    // The default constructor.
    //
    DeriveList ();
    //
    // Determines whether quantity identified by <name> is in the registry.
    //
    bool canDerive (const std::string& name) const;
    //
    // Access the particular record in registry.
    //
    const DeriveRec* get (const std::string& name) const;
    //
    // Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc3D            der_func_3d,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);

    //
    // This version doesn't take a Fortran function.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveRec::DeriveBoxMap box_map = &DeriveRec::TheSameBox );
    //
    //
    // Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
	      Vector<std::string>&     var_names,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
	      Vector<std::string>&     var_names,
              DeriveFunc3D            der_func_3d,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);
    //
    // Adds another StateRange to the DeriveRec identified by <name>.
    //
    void addComponent (const std::string&    name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   start_comp,
                       int                   ncomp);
 
    std::list<DeriveRec>& dlist ();

    void clear () { lst.clear(); }

private:
    //
    // Disallowed.
    //
    DeriveList (const DeriveList&);
    DeriveList& operator= (const DeriveList&);

    std::list<DeriveRec> lst;
};

}

#endif /*_Derive_H_*/
