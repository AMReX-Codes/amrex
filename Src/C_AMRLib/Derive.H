//BL_COPYRIGHT_NOTICE

#ifndef _Derive_H_
#define _Derive_H_

//
// $Id: Derive.H,v 1.2 1997-11-24 18:52:29 lijewski Exp $
//

#include <REAL.H>
#include <List.H>
#include <aString.H>
#include <Box.H>
#include <Interpolater.H>
#include <StateDescriptor.H>

// ------------------------------------------------------------------
// ------------------------------------------------------------------
// --------------------  DERIVED VARIABLES --------------------------

typedef BOX (*DeriveBoxMap)(const BOX&);

/*  Function called by DeriveRec to compute derived quantity.
*/
typedef void (*DeriveFunc)(REAL* data, ARLIM_P(dlo), ARLIM_P(dhi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(compdat_lo), ARLIM_P(compdat_hi),
			   const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec,
                           const int* level, const int* grid_no) ;

/// Computes quantities derived from state data
class DeriveRec
{
  /* DeriveRec is designed to compute quantities which can be derived
    from the state data contained in AmrLevel and its derivatives. Some
    examples might be kinetic energy, vorticity, concentration gradients...
    */
   friend class DeriveList;
   
   protected:

/*    a StateRange is an element of a linked list to point to 
      state quantities in AmrLevels
*/
      struct StateRange {
///
        int               typ;
///
        int               sc;
///
        int               nc;
///
        StateRange       *next;
      };

//    name of derived quantity
      aString  derive_name;

//    type of derived quantity
      IndexType            der_type;

//    number of components in derived quantity
      int                n_derive;

//    function that computes derived quantity from the following state variables
      DeriveFunc        func;

//    interpolater for mapping crse grid derived data to finer levels
      Interpolater      *mapper;

//    type of constituent quantities
      IndexType            rng_type;

//    box mapper that specifies constituent region given derived region
      DeriveBoxMap       bx_map;

//    total number of state variables
      int                n_state;

//    number of state ranges
      int                nsr;

//    list of state data subranges required to derive given quantity
      StateRange        *rng;

//    array of bndry types
      int               *bcr;

   public:

  /// Access functions.
  ///
      const aString& name() const { return derive_name; }
  ///
      IndexType deriveType() const { return der_type; }
  ///
      DeriveFunc derFunc() const { return func; }
  ///
      DeriveBoxMap boxMap() const { return bx_map; }
  ///
      Interpolater* interp() const { return mapper; }
  ///
      IndexType rangeType() const { return rng_type; }
  ///
      int numDerive() const { return n_derive; }
  ///
      int numRange() const { return nsr; }
  ///
      int numState() const { return n_state; }
  ///
      const int* getBC() const { return bcr; }
  ///
      void getRange(int k, int& state_indx, int& src_comp,
                    int& num_comp) const;

   protected:

//    invalid constructor
      DeriveRec();

//    constructor
      DeriveRec(const aString &name, IndexType result_type,
                int nvar_derive, DeriveFunc der_func,
		DeriveBoxMap box_map, IndexType component_type,
		Interpolater *interp = &pc_interp);

//   destructor
     ~DeriveRec();

///
      void addRange(const DescriptorList& d_list, int state_indx,
                    int src_comp, int num_comp);
///
      void buildBC(const DescriptorList& d_list);
};

/// Container class for DeriveRec
class DeriveList
{
/*  DeriveList manages and provides access to the list of DeriveRecs
*/

public:
///  
    DeriveList();
///  
    ~DeriveList();
///  
    void clear();

//  Determines whether quantity identified by <name> is in the registry
    int canDerive(const aString&) const;

//  Accesses to particular record in registry
    const DeriveRec* get(const aString &name) const;

//  Adds another entry to the registry.
    void add(const aString &name, IndexType result_type,
	     int nvar_derive, DeriveFunc der_func,
	     DeriveBoxMap box_map, IndexType component_type,
	     Interpolater *interp = &pc_interp);

//  Adds another Staterange to the DeriveRec identified by <name>
    void addComponent(const aString &name, const DescriptorList& d_list,
		      int state_indx, int start_comp, int ncomp);

    List<DeriveRec*> &dlist() { return lst; }

protected:
///
    List<DeriveRec*> lst;

private:
    DeriveList(const DeriveList&);
    DeriveList& operator = (const DeriveList&);

};


#endif
