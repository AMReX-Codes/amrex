//BL_COPYRIGHT_NOTICE

#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 

//
// $Id: StateDescriptor.H,v 1.2 1997-11-24 18:52:32 lijewski Exp $
//

#include <Array.H>
#include <PArray.H>
#include <aString.H>
#include <REAL.H>
#include <Box.H>
#include <BCRec.H>
#include <ArrayLim.H>

class Interpolater;

/// Function called by BCRec for user-supplied boundary data
typedef void (*BndryFunc)(REAL* data, ARLIM_P(lo), ARLIM_P(hi),
                          const int* dom_lo, const int* dom_hi,
			  const REAL* dx, const REAL* grd_lo,
			  const REAL* time, const int* bc);

enum TimeCenter { Point=0 , Interval };

class StateDescriptor;

/// Attributes of state at this level
class StateDescriptor
{
/** Attributes of state at this level
*/
    friend class DescriptorList;
protected:

//  cell centered, node centered ...
    IndexType        type;

//  temporal centering
    TimeCenter     t_type;

//  unique id
    int            id;

//  number of components
    int            ncomp;

//  grow factor
    int            ngrow;

//  default interpolator
    Interpolater*  mapper;

//  printable names of components
    Array<aString> names;

//  array of bndry types for entire level
    Array<BCRec>   bc;

//  array of pointers to bndry fill functions
    Array<BndryFunc> bc_func;

public:

//  bogus constructor
    StateDescriptor();

//  constructor
    StateDescriptor(IndexType btyp, TimeCenter ttyp, int ident,
                    int nextra, int num_comp,
                    Interpolater *interp);
//  destructor
    ~StateDescriptor();

//  initialize data members
    void define(IndexType btyp, TimeCenter ttyp, int ident,
                int nextra, int num_comp,
                Interpolater *interp);

//  set (define) a state vector component
    void setComponent(int comp, const aString &nm, const BCRec &bc,
                      BndryFunc func);

//  set (define) boundary conditions for a state vector component
    void resetComponentBCs(int comp, const BCRec &bcr,
                      BndryFunc func);

public:

public:
/// access functions
///
    void dumpNames(ostream& os, int start_comp, int num_comp) const;
///
    IndexType getType() const { return type; }
///
    TimeCenter timeType() const { return t_type; }
///
    int nComp() const { return ncomp; }
///
    int nExtra() const { return ngrow; }
///
    Interpolater* interp() const { return mapper; }
///
    const aString& name(int i) const { return names[i]; }
///
    const BCRec& getBC(int i) const { return bc[i]; }
///
    const Array<BCRec>& getBCs() const { return bc; }
///
    BndryFunc bndryFill(int i) const { return bc_func[i]; }
///
    int inRange(int sc, int nc) const { return (sc>=0 && sc+nc<=ncomp); }
///
};

class DescriptorList
{
/** Container class for StateDescriptor
*/
private:
    PArray<StateDescriptor> desc;
    DescriptorList(const DescriptorList&);
    DescriptorList& operator = (const DescriptorList&);

public:

///
    DescriptorList();
///
    ~DescriptorList();
///
    void clear();
///
    int length() const { return desc.length(); }
///
    void addDescriptor(int indx, IndexType typ, TimeCenter ttyp, int nextra,
                       int num_comp, Interpolater *interp);
///
    void resetComponentBCs(int indx, int comp, const BCRec& bc,
                      BndryFunc func);
///
    void setComponent(int indx, int comp, const aString &nm, const BCRec& bc,
                      BndryFunc func);
///
    const StateDescriptor& operator[] (int k) const { return desc[k]; }
};
#endif
