//BL_COPYRIGHT_NOTICE

#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 

//
// $Id: StateDescriptor.H,v 1.6 1998-01-06 23:43:58 lijewski Exp $
//

#include <Array.H>
#include <PArray.H>
#include <aString.H>
#include <REAL.H>
#include <Box.H>
#include <BCRec.H>
#include <ArrayLim.H>

class Interpolater;

//
// Function called by BCRec for user-supplied boundary data.
//
typedef void (*BndryFunc)(Real* data, ARLIM_P(lo), ARLIM_P(hi),
                          const int* dom_lo, const int* dom_hi,
                          const Real* dx, const Real* grd_lo,
                          const Real* time, const int* bc);
//
// Attributes of state at this level.
//

class StateDescriptor
{
    friend class DescriptorList;

public:

    enum TimeCenter { Point = 0 , Interval };

    StateDescriptor ();

    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp);
    ~StateDescriptor ();

    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 Interpolater* interp);

    void setComponent (int            comp,
                       const aString& nm,
                       const BCRec&   bc,
                       BndryFunc      func,
                       Interpolater*  interp=0,
                       int            max_map_start_comp = -1,
                       int            min_map_end_comp   = -1);

    void resetComponentBCs (int          comp,
                            const BCRec& bcr,
                            BndryFunc    func);
    //
    // Set and clean up interpolaters for a subset of the state vector
    // components.
    //
    void setUpMaps (int&                use_default_map,
                    const Interpolater* default_map,
                    int                 start_comp,
                    int                 num_comp,
                    Interpolater**&     maps, 
                    int&                nmaps,
                    int*&               map_start_comp, 
                    int*&               map_num_comp,
                    int*&               max_start_comp,
                    int*&               min_end_comp) const;

    void cleanUpMaps (Interpolater**& maps, 
                      int*&           map_start_comp,
                      int*&           map_num_comp,
                      int*&           max_start_comp,
                      int*&           min_end_comp) const;

    void dumpNames (ostream& os,
                    int      start_comp,
                    int      num_comp) const;

    IndexType getType() const { return type; }

    StateDescriptor::TimeCenter timeType() const { return t_type; }

    int nComp () const { return ncomp; }

    int nExtra () const { return ngrow; }

    Interpolater* interp () const { return mapper; }

    Interpolater* interp (int i) const { return mapper_comp[i]; }

    const aString& name (int i) const { return names[i]; }

    const BCRec& getBC (int i) const { return bc[i]; }

    const Array<BCRec>& getBCs () const { return bc; }

    BndryFunc bndryFill (int i) const { return bc_func[i]; }

    int inRange (int sc, int nc) const { return sc>=0 && sc+nc<=ncomp; }

protected:

    IndexType        type;    // Cell centered, node centered ...
    TimeCenter       t_type;  // Temporal centering
    int              id;      // Unique id
    int              ncomp;   // Number of components
    int              ngrow;   // Grow factor
    Interpolater*    mapper;  // Default interpolator
    Array<aString>   names;   // Printable names of components
    Array<BCRec>     bc;      // Array of bndry types for entire level
    Array<BndryFunc> bc_func; // Array of pointers to bndry fill functions

//  array of pointers to interpolaters
//  If mapper_comp[icomp] != 0, that
//  map is used instead of mapper when
//  interpolating icomp; otherwise, mapper is used.

    Array<Interpolater*> mapper_comp;

//  for component icomp, max_map_start_comp[icomp] and
//  min_map_end_comp[icomp] represent the maximum starting
//  index and minimum ending index with which to use
//  mapper_comp[icomp] when interpolating component icomp. 
//  These are convenient for using
//  with "coupled" multiple component maps like
//  CellConservativeLinear.

//  If there is no need to "couple" a component,
//  max_map_start_comp[icomp] and min_map_end_comp[icomp]
//  simply equal icomp

    Array<int> max_map_start_comp;
    Array<int> min_map_end_comp;
};

class DescriptorList
{
public:

    DescriptorList ();

    ~DescriptorList ();

    void clear ();

    int length () const { return desc.length(); }

    void addDescriptor (int                            indx,
                        IndexType                      typ,
                        StateDescriptor::TimeCenter    ttyp,
                        int                            nextra,
                        int                            num_comp,
                        Interpolater*                  interp);

    void resetComponentBCs (int          indx,
                            int          comp,
                            const BCRec& bc,
                            BndryFunc    func);

    void setComponent (int            indx,
                       int            comp,
                       const aString& nm,
                       const BCRec&   bc,
                       BndryFunc      func,
                       Interpolater*  interp = 0,
                       int            max_map_start_comp = -1,
                       int            min_map_end_comp   = -1);

    const StateDescriptor& operator[] (int k) const { return desc[k]; }

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    PArray<StateDescriptor> desc;
};
#endif
