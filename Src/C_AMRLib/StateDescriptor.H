//BL_COPYRIGHT_NOTICE

#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 

//
// $Id: StateDescriptor.H,v 1.3 1997-11-26 17:46:22 lijewski Exp $
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
typedef void (*BndryFunc)(REAL* data, ARLIM_P(lo), ARLIM_P(hi),
                          const int* dom_lo, const int* dom_hi,
			  const REAL* dx, const REAL* grd_lo,
			  const REAL* time, const int* bc);
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
                       BndryFunc      func);

    void resetComponentBCs (int          comp,
                            const BCRec& bcr,
                            BndryFunc    func);

    void dumpNames (ostream& os,
                    int      start_comp,
                    int      num_comp) const;

    IndexType getType() const { return type; }

    StateDescriptor::TimeCenter timeType() const { return t_type; }

    int nComp () const { return ncomp; }

    int nExtra () const { return ngrow; }

    Interpolater* interp () const { return mapper; }

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
                       BndryFunc      func);

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
