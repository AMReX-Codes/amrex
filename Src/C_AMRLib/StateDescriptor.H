//BL_COPYRIGHT_NOTICE

#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 

//
// $Id: StateDescriptor.H,v 1.13 1999-05-10 18:54:08 car Exp $
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
//@Man:
//@Memo: Attributes of StateData.
/*@Doc:

  Attributes of StateData.
*/

class StateDescriptor
{
    friend class DescriptorList;

public:
    //
    //@ManDoc: Defines the temporal centeredness of data.
    //
    enum TimeCenter { Point = 0 , Interval };
    //
    //@ManDoc: Type of function called by BCRec for user-supplied boundary data.
    //
    typedef void (*BndryFuncDefault)(Real* data, ARLIM_P(lo), ARLIM_P(hi),
                                     const int* dom_lo, const int* dom_hi,
                                     const Real* dx, const Real* grd_lo,
                                     const Real* time, const int* bc);
    //
    // Class wrapper around BndryFuncDefault.
    //
    class BndryFunc
    {
    public:
        //
        // Bogus constructor.
        //
        BndryFunc () : m_func(0) {}
        //
        // Constructor.
        //
        BndryFunc (BndryFuncDefault inFunc) : m_func(inFunc) {}
        //
        // Return a ptr to a clone of this object.
        // It is the responsibility of the caller to delete the result.
        //
        virtual BndryFunc* clone () const
        {
            return new BndryFunc(*this);
        }
        //
        // Destructor.
        //
        virtual ~BndryFunc() {}
        //
        // Fill boundary cells.
        //
        virtual void operator () (REAL* data, const int* lo, const int* hi,
                                  const int* dom_lo, const int* dom_hi,
                                  const REAL* dx, const REAL* grd_lo,
                                  const REAL* time, const int* bc) const
        {
            BL_ASSERT(m_func != 0);
            m_func(data,ARLIM(lo),ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
        }
        
    private:

        BndryFuncDefault m_func;
    };
    //
    //@ManDoc: The default constructor.
    //
    StateDescriptor ();
    //
    //@ManDoc: Constructor that sets all data members.
    //
    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp);
    //
    //@ManDoc: The destructor.
    //
    ~StateDescriptor ();
    //
    //@ManDoc: Define the data members if constructed with default construtor.
    //
    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 Interpolater* interp);
    //
    //@ManDoc: Sets details of `comp' component.
    //
    void setComponent (int              comp,
                       const aString&   nm,
                       const BCRec&     bc,
                       const BndryFunc& func,
                       Interpolater*    interp = 0,
                       int              max_map_start_comp = -1,
                       int              min_map_end_comp   = -1);
    //
    //@ManDoc: Resets boundary condition data for `comp' component.
    //
    void resetComponentBCs (int              comp,
                            const BCRec&     bcr,
                            const BndryFunc& func)
    {
        BL_ASSERT(comp >= 0 && comp < ncomp);
        bc_func.clear(comp);
        bc_func.set(comp,func.clone());
        bc[comp] = bcr;
    }
    //
    //@ManDoc: Set interpolaters for a subset of the state vector components.
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
    //
    //@ManDoc: Cleanup interpolaters for a subset of the state vector components.
    //
    void cleanUpMaps (Interpolater**& maps, 
                      int*&           map_start_comp,
                      int*&           map_num_comp,
                      int*&           max_start_comp,
                      int*&           min_end_comp) const;
    //
    //@ManDoc: Output names of components.
    //
    void dumpNames (ostream& os,
                    int      start_comp,
                    int      num_comp) const;
    //
    //@ManDoc: Returns the IndexType.
    //
    IndexType getType () const { return type; }
    //
    //@ManDoc: Returns StateDescriptor::TimeCenter.
    //
    StateDescriptor::TimeCenter timeType () const { return t_type; }
    //
    //@ManDoc: Returns number of components.
    //
    int nComp () const { return ncomp; }
    //
    //@ManDoc: Returns the grow factor.
    //
    int nExtra () const { return ngrow; }
    //
    //@ManDoc: Returns the interpolater.
    //
    Interpolater* interp () const { return mapper; }
    //
    //@ManDoc: Returns the interpolater of specified component.
    //
    Interpolater* interp (int i) const { return mapper_comp[i]; }
    //
    //@ManDoc: Returns the name of specified component.
    //
    const aString& name (int i) const { return names[i]; }
    //
    //@ManDoc: Returns the BCRec of specified component.
    //
    const BCRec& getBC (int i) const { return bc[i]; }
    //
    //@ManDoc: Returns all BCRecs.
    //
    const Array<BCRec>& getBCs () const { return bc; }
    //
    //@ManDoc: Returns the BndryFunc of specified component.
    //
    const BndryFunc& bndryFill (int i) const { return bc_func[i]; }
    //
    //@ManDoc: Is sc>=0 \&\& sc+nc<=ncomp ?
    //
    int inRange (int sc, int nc) const { return sc>=0 && sc+nc<=ncomp; }

protected:

    IndexType         type;    // Cell centered, node centered ...
    TimeCenter        t_type;  // Temporal centering
    int               id;      // Unique id
    int               ncomp;   // Number of components
    int               ngrow;   // Grow factor
    Interpolater*     mapper;  // Default interpolator
    Array<aString>    names;   // Printable names of components
    Array<BCRec>      bc;      // Array of bndry types for entire level
    PArray<BndryFunc> bc_func; // PArray of pointers to bndry fill functions
    //
    // If mapper_comp[icomp] != 0, that map is used instead of mapper
    // when interpolating icomp; otherwise, mapper is used.
    //
    Array<Interpolater*> mapper_comp;
    //
    // For component icomp, max_map_start_comp[icomp] and
    // min_map_end_comp[icomp] represent the maximum starting index and
    // minimum ending index with which to use mapper_comp[icomp] when
    // interpolating component icomp.  These are convenient for using
    // with "coupled" multiple component maps like CellConservativeLinear.
    //
    // If there is no need to "couple" a component, max_map_start_comp[icomp]
    // and min_map_end_comp[icomp] simply equal icomp.
    //
    Array<int> max_map_start_comp;
    Array<int> min_map_end_comp;
};

//
//@Man:
//@Memo: A List of StateDescriptors.
/*@Doc:

  A container class for StateDescriptors.
*/

class DescriptorList
{
public:
    //
    //@ManDoc: The constructor.
    //
    DescriptorList ();
    //
    //@ManDoc: The destructor.
    //
    ~DescriptorList ();
    //
    //@ManDoc: Set the list to its default state.
    //
    void clear () { desc.clear(); }
    //
    //@ManDoc: Returns number of elements in the list.
    //
    int length () const { return desc.length(); }
    //
    //@ManDoc: Adss new StateDescriptor at index `indx' to list.
    //
    void addDescriptor (int                            indx,
                        IndexType                      typ,
                        StateDescriptor::TimeCenter    ttyp,
                        int                            nextra,
                        int                            num_comp,
                        Interpolater*                  interp);
    //
    //@ManDoc: Calls resetComponentBCs() on StateDescriptor at index `indx'.
    //
    void resetComponentBCs (int                               indx,
                            int                               comp,
                            const BCRec&                      bc,
                            const StateDescriptor::BndryFunc& func)
    {
        desc[indx].resetComponentBCs(comp,bc,func);
    }
    //
    //@ManDoc: Calls setComponent() on StateDescriptor at index `indx'.
    //
    void setComponent (int                               indx,
                       int                               comp,
                       const aString&                    nm,
                       const BCRec&                      bc,
                       const StateDescriptor::BndryFunc& func,
                       Interpolater*                     interp = 0,
                       int                               max_map_start_comp = -1,
                       int                               min_map_end_comp   = -1)
    {
        desc[indx].setComponent(comp,nm,bc,func,interp,max_map_start_comp,
                                min_map_end_comp);
    }
    //
    //@ManDoc: Returns StateDescriptor at index `k'.
    //
    const StateDescriptor& operator[] (int k) const { return desc[k]; }

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    PArray<StateDescriptor> desc;
};

#endif /*_StateDescriptor_H_*/
