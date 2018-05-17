
#ifndef AMREX_StateDescriptor_H_
#define AMREX_StateDescriptor_H_

#include <utility>
#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_PhysBCFunct.H>

namespace amrex {

class Interpolater;

//
// Attributes of StateData.
//

class StateDescriptor
{
    friend class DescriptorList;

public:
    //
    // Defines the temporal centeredness of data.
    //
    enum TimeCenter { Point = 0 , Interval };
    //
    // Type of function called by BCRec for user-supplied boundary data.
    //
    // Class wrapper around BndryFuncDefault.
    //
    class BndryFunc
	: public BndryFunctBase
    {
    public:
        //
        // Bogus constructor.
        //
        BndryFunc ();
        //
        // A Constructor.
        //
        BndryFunc (BndryFuncDefault inFunc);
        BndryFunc (BndryFunc3DDefault inFunc);
        //
        // Another Constructor.
        //
        BndryFunc (BndryFuncDefault inFunc,BndryFuncDefault gFunc);
        BndryFunc (BndryFunc3DDefault inFunc,BndryFunc3DDefault gFunc);
        //
        // Return a ptr to a clone of this object.
        // It is the responsibility of the caller to delete the result.
        //
        virtual BndryFunc* clone () const override;
        //
        // Destructor.
        //
        virtual ~BndryFunc () override;
	//
	// Fill boundary cells using "regular" function.
	//
	virtual void operator () (Real* data, const int* lo, const int* hi,
				  const int* dom_lo, const int* dom_hi,
				  const Real* dx, const Real* grd_lo,
				  const Real* time, const int* bc) const override;
        //
        // Fill boundary cells using "group" function.
        //
        virtual void operator () (Real* data, const int* lo, const int* hi,
                                  const int* dom_lo, const int* dom_hi,
                                  const Real* dx, const Real* grd_lo,
                                  const Real* time, const int* bc, int ng) const;

	void Print() const;

    private:

        BndryFuncDefault   m_gfunc;
        BndryFunc3DDefault m_gfunc3D;
    };
    //
    // The default constructor.
    //
    StateDescriptor ();
    //
    // Constructor that sets all data members.
    //
    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp,
                     bool          extrap = false,
                     bool          store_in_checkpoint = true);
    //
    // The destructor.
    //
    ~StateDescriptor ();
    //
    // Define the data members if constructed with default construtor.
    //
    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 Interpolater* interp,
                 bool          extrap = false,
                 bool          store_in_checkpoint = true);
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp = 0,
                       int                max_map_start_comp = -1,
                       int                min_map_end_comp   = -1);
    //
    // Sets details of comp component.
    //
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp,
                       bool               master_or_slave,
                       int                groupsize);
    //
    // Resets boundary condition data for comp component.
    //
    void resetComponentBCs (int              comp,
                            const BCRec&     bcr,
                            const BndryFunc& func);
    //
    // Set interpolaters for a subset of the state vector components.
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
    // Cleanup interpolaters for a subset of the state vector components.
    //
    void cleanUpMaps (Interpolater**& maps, 
                      int*&           map_start_comp,
                      int*&           map_num_comp,
                      int*&           max_start_comp,
                      int*&           min_end_comp) const;
    //
    // Output names of components.
    //
    void dumpNames (std::ostream& os,
                    int           start_comp,
                    int           num_comp) const;
    //
    // Returns the IndexType.
    //
    IndexType getType () const;
    //
    // Returns StateDescriptor::TimeCenter.
    //
    StateDescriptor::TimeCenter timeType () const;
    //
    // Returns number of components.
    //
    int nComp () const;
    //
    // Returns the grow factor.
    //
    int nExtra () const;
    //
    // Returns the interpolater.
    //
    Interpolater* interp () const;
    //
    // Returns the interpolater of specified component.
    //
    Interpolater* interp (int i) const;
    //
    // Returns the name of specified component.
    //
    const std::string& name (int i) const;
    //
    // Returns the BCRec of specified component.
    //
    const BCRec& getBC (int i) const;
    //
    // Returns all BCRecs.
    //
    const Vector<BCRec>& getBCs () const;
    //
    // Returns the BndryFunc of specified component.
    //
    const BndryFunc& bndryFill (int i) const;
    //
    // Is sc>=0 \&\& sc+nc<=ncomp ?
    //
    int inRange (int sc, int nc) const;
    //
    // Are the interpolaters in the specified range identical?
    //
    bool identicalInterps (int scomp, int ncomp) const;
    //
    // Returns contiguous ranges of comps with identical interpolaters.
    //
    std::vector< std::pair<int,int> > sameInterps (int scomp, int ncomp) const;
    //
    // Can extrapolate in time.
    //
    bool extrap () const;
    //
    // Should store this StateData in the checkpoint file
    //
    bool store_in_checkpoint () const;

    bool master (int i) const { return m_master[i]; }

    int groupsize (int i) const { return m_groupsize[i]; }

    static void setBndryFuncThreadSafety (int ext_dir_safe) { bf_ext_dir_threadsafe = ext_dir_safe; }

    void Print() const;

private:

    IndexType          type;     // Cell centered, node centered ...
    TimeCenter         t_type;   // Temporal centering
    int                id;       // Unique id
    int                ncomp;    // Number of components
    int                ngrow;    // Grow factor
    Interpolater*      mapper;   // Default interpolator
    bool               m_extrap; // Can extrapolate in time?
    bool               m_store_in_checkpoint; // Should store this in the checkpoint file?
    Vector<std::string> names;    // Printable names of components
    Vector<BCRec>       bc;       // Array of bndry types for entire level
    Vector<std::unique_ptr<BndryFunc> >  bc_func;  // Array of pointers to bndry fill functions
    Vector<int>         m_master; // Are we a master or slave? (true or false)
    Vector<int>         m_groupsize; // Groupsize if we're a master
    //
    // If mapper_comp[icomp] != 0, that map is used instead of mapper
    // when interpolating icomp; otherwise, mapper is used.
    //
    Vector<Interpolater*> mapper_comp;
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
    Vector<int> max_map_start_comp;
    Vector<int> min_map_end_comp;

    static int bf_ext_dir_threadsafe;
    static bool bf_thread_safety (const int* lo,const int* hi,
				  const int* dom_lo, const int* dom_hi,
				  const int* bc, int ng);
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
    // The constructor.
    //
    DescriptorList ();
    //
    // Set the list to its default state.
    //
    void clear ();
    //
    // Returns number of elements in the list.
    //
    int size () const;
    //
    // Adds new StateDescriptor at index indx to list.
    //
    void addDescriptor (int                         indx,
                        IndexType                   typ,
                        StateDescriptor::TimeCenter ttyp,
                        int                         nextra,
                        int                         num_comp,
                        Interpolater*               interp,
                        bool                        extrap = false,
                        bool                        store_in_checkpoint = true);
    //
    // Calls resetComponentBCs() on StateDescriptor at index indx.
    //
    void resetComponentBCs (int                               indx,
                            int                               comp,
                            const BCRec&                      bc,
                            const StateDescriptor::BndryFunc& func);
    //
    // Calls setComponent() on StateDescriptor at index indx.
    //
    void setComponent (int                               indx,
                       int                               comp,
                       const std::string&                nm,
                       const BCRec&                      bc,
                       const StateDescriptor::BndryFunc& func,
                       Interpolater*                     interp = 0,
                       int                               max_map_start_comp = -1,
                       int                               min_map_end_comp   = -1);
    //
    // Calls setComponent() on StateDescriptor at index indx on group.
    //
    void setComponent (int                               indx,
                       int                               comp,
                       const Vector<std::string>&         nm,
                       const Vector<BCRec>&               bc,
                       const StateDescriptor::BndryFunc& func,
                       Interpolater*                     interp = 0);
    //
    // Returns StateDescriptor at index k.
    //
    const StateDescriptor& operator[] (int k) const;

    void Print() const;

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    Vector<std::unique_ptr<StateDescriptor> > desc;
};

}

#endif /*_StateDescriptor_H_*/
