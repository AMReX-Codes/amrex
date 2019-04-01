
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

typedef void (*BndryFuncFabDefault) (Box const& bx, FArrayBox& data,
                                     const int dcomp, const int numcomp,
                                     Geometry const& geom, const Real time,
                                     const Vector<BCRec>& bcr, const int bcomp,
                                     const int scomp);

class Interpolater;


/**
* \brief Attributes of StateData.
*/

class StateDescriptor
{
    friend class DescriptorList;

public:

    //! Defines the temporal centeredness of data.
    enum TimeCenter { Point = 0 , Interval };

    /**
    * Type of function called by BCRec for user-supplied boundary data.
    *
    * Class wrapper around BndryFuncDefault.
    */
    class BndryFunc
    {
    public:

        /**
        * \brief Bogus constructor.
        */
        BndryFunc () noexcept {}

        /**
        * \brief A Constructor.
        */
        BndryFunc (BndryFuncDefault inFunc) noexcept : m_func(inFunc) {}

        BndryFunc (BndryFunc3DDefault inFunc) noexcept : m_func3D(inFunc) {}

        BndryFunc (BndryFuncFabDefault inFunc) noexcept : m_funcfab(inFunc) {}

        /**
        * \brief Another Constructor.
        *
        * \param inFunc
        * \param gFunc
        */
        BndryFunc (BndryFuncDefault inFunc,BndryFuncDefault gFunc) noexcept
            : m_func(inFunc), m_gfunc(gFunc) {}

        BndryFunc (BndryFunc3DDefault inFunc,BndryFunc3DDefault gFunc) noexcept
            : m_func3D(inFunc), m_gfunc3D(gFunc) {}

        /**
        * \brief Return a ptr to a clone of this object.
        * It is the responsibility of the caller to delete the result.
        */
        virtual BndryFunc* clone () const;

        /**
        * \brief Destructor.
        */
        virtual ~BndryFunc ();

	/**
	* Fill boundary cells using "regular" function.
    *
	* \param data
	* \param lo
	* \param hi
	* \param dom_lo
	* \param dom_hi
	* \param dx
	* \param grd_lo
	* \param time
	* \param bc
	*/
	virtual void operator () (Real* data, const int* lo, const int* hi,
				  const int* dom_lo, const int* dom_hi,
				  const Real* dx, const Real* grd_lo,
				  const Real* time, const int* bc) const;

        /**
        * \brief Fill boundary cells using "group" function.
        *
        * \param data
        * \param lo
        * \param hi
        * \param dom_lo
        * \param dom_hi
        * \param dx
        * \param grd_lo
        * \param time
        * \param bc
        * \param ng
        */
        virtual void operator () (Real* data, const int* lo, const int* hi,
                                  const int* dom_lo, const int* dom_hi,
                                  const Real* dx, const Real* grd_lo,
                                  const Real* time, const int* bc, int ng) const;

        virtual void operator () (Box const& bx, FArrayBox& data,
                                  const int dcomp, const int numcomp,
                                  Geometry const& geom, const Real time,
                                  const Vector<BCRec>& bcr, const int bcomp,
                                  const int scomp) const;

        bool RunOnGPU () const noexcept { return m_run_on_gpu; }

        void setRunOnGPU (bool b) noexcept { m_run_on_gpu = b; }

        bool hasFabVersion () const noexcept { return m_funcfab != nullptr; }

    private:

        BndryFuncDefault   m_func = nullptr;
        BndryFunc3DDefault m_func3D = nullptr;
        BndryFuncDefault   m_gfunc = nullptr;
        BndryFunc3DDefault m_gfunc3D = nullptr;
        BndryFuncFabDefault m_funcfab = nullptr;
        bool m_run_on_gpu = false;
    };

    /**
    * \brief The default constructor.
    */
    StateDescriptor () noexcept;

    /**
    * \brief Constructor that sets all data members.
    *
    * \param btyp
    * \param ttyp
    * \param ident
    * \param nextra
    * \param num_comp
    * \param interp
    * \param extrap
    * \param store_in_checkpoint
    */
    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     Interpolater* interp,
                     bool          extrap = false,
                     bool          store_in_checkpoint = true);

    /**
    * \brief The destructor.
    */
    ~StateDescriptor ();

    /**
    * \brief Define the data members if constructed with default construtor.
    *
    * \param btyp
    * \param ttyp
    * \param ident
    * \param nextra
    * \param num_comp
    * \param interp
    * \param extrap
    * \param store_in_checkpoint
    */
    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 Interpolater* interp,
                 bool          extrap = false,
                 bool          store_in_checkpoint = true);

    /**
    * \brief Sets details of comp component.
    *
    * \param comp
    * \param nm
    * \param bc
    * \param func
    * \param interp
    * \param max_map_start_comp
    * \param min_map_end_comp
    */
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp = 0,
                       int                max_map_start_comp = -1,
                       int                min_map_end_comp   = -1);

    /**
    * \brief Sets details of comp component.
    *
    * \param comp
    * \param nm
    * \param bc
    * \param func
    * \param interp
    * \param master_or_slave
    * \param groupsize
    */
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       Interpolater*      interp,
                       bool               master_or_slave,
                       int                groupsize);

    /**
    * \brief Resets boundary condition data for comp component.
    *
    * \param comp
    * \param bcr
    * \param func
    */
    void resetComponentBCs (int              comp,
                            const BCRec&     bcr,
                            const BndryFunc& func);

    /**
    * \brief Set interpolaters for a subset of the state vector components.
    *
    * \param use_default_map
    * \param default_map
    * \param start_comp
    * \param num_comp
    * \param maps
    * \param nmaps
    * \param map_start_comp
    * \param map_num_comp
    * \param max_start_comp
    * \param min_end_comp
    */
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

    /**
    * \brief Cleanup interpolaters for a subset of the state vector components.
    *
    * \param maps
    * \param map_start_comp
    * \param map_num_comp
    * \param max_start_comp
    * \param min_end_comp
    */
    void cleanUpMaps (Interpolater**& maps,
                      int*&           map_start_comp,
                      int*&           map_num_comp,
                      int*&           max_start_comp,
                      int*&           min_end_comp) const;

    /**
    * \brief Output names of components.
    *
    * \param os
    * \param start_comp
    * \param num_comp
    */
    void dumpNames (std::ostream& os,
                    int           start_comp,
                    int           num_comp) const;

    /**
    * \brief Returns the IndexType.
    */
    IndexType getType () const noexcept;

    /**
    * \brief Returns StateDescriptor::TimeCenter.
    */
    StateDescriptor::TimeCenter timeType () const noexcept;

    /**
    * \brief Returns number of components.
    */
    int nComp () const noexcept;

    /**
    * \brief Returns the grow factor.
    */
    int nExtra () const noexcept;

    /**
    * \brief Returns the interpolater.
    */
    Interpolater* interp () const noexcept;

    /**
    * \brief Returns the interpolater of specified component.
    *
    * \param i
    */
    Interpolater* interp (int i) const noexcept;

    /**
    * \brief Returns the name of specified component.
    *
    * \param i
    */
    const std::string& name (int i) const noexcept;

    /**
    * \brief Returns the BCRec of specified component.
    *
    * \param i
    */
    const BCRec& getBC (int i) const noexcept;

    /**
    * \brief Returns all BCRecs.
    */
    const Vector<BCRec>& getBCs () const noexcept;

    /**
    * \brief Returns the BndryFunc of specified component.
    *
    * \param i
    */
    const BndryFunc& bndryFill (int i) const noexcept;

    /**
    * \brief Is sc\>=0 \&\& sc+nc\<=ncomp ?
    *
    * \param sc
    * \param nc
    */
    int inRange (int sc, int nc) const noexcept;

    /**
    * \brief Are the interpolaters in the specified range identical?
    *
    * \param scomp
    * \param ncomp
    */
    bool identicalInterps (int scomp, int ncomp) const noexcept;
    //
    // Returns contiguous ranges of comps with identical interpolaters.
    //
    std::vector< std::pair<int,int> > sameInterps (int scomp, int ncomp) const;

    /**
    * \brief Can extrapolate in time.
    */
    bool extrap () const noexcept;

    /**
    * \brief Should store this StateData in the checkpoint file
    */
    bool store_in_checkpoint () const noexcept;

    bool master (int i) const noexcept { return m_master[i]; }

    int groupsize (int i) const noexcept { return m_groupsize[i]; }


    /**
    * \brief will it run on gpu?
    */
    bool RunOnGPU () const noexcept { return bc_func[0]->RunOnGPU(); }

    /**
    * \brief has new fab version of BndryFunc?
    */
    bool hasBndryFuncFab () const noexcept { return bc_func[0]->hasFabVersion(); }

    static void setBndryFuncThreadSafety (int ext_dir_safe) noexcept {
        bf_ext_dir_threadsafe = ext_dir_safe;
    }

private:

    IndexType          type;         //!< Cell centered, node centered ...
    TimeCenter         t_type;       //!< Temporal centering
    int                id;           //!< Unique id
    int                ncomp;        //!< Number of components
    int                ngrow;        //!< Grow factor
    Interpolater*      mapper;       //!< Default interpolator
    bool               m_extrap;     //!< Can extrapolate in time?
    bool               m_store_in_checkpoint;     //!< Should store this in the checkpoint file?
    Vector<std::string> names;        //!< Printable names of components
    Vector<BCRec>       bc;           //!< Array of bndry types for entire level
    Vector<std::unique_ptr<BndryFunc> >  bc_func;      //!< Array of pointers to bndry fill functions
    Vector<int>         m_master;     //!< Are we a master or slave? (true or false)
    Vector<int>         m_groupsize;     //!< Groupsize if we're a master

    /**
    * \brief If mapper_comp[icomp] != 0, that map is used instead of mapper
    * when interpolating icomp; otherwise, mapper is used.
    */
    Vector<Interpolater*> mapper_comp;

    /**
    * \brief For component icomp, max_map_start_comp[icomp] and
    * min_map_end_comp[icomp] represent the maximum starting index and
    * minimum ending index with which to use mapper_comp[icomp] when
    * interpolating component icomp.  These are convenient for using
    * with "coupled" multiple component maps like CellConservativeLinear.
    *
    * If there is no need to "couple" a component, max_map_start_comp[icomp]
    * and min_map_end_comp[icomp] simply equal icomp.
    */
    Vector<int> max_map_start_comp;
    Vector<int> min_map_end_comp;

    static int bf_ext_dir_threadsafe;

    static bool bf_thread_safety (const int* lo,const int* hi,
				  const int* dom_lo, const int* dom_hi,
				  const int* bc, int ng);
};

/**
* A List of StateDescriptors.
*
* A container class for StateDescriptors.
*/

class DescriptorList
{
public:

    /**
    * \brief The constructor.
    */
    DescriptorList () noexcept;

    DescriptorList (const DescriptorList&) = delete;
    DescriptorList& operator= (const DescriptorList&) = delete;

    /**
    * \brief Set the list to its default state.
    */
    void clear ();

    /**
    * \brief Returns number of elements in the list.
    */
    int size () const noexcept;

    /**
    * \brief Adds new StateDescriptor at index indx to list.
    *
    * \param indx
    * \param typ
    * \param ttyp
    * \param nextra
    * \param num_comp
    * \param interp
    * \param extrap
    * \param store_in_checkpoint
    */
    void addDescriptor (int                         indx,
                        IndexType                   typ,
                        StateDescriptor::TimeCenter ttyp,
                        int                         nextra,
                        int                         num_comp,
                        Interpolater*               interp,
                        bool                        extrap = false,
                        bool                        store_in_checkpoint = true);

    /**
    * \brief Calls resetComponentBCs() on StateDescriptor at index indx.
    *
    * \param indx
    * \param comp
    * \param bc
    * \param func
    */
    void resetComponentBCs (int                               indx,
                            int                               comp,
                            const BCRec&                      bc,
                            const StateDescriptor::BndryFunc& func);

    /**
    * \brief Calls setComponent() on StateDescriptor at index indx.
    *
    * \param indx
    * \param comp
    * \param nm
    * \param bc
    * \param func
    * \param interp
    * \param max_map_start_comp
    * \param min_map_end_comp
    */
    void setComponent (int                               indx,
                       int                               comp,
                       const std::string&                nm,
                       const BCRec&                      bc,
                       const StateDescriptor::BndryFunc& func,
                       Interpolater*                     interp = 0,
                       int                               max_map_start_comp = -1,
                       int                               min_map_end_comp   = -1);

    /**
    * \brief Calls setComponent() on StateDescriptor at index indx on group.
    *
    * \param indx
    * \param comp
    * \param nm
    * \param bc
    * \param func
    * \param interp
    */
    void setComponent (int                               indx,
                       int                               comp,
                       const Vector<std::string>&         nm,
                       const Vector<BCRec>&               bc,
                       const StateDescriptor::BndryFunc& func,
                       Interpolater*                     interp = 0);
    /**
    * Returns StateDescriptor at index k.
    */
    const StateDescriptor& operator[] (int k) const noexcept;

private:

    Vector<std::unique_ptr<StateDescriptor> > desc;
};

}

#endif /*_StateDescriptor_H_*/
