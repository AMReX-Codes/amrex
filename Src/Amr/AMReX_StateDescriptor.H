
#ifndef AMREX_StateDescriptor_H_
#define AMREX_StateDescriptor_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterpolater.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBMFInterpolater.H>
#endif

#include <utility>
#include <memory>

namespace amrex {

    using BndryFuncFabDefault = std::function<void(Box const& bx, FArrayBox& data,
                                                   int dcomp, int numcomp,
                                                   Geometry const& geom, Real time,
                                                   const Vector<BCRec>& bcr, int bcomp,
                                                   int scomp)>;

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
        BndryFunc () noexcept = default;

        /**
         * \brief Construct from a scalar boundary-fill callback \p inFunc.
         *
         * \param inFunc Function that fills ghost cells one component at a time.
         */
        BndryFunc (BndryFuncDefault inFunc) noexcept : m_func(inFunc) {}

        /**
         * \brief Construct from a 3-D boundary-fill callback \p inFunc.
         *
         * \param inFunc Function that fills ghost cells using 3-D array arguments.
         */
        BndryFunc (BndryFunc3DDefault inFunc) noexcept : m_func3D(inFunc) {}

        /**
         * \brief Construct from a FAB-aware functor \p inFunc.
         *
         * \param inFunc Callable that operates directly on FArrayBox data.
         */
        BndryFunc (BndryFuncFabDefault inFunc) noexcept : m_funcfab(std::move(inFunc)) {}

        /**
         * \brief Constructor that installs both regular and group boundary callbacks.
         *
         * \param inFunc Function used for per-cell BC fills.
         * \param gFunc  Function used when treating grown regions as a whole.
         */
        BndryFunc (BndryFuncDefault inFunc,BndryFuncDefault gFunc) noexcept
            : m_func(inFunc), m_gfunc(gFunc) {}

        //! 3D counterpart installing per-cell function \p inFunc and group function \p gFunc.
        BndryFunc (BndryFunc3DDefault inFunc,BndryFunc3DDefault gFunc) noexcept
            : m_func3D(inFunc), m_gfunc3D(gFunc) {}

        /**
         * \brief Fill boundary cells using the "regular" scalar function pointer.
         *
         * \param data   Pointer to the raw FAB data.
         * \param lo     Lower bounds of the region being filled.
         * \param hi     Upper bounds of the region being filled.
         * \param dom_lo Lower bounds of the domain.
         * \param dom_hi Upper bounds of the domain.
         * \param dx     Cell spacing array.
         * \param grd_lo Index-space location of the owning grid.
         * \param time   Physical time (as an array for compatibility).
         * \param bc     Encoded boundary condition types.
         */
        void operator () (Real* data, const int* lo, const int* hi,
                          const int* dom_lo, const int* dom_hi,
                          const Real* dx, const Real* grd_lo,
                          const Real* time, const int* bc) const;

        /**
         * \brief Fill boundary cells using the "group" function (supports grown regions).
         *
         * \param data   Pointer to the raw FAB data.
         * \param lo     Lower bounds of the region being filled.
         * \param hi     Upper bounds of the region being filled.
         * \param dom_lo Lower bounds of the domain.
         * \param dom_hi Upper bounds of the domain.
         * \param dx     Cell spacing array.
         * \param grd_lo Index-space location of the owning grid.
         * \param time   Physical time (as an array for compatibility).
         * \param bc     Encoded boundary condition types.
         * \param ng     Number of ghost cells to treat as a group.
         */
        void operator () (Real* data, const int* lo, const int* hi,
                          const int* dom_lo, const int* dom_hi,
                          const Real* dx, const Real* grd_lo,
                          const Real* time, const int* bc, int ng) const;
        /**
         * \brief FAB-aware overload that operates directly on \p data over \p bx at time \p time.
         *
         * \param bx     Region to fill.
         * \param data   FAB storing the components.
         * \param dcomp  Destination component offset in \p data.
         * \param numcomp Number of components to update.
         * \param geom   Geometry describing the level.
         * \param time   Physical time at which BCs are evaluated.
         * \param bcr    Array of BCRecs (one per component).
         * \param bcomp  Starting component offset into \p bcr.
         * \param scomp  Source component offset when referencing BC metadata.
         */
        void operator () (Box const& bx, FArrayBox& data,
                          int dcomp, int numcomp,
                          Geometry const& geom, Real time,
                          const Vector<BCRec>& bcr, int bcomp,
                          int scomp) const;

        [[nodiscard]] bool RunOnGPU () const noexcept { return m_run_on_gpu; }

        void setRunOnGPU (bool b) noexcept { m_run_on_gpu = b; }

        [[nodiscard]] bool hasFabVersion () const noexcept { return m_funcfab != nullptr; }

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
    StateDescriptor () noexcept {} // NOLINT /* = default does not work with old clang */

    /**
    * \brief Constructor that sets all data members.
    *
    * \param btyp    Index type describing cell/nodal staggering.
    * \param ttyp    Time-centering (Point or Interval).
    * \param ident   Identifier assigned to this descriptor.
    * \param nextra  Number of ghost cells.
    * \param num_comp Number of components stored in the state.
    * \param interp  Default interpolater to use between levels.
    * \param extrap  Whether to allow time extrapolation outside stored intervals.
    * \param store_in_checkpoint Should this descriptor be persisted in checkpoints?
    */
    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp,
                     InterpBase*   interp,
                     bool          extrap = false,
                     bool          store_in_checkpoint = true);

    /**
    * \brief Define the data members if constructed with default constructor.
    *
    * \param btyp    Index type describing cell/nodal staggering.
    * \param ttyp    Time-centering (Point or Interval).
    * \param ident   Identifier assigned to this descriptor.
    * \param nextra  Number of ghost cells.
    * \param num_comp Number of components stored in the state.
    * \param interp  Default interpolater to use between levels.
    * \param extrap  Whether to allow time extrapolation outside stored intervals.
    * \param store_in_checkpoint Should this descriptor be persisted in checkpoints?
    */
    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp,
                 InterpBase*   interp,
                 bool          extrap = false,
                 bool          store_in_checkpoint = true);

    /**
    * \brief Configure metadata for component \p comp.
    *
    * \param comp   Component index to configure.
    * \param nm     Name of the component.
    * \param bc     Boundary condition recipe for the component.
    * \param func   Boundary fill functor for physical BCs.
    * \param interp Optional per-component interpolater override.
    * \param max_map_start_comp If >=0, restricts upper bound of remapping region.
    * \param min_map_end_comp   If >=0, restricts lower bound of remapping region.
    */
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       InterpBase*        interp = nullptr,
                       int                max_map_start_comp = -1,
                       int                min_map_end_comp   = -1);

    /**
     * \brief Configure component \p comp using explicit primary/groupsize flags.
     *
     * \param comp      Component index to configure.
     * \param nm        Name of the component.
     * \param bc        Boundary condition recipe for the component.
     * \param func      Boundary fill functor for physical BCs.
     * \param interp    Optional per-component interpolater override.
     * \param a_primary Whether this component is the primary entry for a grouped set.
     * \param groupsize How many coupled components participate in the group anchored here.
     */
    void setComponent (int                comp,
                       const std::string& nm,
                       const BCRec&       bc,
                       const BndryFunc&   func,
                       InterpBase*        interp,
                       bool               a_primary,
                       int                groupsize);

    /**
     * \brief Resets boundary condition data for component \p comp.
     *
     * \param comp Component index to update.
     * \param bcr  Replacement BCRec.
     * \param func Replacement physical BC functor.
     */
    void resetComponentBCs (int              comp,
                            const BCRec&     bcr,
                            const BndryFunc& func);

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
    [[nodiscard]] IndexType getType () const noexcept;

    /**
    * \brief Returns StateDescriptor::TimeCenter.
    */
    [[nodiscard]] StateDescriptor::TimeCenter timeType () const noexcept;

    /**
    * \brief Returns number of components.
    */
    [[nodiscard]] int nComp () const noexcept;

    /**
    * \brief Returns the grow factor.
    */
    [[nodiscard]] int nExtra () const noexcept;

    /**
    * \brief Returns the interpolater.
    */
    [[nodiscard]] InterpBase* interp () const noexcept;

    /**
     * \brief Returns the interpolater used by component \p i.
     *
     * \param i Component index.
     */
    [[nodiscard]] InterpBase* interp (int i) const noexcept;

    /**
     * \brief Returns the name of component \p i.
     *
     * \param i Component index.
     */
    [[nodiscard]] const std::string& name (int i) const noexcept;

    /**
     * \brief Returns the BCRec of component \p i.
     *
     * \param i Component index.
     */
    [[nodiscard]] const BCRec& getBC (int i) const noexcept;

    /**
    * \brief Returns all BCRecs.
    */
    [[nodiscard]] const Vector<BCRec>& getBCs () const noexcept;

    /**
     * \brief Returns the BndryFunc associated with component \p i.
     *
     * \param i Component index.
     */
    [[nodiscard]] const BndryFunc& bndryFill (int i) const noexcept;

    /**
     * \brief Return true if the component range [\p sc, sc+\p nc) lies within this descriptor.
     *
     * \param sc Starting component index.
     * \param nc Number of components in the range.
     */
    [[nodiscard]] int inRange (int sc, int nc) const noexcept;

    /**
     * \brief Return true if all components in [\p scomp, scomp+\p ncomp) share the same interpolater.
     *
     * \param scomp Starting component index.
     * \param ncomp Number of components inspected.
     */
    [[nodiscard]] bool identicalInterps (int scomp, int ncomp) const noexcept;
    //! Returns contiguous component ranges (start,count) that share identical interpolaters within [\p scomp, scomp+\p ncomp).
    [[nodiscard]] std::vector< std::pair<int,int> > sameInterps (int scomp, int ncomp) const;

    /**
    * \brief Can extrapolate in time.
    */
    [[nodiscard]] bool extrap () const noexcept;

    /**
    * \brief Should store this StateData in the checkpoint file
    */
    [[nodiscard]] bool store_in_checkpoint () const noexcept;

    //! Return whether component \p i is the primary entry of a grouped set.
    [[nodiscard]] bool primary (int i) const noexcept { return m_primary[i]; }

    //! Return the size of the component group anchored at primary component \p i.
    [[nodiscard]] int groupsize (int i) const noexcept { return m_groupsize[i]; }


    /**
    * \brief will it run on gpu?
    */
    [[nodiscard]] bool RunOnGPU () const noexcept { return bc_func[0]->RunOnGPU(); }

    /**
    * \brief has new fab version of BndryFunc?
    */
    [[nodiscard]] bool hasBndryFuncFab () const noexcept { return bc_func[0]->hasFabVersion(); }

    //! Control thread-safety flag for boundary fill functors (nonzero \p ext_dir_safe triggers thread-safe handling).
    static void setBndryFuncThreadSafety (int ext_dir_safe) noexcept {
        bf_ext_dir_threadsafe = ext_dir_safe;
    }

private:

    IndexType           type;            //!< Cell centered, node centered ...
    TimeCenter          t_type{Point};   //!< Temporal centering
    int                 id{-1};          //!< Unique id
    int                 ncomp{0};        //!< Number of components
    int                 ngrow{0};        //!< Grow factor
    InterpBase*         mapper{nullptr}; //!< Default interpolator
    bool                m_extrap{false}; //!< Can extrapolate in time?
    bool                m_store_in_checkpoint{true}; //!< Should store this in the checkpoint file?
    Vector<std::string> names;         //!< Printable names of components
    Vector<BCRec>       bc;            //!< Array of bndry types for entire level
    Vector<std::unique_ptr<BndryFunc> >  bc_func; //!< Array of pointers to bndry fill functions
    Vector<int>         m_primary;     //!< Are we a primary or secondary? (true or false)
    Vector<int>         m_groupsize;   //!< Groupsize if we're a primary

    /**
    * \brief If mapper_comp[icomp] != 0, that map is used instead of mapper
    * when interpolating icomp; otherwise, mapper is used.
    */
    Vector<InterpBase*> mapper_comp;

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

    //! Remove all descriptors and reset the list to an empty state.
    void clear ();

    //! Returns number of descriptors currently stored.
    [[nodiscard]] int size () const noexcept;

    /**
     * \brief Adds a new StateDescriptor at index \p indx.
     *
     * \param indx                Slot within the descriptor list (usually StateData index).
     * \param typ                 Index type (cell/nodal).
     * \param ttyp                Time-centering (Point or Interval).
     * \param nextra              Number of ghost cells.
     * \param num_comp            Number of components stored.
     * \param interp              Default interpolater between levels.
     * \param extrap              Whether to allow time extrapolation.
     * \param store_in_checkpoint Whether to persist this descriptor in checkpoints.
     */
    void addDescriptor (int                         indx,
                        IndexType                   typ,
                        StateDescriptor::TimeCenter ttyp,
                        int                         nextra,
                        int                         num_comp,
                        InterpBase*                 interp,
                        bool                        extrap = false,
                        bool                        store_in_checkpoint = true);

    /**
     * \brief Forward resetComponentBCs() to descriptor \p indx.
     *
     * \param indx Descriptor index into the list.
     * \param comp Component to update.
     * \param bc   Replacement BCRec.
     * \param func Replacement physical BC functor.
     */
    void resetComponentBCs (int                               indx,
                            int                               comp,
                            const BCRec&                      bc,
                            const StateDescriptor::BndryFunc& func);

    /**
     * \brief Forward the basic setComponent() overload to descriptor \p indx.
     *
     * \param indx                Descriptor index.
     * \param comp                Component to configure.
     * \param nm                  Component name.
     * \param bc                  Boundary condition recipe.
     * \param func                Physical BC functor.
     * \param interp              Optional per-component interpolater override.
     * \param max_map_start_comp  Optional coupling bound (see StateDescriptor::setComponent()).
     * \param min_map_end_comp    Optional coupling bound.
     */
    void setComponent (int                               indx,
                       int                               comp,
                       const std::string&                nm,
                       const BCRec&                      bc,
                       const StateDescriptor::BndryFunc& func,
                       InterpBase*                       interp = nullptr,
                       int                               max_map_start_comp = -1,
                       int                               min_map_end_comp   = -1);

    /**
     * \brief Batch set multiple grouped components on descriptor \p indx.
     *
     * \param indx   Descriptor index.
     * \param comp   Starting component index for the group.
     * \param nm     Names assigned to each group entry.
     * \param bc     Boundary condition records per entry.
     * \param func   Physical BC functor shared across entries.
     * \param interp Optional interpolater override applied to every component in the group.
     */
    void setComponent (int                               indx,
                       int                               comp,
                       const Vector<std::string>&         nm,
                       const Vector<BCRec>&               bc,
                       const StateDescriptor::BndryFunc& func,
                       InterpBase*                       interp = nullptr);
    //! Returns StateDescriptor stored at index \p k.
    const StateDescriptor& operator[] (int k) const noexcept;

private:

    Vector<std::unique_ptr<StateDescriptor> > desc;
};

}

#endif /*_StateDescriptor_H_*/
