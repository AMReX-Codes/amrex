#ifndef AMREX_Derive_H_
#define AMREX_Derive_H_
#include <AMReX_Config.H>

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

/**
 * \file AMReX_Derive.H
 * \brief Infrastructure for registering and computing derived quantities from state data.
 */

#include <list>
#include <string>

namespace amrex {

extern "C"
{

    /**
    * \brief Type of extern "C" function called by DeriveRec to compute derived quantity.
    *
    * Note that AMREX_ARLIM_P will be preprocessed into DIM const int&'s.
    *
    * \param data Pointer to the destination data slab that will store the derived quantity.
    * \param AMREX_ARLIM_P(dlo) Lower index bounds of \p data (macro expands to i/j/k references).
    * \param AMREX_ARLIM_P(dhi) Upper index bounds of \p data.
    * \param nvar Number of components laid out in \p data (destination FAB).
    * \param compdat Pointer to the source state data for the requested component set.
    * \param AMREX_ARLIM_P(compdat_lo) Lower bounds of the source data buffer.
    * \param AMREX_ARLIM_P(compdat_hi) Upper bounds of the source data buffer.
    * \param ncomp Number of source components provided through \p compdat.
    * \param lo Cell-space lower bounds of the target region currently being filled.
    * \param hi Cell-space upper bounds (inclusive) of that region.
    * \param domain_lo Lower bounds of the problem domain (helps apply BCs).
    * \param domain_hi Upper bounds of the problem domain.
    * \param delta Cell spacing vector (dx,dy,dz).
    * \param xlo Physical-space origin of the owning grid.
    * \param time Physical time at which the derived quantity should be evaluated.
    * \param dt Time step associated with the state data (old-to-new interval).
    * \param bcrec Encoded boundary condition types for each component.
    * \param level AMR level index containing the data.
    * \param grid_no Grid identifier within that level.
    */
    using DeriveFunc = void (*)(amrex::Real* data, AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
                                const int* nvar, const amrex::Real* compdat,
                                AMREX_ARLIM_P(compdat_lo), AMREX_ARLIM_P(compdat_hi),
                                const int* ncomp,
                                const int* lo, const int* hi,
                                const int* domain_lo, const int* domain_hi,
                                const amrex::Real* delta, const amrex::Real* xlo,
                                const amrex::Real* time, const amrex::Real* dt,
                                const int* bcrec,
                                const int* level, const int* grid_no) ;

    /**
    * \brief This is dimension agnostic.  For example, dlo always has three elements.
    *
    * \param data Pointer to the destination FAB data.
    * \param dlo Lower index bounds of the destination FAB.
    * \param dhi Upper index bounds of the destination FAB.
    * \param nvar Number of destination components.
    * \param compdat Pointer to the source state data.
    * \param clo Lower index bounds of the source buffer.
    * \param chi Upper index bounds of the source buffer.
    * \param ncomp Number of source components available.
    * \param lo Cell-space lower bounds of the region being filled.
    * \param hi Cell-space upper bounds of that region.
    * \param domain_lo Lower bounds of the domain (for BC handling).
    * \param domain_hi Upper bounds of the domain.
    * \param delta Cell spacing vector (dx,dy,dz).
    * \param xlo Physical origin of the owning grid.
    * \param time Physical time at which the derived quantity is evaluated.
    * \param dt Time step between the stored state slots.
    * \param bcrec Encoded boundary condition identifiers.
    * \param level AMR level index.
    * \param grid_no Grid identifier within that level.
    */
    using DeriveFunc3D = void (*)(amrex::Real* data, const int* dlo, const int* dhi, const int* nvar,
                                  const amrex::Real* compdat, const int* clo, const int* chi, const int* ncomp,
                                  const int* lo, const int* hi,
                                  const int* domain_lo, const int* domain_hi,
                                  const amrex::Real* delta, const amrex::Real* xlo,
                                  const amrex::Real* time, const amrex::Real* dt,
                                  const int* bcrec,
                                  const int* level, const int* grid_no) ;
}

using DeriveFuncFab = std::function<void(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                                         const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                                         amrex::Real time, const int* bcrec, int level)>;

using DeriveFuncMF = std::function<void(amrex::MultiFab& der_mf, int dcomp, int ncomp,
                                        const amrex::MultiFab& data_mf, const amrex::Geometry& geomdata,
                                        amrex::Real time, const int* bcrec, int level)>;

class DescriptorList;


/**
* \brief Derived Type Record
*
* Computes quantities derived from state data.
*
* DeriveRec is designed to compute quantities which can be derived
* from the state data contained in AmrLevel and its derivatives. Some
* examples might be kinetic energy, vorticity, concentration gradients ...
*/
class DeriveRec
{
   friend class DeriveList;

public:

    /**
    * \brief A pointer to function taking and returning a Box.
    *
    */
    using DeriveBoxMap = std::function<Box(const Box&)>;

    /**
     * \brief Identity mapping from destination to source boxes.
     *
     * \param box Destination box passed to the callback.
     */
    static Box TheSameBox (const Box& box) noexcept;

    /**
     * \brief Convenience mapper that grows the supplied box by one cell in every direction.
     *
     * \param box Destination box; the return value is `grow(box,1)`.
     */
    static Box GrowBoxByOne (const Box& box) noexcept;

    /**
    * \brief The destructor.
    */
    ~DeriveRec ();

    //! Non-copyable; derived metadata owns non-trivial resources.
    DeriveRec (DeriveRec const&) = delete;
    //! Non-movable to keep internal pointers/registrations stable.
    DeriveRec (DeriveRec &&) = delete;
    //! Copy assignment disabled for the same reason as the copy constructor.
    DeriveRec& operator= (DeriveRec const&) = delete;
    //! Move assignment disabled to avoid dangling references inside registries.
    DeriveRec& operator= (DeriveRec &&) = delete;

    /**
    * \brief The name of the derived type.
    */
    [[nodiscard]] const std::string& name () const noexcept;

    /**
    * \brief Return the name of derived component \p comp.
    */
    [[nodiscard]] const std::string& variableName (int comp) const noexcept;

    /**
    * \brief The IndexType of the derived type.
    */
    [[nodiscard]] IndexType deriveType () const noexcept;

    /**
    * \brief The DeriveFunc used to calculate the derived type.
    */
    [[nodiscard]] DeriveFunc    derFunc    () const noexcept;
    //! Dimension-agnostic version of derFunc() that always takes 3-component bounds.
    [[nodiscard]] DeriveFunc3D  derFunc3D  () const noexcept;
    //! Callable that operates on a FAB-sized box instead of raw pointers.
    [[nodiscard]] DeriveFuncFab derFuncFab () const noexcept;
    //! MultiFab variant that fills directly without per-box callbacks.
    [[nodiscard]] DeriveFuncMF  derFuncMF () const noexcept;

    /**
    * \brief Maps state data box to derived data box.
    */
    [[nodiscard]] DeriveBoxMap boxMap () const noexcept;

    /**
    * \brief Type of interpolater to use in computing derived type.
    */
    [[nodiscard]] Interpolater* interp () const noexcept;

    /**
    * \brief Number of components in the derived type.
    */
    [[nodiscard]] int numDerive () const noexcept;

    /**
    * \brief Number of different chunks of state data needed for derived type.
    */
    [[nodiscard]] int numRange () const noexcept;

    /**
    * \brief Total number of state variables needed for derived type.
    */
    [[nodiscard]] int numState () const noexcept;

    /**
    * \brief The boundary conditions.
    */
    [[nodiscard]] const int* getBC () const noexcept;
    //! Boundary-condition array padded/truncated to three entries.
    [[nodiscard]] const int* getBC3D () const noexcept;

    /**
    * \brief Describe which chunk of state data is needed for derived component range \p k.
    *
    * \param k           Range index in [0,numRange()).
    * \param state_indx  Output: index into the state descriptor array.
    * \param src_comp    Output: starting component in the chosen state.
    * \param num_comp    Output: number of components to read.
    */
    void getRange (int  k,
                   int& state_indx,
                   int& src_comp,
                   int& num_comp) const;

    /**
     * \brief Construct a record driven by a Fortran-style DeriveFunc.
     *
     * \param name         Derived quantity name (unique registry key).
     * \param result_type  Index type describing the derived FAB layout.
     * \param nvar_derive  Number of components produced by the derive function.
     * \param der_func     Function pointer invoked to compute the quantity.
     * \param box_map      Mapping from destination boxes back to source boxes.
     * \param interp       Interpolater used when populating finer levels (defaults to piecewise-constant).
     */
    DeriveRec (std::string    name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc     der_func,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);

    /**
     * \brief Construct a record using a dimension-agnostic DeriveFunc3D callback.
     *
     * Parameters mirror the DeriveFunc constructor but \p der_func_3d follows the 3-D argument shape.
     */
    DeriveRec (std::string    name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc3D   der_func_3d,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);

    /**
     * \brief Construct a record whose implementation operates directly on FABs.
     *
     * Parameters match the scalar constructor but \p der_func_fab is a callable accepting FABs.
     */
    DeriveRec (std::string    name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFuncFab  der_func_fab,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);

    /**
     * \brief Construct a record whose implementation operates on MultiFabs.
     *
     * Parameters match the scalar constructor but \p der_func_mf takes MultiFabs.
     */
    DeriveRec (std::string    name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFuncMF   der_func_mf,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);


    /**
     * \brief Construct a placeholder record without attaching a callable yet.
     *
     * \param name        Derived quantity name.
     * \param result_type Index type describing the derived FAB layout.
     * \param nvar_derive Number of components produced.
     * \param box_map     Mapping from destination boxes back to source boxes.
     */
    DeriveRec (std::string    name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveRec::DeriveBoxMap box_map);

    /**
     * \brief Same as the basic constructor but labels each derived component.
     *
     * \param name         Derived quantity name.
     * \param result_type  Index type describing the derived FAB layout.
     * \param nvar_derive  Number of components produced.
     * \param var_names    Human-readable names for each component.
     * \param der_func     Fortran-style callback.
     * \param box_map      Mapping from destination to source boxes.
     * \param interp       Interpolater override (defaults to piecewise-constant).
     */
    DeriveRec (std::string         name,
               IndexType           result_type,
               int                 nvar_derive,
               Vector<std::string> const& var_names,
               DeriveFunc          der_func,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

    //! 3-D callback variant (parameters mirror the previous constructor).
    DeriveRec (std::string         name,
               IndexType           result_type,
               int                 nvar_derive,
               Vector<std::string> const& var_names,
               DeriveFunc3D        der_func_3d,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

    //! FAB callback variant (parameters mirror the previous constructor).
    DeriveRec (std::string         name,
               IndexType           result_type,
               int                 nvar_derive,
               Vector<std::string> const& var_names,
               DeriveFuncFab       der_func_fab,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

    //! MultiFab callback variant (parameters mirror the previous constructor).
    DeriveRec (std::string         name,
               IndexType           result_type,
               int                 nvar_derive,
               Vector<std::string> const& var_names,
               DeriveFuncMF        der_func_mf,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

    /**
     * \brief Append another chunk of state data requirements.
     *
     * \param d_list      Descriptor registry providing metadata.
     * \param state_indx  Index into the descriptor list.
     * \param src_comp    Starting component inside the chosen state descriptor.
     * \param num_comp    Number of components to sample for this range.
     */
    void addRange (const DescriptorList& d_list,
                   int                   state_indx,
                   int                   src_comp,
                   int                   num_comp);

    //! Build the 1D/2D BC arrays from the registered state descriptors \p d_list.
    void buildBC (const DescriptorList& d_list);
    //! Build the 3D BC arrays from the registered state descriptors \p d_list.
    void buildBC3D (const DescriptorList& d_list);

    //! An element of a linked list to point to state quantities in AmrLevels.
    struct StateRange
    {
        int         typ;
        int         sc;
        int         nc;
        StateRange* next;
    };

private:

    //! Name of derived quantity.
    std::string derive_name;

    //! Names of derived variables
    Vector<std::string> variable_names;

    //! Type of derived quantity.
    IndexType der_type;

    //! Number of components in derived quantity.
    int n_derive = 0;

    //! Function that computes derived quantity from state variables.
    DeriveFunc    func = nullptr;
    DeriveFunc3D  func_3d = nullptr;
    DeriveFuncFab func_fab = nullptr;
    DeriveFuncMF  func_mf = nullptr;

    //! Interpolater for mapping crse grid derived data to finer levels.
    Interpolater* mapper = nullptr;

    //! Box mapper that specifies constituent region given derived region.
    DeriveBoxMap bx_map = nullptr;

    //! Total number of state variables.
    int n_state = 0;

    //! Number of state ranges.
    int nsr = 0;

    //! List of state data subranges required to derive given quantity.
    StateRange* rng = nullptr;

    //! Array of bndry types.
    int* bcr = nullptr;
    int* bcr3D = nullptr;
};


/**
* \brief A list of DeriveRecs.
*
* DeriveList manages and provides access to the list of DeriveRecs.
*/
class DeriveList
{
public:

    /**
    * \brief The default constructor.
    */
    DeriveList () = default;

    //! Destructor (defaulted to keep header-only implementation trivial).
    ~DeriveList () = default;
    //! Non-copyable registry; it owns pointers to callback records.
    DeriveList (const DeriveList&) = delete;
    //! Non-movable because the registry’s entries self-register with static globals.
    DeriveList (DeriveList &&) = delete;
    //! Copy assignment disabled for the same reason as the copy constructor.
    DeriveList& operator= (const DeriveList&) = delete;
    //! Move assignment disabled to avoid dangling pointers inside DeriveRec entries.
    DeriveList& operator= (DeriveList &&) = delete;

    /**
     * \brief Determines whether quantity identified by \p name is registered.
     *
     * \param name Derived quantity name.
     */
    [[nodiscard]] bool canDerive (const std::string& name) const;

    /**
     * \brief Access the particular record in registry.
     *
     * \param name Derived quantity name.
     */
    [[nodiscard]] const DeriveRec* get (const std::string& name) const;

    /**
     * \brief Register a new derived quantity backed by a DeriveFunc.
     *
     * \param name        Derived quantity name.
     * \param result_type Index type describing the derived FAB layout.
     * \param nvar_derive Number of components produced.
     * \param der_func    Function pointer to invoke.
     * \param bx_map      Mapping from destination boxes to required source boxes.
     * \param interp      Interpolater used on finer levels (defaults to piecewise-constant).
     */
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc              der_func,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! DeriveFunc3D variant (parameters mirror the DeriveFunc version).
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc3D            der_func_3d,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! DeriveFuncFab variant (parameters mirror the DeriveFunc version).
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              const DeriveFuncFab&           der_func_fab,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! DeriveFuncMF variant (parameters mirror the DeriveFunc version).
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              const DeriveFuncMF&           der_func_mf,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    /**
     * \brief Register a placeholder record without attaching a callable yet.
     *
     * \param name        Derived quantity name.
     * \param result_type Index type describing the derived FAB layout.
     * \param nvar_derive Number of components produced.
     * \param box_map     Mapping from destination boxes to required source boxes.
     */
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              const DeriveRec::DeriveBoxMap& box_map = &DeriveRec::TheSameBox );

    /**
     * \brief Register a derived quantity and supply per-component names up front.
     *
     * \param name        Derived quantity name.
     * \param res_typ     Index type describing the derived FAB layout.
     * \param nvar_derive Number of components produced.
     * \param vars        Component names (size \p nvar_derive).
     * \param der_func    Function pointer to invoke.
     * \param bx_map      Mapping from destination boxes to source boxes.
     * \param interp      Interpolater used on finer levels.
     */
    void add (const std::string&      name,
              IndexType               res_typ,
              int                     nvar_derive,
              Vector<std::string> const&    vars,
              DeriveFunc              der_func,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! Named-component variant using DeriveFunc3D callbacks.
    void add (const std::string&      name,
              IndexType               res_typ,
              int                     nvar_derive,
              Vector<std::string> const&    vars,
              DeriveFunc3D            der_func_3d,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! Named-component variant using DeriveFuncFab callbacks.
    void add (const std::string&      name,
              IndexType               res_typ,
              int                     nvar_derive,
              Vector<std::string> const&    vars,
              const DeriveFuncFab&           der_func_fab,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    //! Named-component variant using DeriveFuncMF callbacks.
    void add (const std::string&      name,
              IndexType               res_typ,
              int                     nvar_derive,
              Vector<std::string> const&    vars,
              const DeriveFuncMF&           der_func_mf,
              const DeriveRec::DeriveBoxMap& bx_map,
              Interpolater*           interp = &pc_interp);

    /**
     * \brief Adds another StateRange to the DeriveRec identified by \p name.
     *
     * \param name       Derived quantity being augmented.
     * \param d_list     Descriptor registry providing metadata.
     * \param state_indx State descriptor index supplying data.
     * \param s_comp     Starting component within that descriptor.
     * \param n_comp     Number of components included in the range.
     */
    void addComponent (const std::string&    name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   s_comp,
                       int                   n_comp);

    //! Mutable access to the underlying list of records.
    [[nodiscard]] std::list<DeriveRec>& dlist ();

    //! Remove all registered derived records.
    void clear () { lst.clear(); }

private:

    std::list<DeriveRec> lst;
};

}

#endif /*_Derive_H_*/
