#ifndef AMREX_Derive_H_
#define AMREX_Derive_H_

#include <list>
#include <string>

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Interpolater.H>

namespace amrex {

extern "C"
{

    /**
    * \brief Type of extern "C" function called by DeriveRec to compute derived quantity.
    *
    * Note that AMREX_ARLIM_P will be preprocessed into DIM const int&'s.
    *
    * \param data
    * \param AMREX_ARLIM_P(dlo)
    * \param AMREX_ARLIM_P(dhi)
    * \param nvar
    * \param compdat
    * \param AMREX_ARLIM_P(compdat_lo)
    * \param AMREX_ARLIM_P(compdat_hi)
    * \param ncomp
    * \param lo
    * \param hi
    * \param domain_lo
    * \param domain_hi
    * \param delta
    * \param xlo
    * \param time
    * \param dt
    * \param bcrec
    * \param level
    * \param grid_no
    */
    typedef void (*DeriveFunc)(amrex::Real* data, AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
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
    * \param data
    * \param dlo
    * \param dhi
    * \param nvar
    * \param compdat
    * \param clo
    * \param chi
    * \param ncomp
    * \param lo
    * \param hi
    * \param domain_lo
    * \param domain_hi
    * \param delta
    * \param xlo
    * \param time
    * \param dt
    * \param bcrec
    * \param level
    * \param grid_no
    */
    typedef void (*DeriveFunc3D)(amrex::Real* data, const int* dlo, const int* dhi, const int* nvar,
				 const amrex::Real* compdat, const int* clo, const int* chi, const int* ncomp,
				 const int* lo, const int* hi,
				 const int* domain_lo, const int* domain_hi,
				 const amrex::Real* delta, const amrex::Real* xlo,
				 const amrex::Real* time, const amrex::Real* dt,
				 const int* bcrec,
				 const int* level, const int* grid_no) ;
}

typedef void (*DeriveFuncFab) (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                               const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                               amrex::Real time, const int* bcrec, int level);

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
    typedef Box (*DeriveBoxMap)(const Box&);

    static Box TheSameBox (const Box& box) noexcept;

    static Box GrowBoxByOne (const Box& box) noexcept;

    /**
    * \brief The destructor.
    */
    ~DeriveRec ();

    /**
    * \brief The name of the derived type.
    */
    const std::string& name () const noexcept;

    /**
    * \brief The names of components
    *
    * \param comp
    */
    const std::string& variableName (int comp) const noexcept;

    /**
    * \brief The IndexType of the derived type.
    */
    IndexType deriveType () const noexcept;

    /**
    * \brief The DeriveFunc used to calculate the derived type.
    */
    DeriveFunc    derFunc    () const noexcept;
    DeriveFunc3D  derFunc3D  () const noexcept;
    DeriveFuncFab derFuncFab () const noexcept;

    /**
    * \brief Maps state data box to derived data box.
    */
    DeriveBoxMap boxMap () const noexcept;

    /**
    * \brief Type of interpolater to use in computing derived type.
    */
    Interpolater* interp () const noexcept;

    /**
    * \brief Number of components in the derived type.
    */
    int numDerive () const noexcept;

    /**
    * \brief Number of different chunks of state data needed for derived type.
    */
    int numRange () const noexcept;

    /**
    * \brief Total number of state variables needed for derived type.
    */
    int numState () const noexcept;

    /**
    * \brief The boundary conditions.
    */
    const int* getBC () const noexcept;

    /**
    * \brief Sets state_indx, src_comp and num_comp for the kth
    * range (or chunk) of state data components needed to calculate
    * the derived quantity.
    *
    * \param k
    * \param state_indx
    * \param src_comp
    * \param num_comp
    */
    void getRange (int  k,
                   int& state_indx,
                   int& src_comp,
                   int& num_comp) const;
protected:

    /**
    * \brief Constructor.
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param der_func
    * \param box_map
    * \param interp
    */
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

    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFuncFab  der_func_fab,
               DeriveBoxMap   box_map,
               Interpolater*  interp = &pc_interp);


    /**
    * \brief Constructor without a Fortran function
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param box_map
    */
    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveRec::DeriveBoxMap box_map);

    /**
    * \brief Constructor.
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param var_names
    * \param der_func
    * \param box_map
    * \param interp
    */
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

    DeriveRec (const std::string&  name,
               IndexType           result_type,
               int                 nvar_derive,
	       Vector<std::string>& var_names,
               DeriveFuncFab       der_func_fab,
               DeriveBoxMap        box_map,
               Interpolater*       interp = &pc_interp);

    void addRange (const DescriptorList& d_list,
                   int                   state_indx,
                   int                   src_comp,
                   int                   num_comp);

    void buildBC (const DescriptorList& d_list);

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

    //! Interpolater for mapping crse grid derived data to finer levels.
    Interpolater* mapper = nullptr;

    //! Box mapper that specifies constituent region given derived region.
    DeriveBoxMap bx_map;

    //! Total number of state variables.
    int n_state = 0;

    //! Number of state ranges.
    int nsr = 0;

    //! List of state data subranges required to derive given quantity.
    StateRange* rng = nullptr;

    //! Array of bndry types.
    int* bcr = nullptr;
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
    DeriveList ();

    DeriveList (const DeriveList&) = delete;
    DeriveList& operator= (const DeriveList&) = delete;

    /**
    * \brief Determines whether quantity identified by \<name\> is in the registry.
    *
    * \param name
    */
    bool canDerive (const std::string& name) const;

    /**
    * \brief Access the particular record in registry.
    *
    * \param name
    */
    const DeriveRec* get (const std::string& name) const;

    /**
    * \brief Adds another entry to the registry.
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param der_func
    * \param box_map
    * \param interp
    */
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

    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFuncFab           der_func_fab,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);


    /**
    * \brief This version doesn't take a Fortran function.
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param box_map
    */
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveRec::DeriveBoxMap box_map = &DeriveRec::TheSameBox );

    /**
    * \brief Adds another entry to the registry.
    *
    * \param name
    * \param result_type
    * \param nvar_derive
    * \param var_names
    * \param der_func
    * \param box_map
    * \param interp
    */
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              Vector<std::string>&    var_names,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);

    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              Vector<std::string>&    var_names,
              DeriveFunc3D            der_func_3d,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);

    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              Vector<std::string>&    var_names,
              DeriveFuncFab           der_func_fab,
              DeriveRec::DeriveBoxMap box_map,
              Interpolater*           interp = &pc_interp);

    /**
    * \brief Adds another StateRange to the DeriveRec identified by \<name\>.
    *
    * \param name
    * \param d_list
    * \param state_indx
    * \param start_comp
    * \param ncomp
    */
    void addComponent (const std::string&    name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   start_comp,
                       int                   ncomp);

    std::list<DeriveRec>& dlist ();

    void clear () { lst.clear(); }

private:

    std::list<DeriveRec> lst;
};

}

#endif /*_Derive_H_*/
