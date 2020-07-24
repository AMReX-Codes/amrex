#ifndef AMREX_ErrorList_H_
#define AMREX_ErrorList_H_

#include <string>
#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_MultiFab.H>
#include <AMReX_TagBox.H>
#include <AMReX_Geometry.H>

namespace amrex {


extern "C"
{

    /**
    * \brief Type of extern "C" function called by ErrorRec to do tagging of cells for refinement.
    *
    * \param tag
    * \param tlo
    * \param thi
    * \param tagval
    * \param clearval
    * \param data
    * \param data_lo
    * \param data_hi
    * \param lo
    * \param hi
    * \param nvar
    * \param domain_lo
    * \param domain_hi
    * \param dx
    * \param xlo
    * \param prob_lo
    * \param time
    * \param level
    */
    typedef void (*ErrorFuncDefault)(int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                     const int* tagval, const int* clearval,
                                     amrex::Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                     const int* lo, const int * hi, const int* nvar,
                                     const int* domain_lo, const int* domain_hi,
                                     const amrex::Real* dx, const amrex::Real* xlo,
                                     const amrex::Real* prob_lo, const amrex::Real* time,
                                     const int* level);

    typedef void (*ErrorFunc2Default)(int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                      const int* tagval, const int* clearval,
                                      amrex::Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                      const int* lo, const int * hi, const int* nvar,
                                      const int* domain_lo, const int* domain_hi,
                                      const amrex::Real* dx, const int* level, const amrex::Real* avg);


    /**
    * \brief Dimension agnostic version that always has three elements.
    * Note that this is only implemented for the ErrorFunc class, not ErrorFunc2.
    *
    * \param tag
    * \param tlo
    * \param thi
    * \param tagval
    * \param clearval
    * \param data
    * \param data_lo
    * \param data_hi
    * \param lo
    * \param hi
    * \param nvar
    * \param domain_lo
    * \param domain_hi
    * \param dx
    * \param xlo
    * \param prob_lo
    * \param time
    * \param level
    */
    typedef void (*ErrorFunc3DDefault)(int* tag, const int* tlo, const int* thi,
                                       const int* tagval, const int* clearval,
         			       amrex::Real* data, const int* data_lo, const int* data_hi,
                                       const int* lo, const int * hi, const int* nvar,
                                       const int* domain_lo, const int* domain_hi,
                                       const amrex::Real* dx, const amrex::Real* xlo,
                                       const amrex::Real* prob_lo, const amrex::Real* time,
                                       const int* level);

}


/**
* \brief Error Record.
*
*
* ErrorRec is designed to tag cells for regridding based on the state
* data contained in AmrLevel and its derivatives.  It is conceivable that
* actual error tagging will be through derivation, so provision is made
* for this as well.
*/

class ErrorRec
{
public:
    //
    // Error types.
    //
    enum ErrorType { Special=0, Standard, UseAverage };
    //
    // Class wrapper around ErrorFuncDefault.
    //
    class ErrorFunc
    {
    public:

        /**
        * \brief Bogus constructor.
        */
        ErrorFunc ();

        /**
        * \brief A Constructor.
        *
        * \param inFunc
        */
        // No explicit because some codes rely on it.
        ErrorFunc (ErrorFuncDefault inFunc);

        /**
        * \brief A Constructor.
        *
        * \param inFunc
        */
        // No explicit because some codes rely on it.
        ErrorFunc (ErrorFunc3DDefault inFunc);

        /**
        * \brief Return a ptr to a clone of this object.
        * It is the responsibility of the caller to delete the result.
        */
        virtual ErrorFunc* clone () const;

        /**
        * \brief Destructor.
        */
        virtual ~ErrorFunc ();

        /**
        * \brief Tag cells using "regular" function.
        *
        * \param tag
        * \param tlo
        * \param thi
        * \param tagval
        * \param clearval
        * \param data
        * \param data_lo
        * \param data_hi
        * \param lo
        * \param hi
        * \param nvar
        * \param domain_lo
        * \param domain_hi
        * \param dx
        * \param xlo
        * \param prob_lo
        * \param time
        * \param level
        */
        virtual void operator () (int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                  const int* tagval, const int* clearval,
                                  Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const;

        /**
        * \brief Tag cells using dimension-agnostic "regular" function.
        *
        * \param tag
        * \param tlo
        * \param thi
        * \param tagval
        * \param clearval
        * \param data
        * \param data_lo
        * \param data_hi
        * \param lo
        * \param hi
        * \param nvar
        * \param domain_lo
        * \param domain_hi
        * \param dx
        * \param xlo
        * \param prob_lo
        * \param time
        * \param level
        */
        virtual void operator () (int* tag, const int* tlo, const int* thi,
                                  const int* tagval, const int* clearval,
                                  Real* data, const int* data_lo, const int* data_hi,
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const;
    protected:

        ErrorFuncDefault   m_func;
        ErrorFunc3DDefault m_func3D;
    };
    class ErrorFunc2
    {
    public:

        /**
        * \brief Bogus constructor.
        */
        ErrorFunc2 ();

        /**
        * \brief A Constructor.
        *
        * \param inFunc
        */
        // No explicit because some codes rely on it.
        ErrorFunc2 (ErrorFunc2Default inFunc);

        /**
        * \brief Return a ptr to a clone of this object.
        * It is the responsibility of the caller to delete the result.
        */
        virtual ErrorFunc2* clone () const;

        /**
        * \brief Destructor.
        */
        virtual ~ErrorFunc2 ();

        /**
        * \brief Tag cells cells using "v2" interface
        *
        * \param tag
        * \param tlo
        * \param thi
        * \param tagval
        * \param clearval
        * \param data
        * \param data_lo
        * \param data_hi
        * \param lo
        * \param hi
        * \param nvar
        * \param domain_lo
        * \param domain_hi
        * \param dx
        * \param level
        * \param avg
        */
        virtual void operator () (int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                  const int* tagval, const int* clearval,
                                  Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const int* level, const Real* avg) const;
    protected:

        ErrorFunc2Default  m_func;
    };

    /**
    * \brief The constructors.
    *
    * \param nm
    * \param ng
    * \param etyp
    * \param f2
    */
    ErrorRec (const std::string&          nm,
              int                         ng,
              ErrorType                   etyp,
              const ErrorRec::ErrorFunc2& f2);

    ErrorRec (const std::string&         nm,
              int                        ng,
              ErrorType                  etyp,
              const ErrorRec::ErrorFunc& f);
    virtual ~ErrorRec();

    /**
    * \brief The name of the quantity to derive.
    */
    const std::string& name () const noexcept;

    /**
    * \brief The number of extra zones needed for derivation.
    */
    int nGrow () const noexcept;

    /**
    * \brief The type of the error tagging.
    */
    ErrorType errType () const noexcept;

    /**
    * \brief The extern "C" functions to do the error tagging.
    */
    virtual const ErrorRec::ErrorFunc&  errFunc () const;
    virtual const ErrorRec::ErrorFunc2& errFunc2() const;

private:

    //! Name of quantity to derive.
    std::string derive_name;

    //! Number of extra zones.
    int ngrow;

    //! The type of Error.
    ErrorType err_type;

    //! Functions to do error estimation.
    ErrorFunc*  err_func;
    ErrorFunc2* err_func2;
};


/**
* \brief A List of ErrorRecs.
*
* Container class for ErrorRecs.
*/

class ErrorList
{
public:
    ErrorList() noexcept {}

    /**
    * \brief The number of ErrorRecs in the list.
    */
    int size () const noexcept;

    /**
    * \brief Append a new ErrorRec to the list.
    *
    * \param name
    * \param nextra
    * \param typ
    * \param func
    */
    void add (const std::string&  name,
              int                 nextra,
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc&    func);

    void add (const std::string&  name,
              int                 nextra,
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc2&   func);

    //! The kth ErrorRec.
    const ErrorRec& operator[] (int k) const noexcept;

    void clear (bool rs0 = false) { vec.clear(); if(rs0) { vec.resize(0); } }

private:

    Vector<std::unique_ptr<ErrorRec> > vec;
};

std::ostream& operator << (std::ostream& os, const ErrorList& elst);

  struct AMRErrorTagInfo
  {
    int m_max_level = 100000;
    Real m_min_time = std::numeric_limits<Real>::lowest();
    Real m_max_time = std::numeric_limits<Real>::max();
    RealBox m_realbox;

    AMRErrorTagInfo& SetMaxLevel (int max_level) noexcept {
      m_max_level = max_level;
      return *this;
    }
    AMRErrorTagInfo& SetMinTime (amrex::Real min_time) noexcept {
      m_min_time = min_time;
      return *this;
    }
    AMRErrorTagInfo& SetMaxTime (amrex::Real max_time) noexcept {
      m_max_time = max_time;
      return *this;
    }
    AMRErrorTagInfo& SetRealBox (const amrex::RealBox& realbox) noexcept {
      m_realbox = realbox;
      return *this;
    }
  };

  class AMRErrorTag
  {
  public:

    enum TEST {GRAD=0, LESS, GREATER, VORT, BOX, USER};

    struct UserFunc
    {
      virtual void operator() (const amrex::Box&                       bx,
                               amrex::Array4<const amrex::Real> const& dat,
                               amrex::Array4<char> const&              tag,
                               amrex::Real                             time,
                               int                                     level,
                               char                                    tagval,
                               char                                    clearval) = 0;
    };

    explicit AMRErrorTag (const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_test(BOX), m_field(std::string()), m_info(info) {m_ngrow = SetNGrow();}

    AMRErrorTag (amrex::Real            value,
                 AMRErrorTag::TEST      test,
                 const std::string&     field,
                 const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_value(value), m_test(test), m_field(field), m_info(info) {m_ngrow = SetNGrow();}

    AMRErrorTag (AMRErrorTag::UserFunc* userfunc,
                 const std::string&     field,
                 int                    ngrow,
                 const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_userfunc(userfunc), m_field(field), m_info(info), m_ngrow(ngrow) {}

    virtual void operator() (amrex::TagBoxArray&    tb,
                             const amrex::MultiFab* mf,
                             char                   clearval,
                             char                   tagval,
                             amrex::Real            time,
                             int                    level,
                             const amrex::Geometry& geom) const noexcept;

    int NGrow() const noexcept {return m_ngrow;}
    const std::string& Field () const noexcept {return m_field;}

  protected:
    int SetNGrow () const noexcept;

    Real m_value;    
    TEST m_test;
    UserFunc* m_userfunc = nullptr;
    std::string m_field;
    AMRErrorTagInfo m_info;
    int m_ngrow;
  };
}

#endif
