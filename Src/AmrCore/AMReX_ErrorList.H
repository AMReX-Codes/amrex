#ifndef AMREX_ErrorList_H_
#define AMREX_ErrorList_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>
#include <AMReX_Vector.H>

#include <string>
#include <memory>

namespace amrex {


extern "C"
{

    /**
    * \brief Type of extern "C" function called by ErrorRec to do tagging of cells for refinement.
    */
    using ErrorFuncDefault = void (*)(int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                      const int* tagval, const int* clearval,
                                      amrex::Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                      const int* lo, const int * hi, const int* nvar,
                                      const int* domain_lo, const int* domain_hi,
                                      const amrex::Real* dx, const amrex::Real* xlo,
                                      const amrex::Real* prob_lo, const amrex::Real* time,
                                      const int* level);

    using ErrorFunc2Default = void (*)(int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                       const int* tagval, const int* clearval,
                                       amrex::Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                       const int* lo, const int * hi, const int* nvar,
                                       const int* domain_lo, const int* domain_hi,
                                       const amrex::Real* dx, const int* level, const amrex::Real* avg);


    /**
    * \brief Dimension agnostic version that always has three elements.
    * Note that this is only implemented for the ErrorFunc class, not ErrorFunc2.
    *
    * \param tag      Cell-status array to update (tagged or clear).
    * \param tlo      Lower bounds of the tag array.
    * \param thi      Upper bounds of the tag array.
    * \param tagval   Value that indicates a tagged cell.
    * \param clearval Value that indicates an untagged cell.
    * \param data     Pointer to data.
    * \param data_lo  Lower bounds of the data array.
    * \param data_hi  Upper bounds of the data array.
    * \param lo       Lower bounds of the tile being tagged.
    * \param hi       Upper bounds of the tile being tagged.
    * \param nvar     Number of data components supplied.
    * \param domain_lo Lower bounds of the problem domain.
    * \param domain_hi Upper bounds of the problem domain.
    * \param dx       Cell spacing vector.
    * \param xlo      Physical location of the lower corner of the tile.
    * \param prob_lo  Physical location of the lower domain corner.
    * \param time     Current simulation time.
    * \param level    AMR level index being tagged.
    */
    using ErrorFunc3DDefault = void (*)(int* tag, const int* tlo, const int* thi,
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
        [[nodiscard]] virtual ErrorFunc* clone () const;

        /**
        * \brief Destructor.
        */
        virtual ~ErrorFunc () = default;

        ErrorFunc (ErrorFunc const&) = default;
        ErrorFunc (ErrorFunc &&) = delete;
        ErrorFunc& operator= (ErrorFunc const&) = default;
        ErrorFunc& operator= (ErrorFunc &&) = delete;

        /**
        * \brief Tag cells using "regular" function.
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
        * \param tag      Cell-status array to update (tagged or clear).
        * \param tlo      Lower bounds of the tag array.
        * \param thi      Upper bounds of the tag array.
        * \param tagval   Value that indicates a tagged cell.
        * \param clearval Value that indicates an untagged cell.
        * \param data     Pointer to data.
        * \param dlo      Lower bounds of the data array.
        * \param dhi      Upper bounds of the data array.
        * \param lo       Lower bounds of the tile being tagged.
        * \param hi       Upper bounds of the tile being tagged.
        * \param nvar     Number of data components supplied.
        * \param domain_lo Lower bounds of the problem domain.
        * \param domain_hi Upper bounds of the problem domain.
        * \param dx       Cell spacing vector.
        * \param xlo      Physical location of the lower corner of the tile.
        * \param prob_lo  Physical location of the lower domain corner.
        * \param time     Current simulation time.
        * \param level    AMR level index being tagged.
        */
        virtual void operator () (int* tag, const int* tlo, const int* thi,
                                  const int* tagval, const int* clearval,
                                  Real* data, const int* dlo, const int* dhi,
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const;
    protected:

        ErrorFuncDefault   m_func{nullptr};
        ErrorFunc3DDefault m_func3D{nullptr};
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
        [[nodiscard]] virtual ErrorFunc2* clone () const;

        /**
        * \brief Destructor.
        */
        virtual ~ErrorFunc2 () = default;

        ErrorFunc2 (ErrorFunc2 const&) = default;
        ErrorFunc2 (ErrorFunc2 &&) = delete;
        ErrorFunc2& operator= (ErrorFunc2 const&) = default;
        ErrorFunc2& operator= (ErrorFunc2 &&) = delete;

        /// \cond DOXYGEN_IGNORE
        /**
        * \brief Tag cells using the alternative interface that supports average-based tagging.
        *
        * \param tag      Cell-status array to update (tagged or clear).
        * \param tlo      Lower bounds of the tag array.
        * \param thi      Upper bounds of the tag array.
        * \param tagval   Value indicating a tagged cell.
        * \param clearval Value indicating an untagged cell.
        * \param data     Pointer to data.
        * \param data_lo  Lower bounds of the data array.
        * \param data_hi  Upper bounds of the data array.
        * \param lo       Lower bounds of the tile being tagged.
        * \param hi       Upper bounds of the tile being tagged.
        * \param nvar     Number of state components supplied.
        * \param domain_lo Lower bounds of the problem domain.
        * \param domain_hi Upper bounds of the problem domain.
        * \param dx       Cell spacing vector.
        * \param level    AMR level index.
        * \param avg      Optional coarse-level averages for comparison.
        */
        virtual void operator () (int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
                                  const int* tagval, const int* clearval,
                                  Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const int* level, const Real* avg) const;
        /// \endcond
    protected:

        ErrorFunc2Default  m_func{nullptr};
    };

    /**
    * \brief Construct an error record that uses the ErrorFunc2 interface.
    *
    * \param nm   Derived quantity name used for tagging.
    * \param ng   Number of grow cells required by the tagging functor.
    * \param etyp Error tagging type (Special/Standard/UseAverage).
    * \param f2   Functor describing how to tag cells.
    */
    ErrorRec (std::string nm, int ng, ErrorType etyp,
              const ErrorRec::ErrorFunc2& f2);

    //! Construct an error record that uses the ErrorFunc interface (same parameters as above).
    ErrorRec (std::string nm, int ng, ErrorType etyp,
              const ErrorRec::ErrorFunc& f);

    virtual ~ErrorRec ();

    ErrorRec (ErrorRec const&) = delete;
    ErrorRec (ErrorRec &&) = delete;
    ErrorRec& operator= (ErrorRec const&) = delete;
    ErrorRec& operator= (ErrorRec &&) = delete;

    /**
    * \brief The name of the quantity to derive.
    */
    [[nodiscard]] const std::string& name () const noexcept;

    /**
    * \brief The number of extra zones needed for derivation.
    */
    [[nodiscard]] int nGrow () const noexcept;

    /**
    * \brief The type of the error tagging.
    */
    [[nodiscard]] ErrorType errType () const noexcept;

    /**
    * \brief The extern "C" functions to do the error tagging.
    */
    [[nodiscard]] virtual const ErrorRec::ErrorFunc&  errFunc () const;
    [[nodiscard]] virtual const ErrorRec::ErrorFunc2& errFunc2() const;

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
    ErrorList() noexcept = default;

    /**
    * \brief The number of ErrorRecs in the list.
    */
    [[nodiscard]] int size () const noexcept;

    /**
    * \brief Append a new ErrorRec to the list.
    *
    * \param name   Derived quantity name used for tagging.
    * \param nextra Number of grow cells required by the functor.
    * \param typ    Tagging type (Special/Standard/UseAverage).
    * \param func   ErrorFunc functor describing tagging logic.
    */
    void add (const std::string&  name,
              int                 nextra,
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc&    func);

    //! Overload that accepts ErrorFunc2 instead of ErrorFunc.
    void add (const std::string&  name,
              int                 nextra,
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc2&   func);

    //! The kth ErrorRec.
    [[nodiscard]] const ErrorRec& operator[] (int k) const noexcept;

    void clear (bool rs0 = false) { vec.clear(); if(rs0) { vec.resize(0); } }

private:

    Vector<std::unique_ptr<ErrorRec> > vec;
};

std::ostream& operator << (std::ostream& os, const ErrorList& elst);

  //! Optional knobs that control AMRErrorTag behavior (time window, level limits, etc.).
  struct AMRErrorTagInfo
  {
    int m_max_level = 1000;
    Real m_min_time = std::numeric_limits<Real>::lowest();
    Real m_max_time = std::numeric_limits<Real>::max();
    int m_volume_weighting = 0;
    int m_derefine = 0;
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
    AMRErrorTagInfo& SetVolumeWeighting (int volume_weighting) noexcept {
      m_volume_weighting = volume_weighting;
      return *this;
    }
    AMRErrorTagInfo& SetDerefine (int derefine) noexcept {
      m_derefine = derefine;
      return *this;
    }
  };

  /**
   * \brief Flexible tagging rule that supports gradients, thresholds, boxes, and user functors.
   */
  class AMRErrorTag
  {
  public:

    enum TEST {GRAD=0, RELGRAD, LESS, GREATER, VORT, BOX, USER, PARSER};

    struct UserFunc
    {
      virtual ~UserFunc () = default;

      UserFunc (UserFunc const&) = default;
      UserFunc (UserFunc &&) = default;
      UserFunc& operator= (UserFunc const&) = default;
      UserFunc& operator= (UserFunc &&) = default;

      virtual void operator() (const amrex::Box&                       bx,
                               amrex::Array4<const amrex::Real> const& dat,
                               amrex::Array4<char> const&              tag,
                               amrex::Real                             time,
                               int                                     level,
                               char                                    tagval,
                               char                                    clearval) = 0;
    };

    //! Construct a tag with default BOX test and the supplied metadata.
    explicit AMRErrorTag (const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_info(info), m_ngrow(SetNGrow()) {}

    /**
     * \brief Construct a tag using a scalar threshold across all levels.
     *
     * \param value Threshold (or reference) value.
     * \param test  Test mode (e.g., GREATER, LESS).
     * \param field MultiFab component name to examine.
     * \param info  Optional metadata describing min/max time, level limits, etc.
     */
    AMRErrorTag (amrex::Real            value,
                 AMRErrorTag::TEST      test,
                 std::string            field,
                 const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_value(info.m_max_level, value), m_test(test), m_field(std::move(field)), m_info(info),
        m_ngrow(SetNGrow())
      {
      }

    /**
     * \brief Construct a tag with per-level threshold values.
     *
     * \param value Vector of thresholds (last entry reused for higher levels if needed).
     * \param test  Test mode.
     * \param field MultiFab component name to examine.
     * \param info  Optional metadata (time window, level window, etc.).
     */
    AMRErrorTag (amrex::Vector<amrex::Real>  value,
                 AMRErrorTag::TEST           test,
                 std::string                 field,
                 const AMRErrorTagInfo&      info = AMRErrorTagInfo()) noexcept
      : m_test(test), m_field(std::move(field)), m_info(info), m_ngrow(SetNGrow())
      {
          AMREX_ASSERT(!value.empty());
          m_value.resize(info.m_max_level);
          for (int i = 0; i < m_value.size() && i < value.size(); ++i) {
              m_value[i] = value[i];
          }
          // If the user didn't provided a value for every level,
          // assume the last value holds for all higher levels.
          for (auto i = int(value.size()); i < m_value.size(); ++i) {
              m_value[i] = value[value.size()-1];
          }
      }

    /**
     * \brief Construct a tag that delegates to a user functor.
     *
     * \param userfunc User-defined functor (ownership retained by caller).
     * \param field    MultiFab component name to examine.
     * \param ngrow    Number of grow cells the functor requires.
     * \param info     Optional metadata.
     */
    AMRErrorTag (AMRErrorTag::UserFunc* userfunc,
                 std::string            field,
                 int                    ngrow,
                 const AMRErrorTagInfo& info = AMRErrorTagInfo()) noexcept
      : m_userfunc(userfunc), m_field(std::move(field)), m_info(info), m_ngrow(ngrow) {}

    //! Construct AMRErrorTag from a Parser that returns Real for given
    //! (x,y,z,t). This class's operator() can then be used to tag regions
    //! where the Parser's return value is greater than 0.
    explicit AMRErrorTag (Parser parser,
                          const AMRErrorTagInfo& info = AMRErrorTagInfo());

    /**
     * \brief Apply the tagging rule to a TagBoxArray/MultiFab pair.
     *
     * \param tb      Destination TagBoxArray to update.
     * \param mf      Optional pointer to the MultiFab whose data drive the tagging (can be nullptr for USER mode).
     * \param clearval Value indicating an untagged cell.
     * \param tagval   Value indicating a tagged cell.
     * \param time     Simulation time.
     * \param level    AMR level index.
     * \param geom     Geometry describing the level (for coordinate-dependent tests).
     */
    void operator() (amrex::TagBoxArray&    tb,
                     const amrex::MultiFab* mf,
                     char                   clearval,
                     char                   tagval,
                     amrex::Real            time,
                     int                    level,
                     const amrex::Geometry& geom) const noexcept;

    [[nodiscard]] int NGrow() const noexcept {return m_ngrow;}
    [[nodiscard]] const std::string& Field () const noexcept {return m_field;}

    [[nodiscard]] AMRErrorTagInfo& GetInfo () noexcept {return m_info;}
    [[nodiscard]] AMRErrorTagInfo const& GetInfo () const noexcept {return m_info;}
    void SetInfo (AMRErrorTagInfo const& info) noexcept {m_info = info;}

  protected:
    [[nodiscard]] int SetNGrow () const noexcept;

    Vector<Real> m_value;
    TEST m_test{BOX};
    UserFunc* m_userfunc = nullptr;
    std::unique_ptr<Parser> m_parser;
    ParserExecutor<4> m_parser_exe; // 4: (x,y,z,t)
    std::string m_field;
    AMRErrorTagInfo m_info;
    int m_ngrow = 0;
  };
}

#endif
