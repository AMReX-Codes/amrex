#ifndef AMREX_ANY_H_
#define AMREX_ANY_H_
#include <AMReX_Config.H>

#include <memory>
#include <typeinfo>

namespace amrex {

/**
 * This class is similar to std::any.  However, amrex::Any is for storing
 * move-only types (e.g., amrex::MultiFab), whereas std::any is for copy
 * constructible types.  An object of this class is non-copyable.
 */
class Any
{
public:

    //! Default constructor.  By default it stores an int.
    Any ()
        : m_ptr(std::make_unique<innards<int> >(0))
        {}

    ~Any () = default;

    Any (Any const& rhs) = delete;
    Any& operator= (Any const& rhs) = delete;

    Any (Any && rhs) = default;
    Any& operator= (Any && rhs) = default;

    //! Constructs by moving the given object.
    template <typename MF>
    Any (MF && mf)
        : m_ptr(std::make_unique<innards<MF> >(std::forward<MF>(mf)))
        {}

    //! Assigns by moving the given object.
    template <typename MF>
    void operator= (MF && mf) {
        m_ptr = std::make_unique<innards<MF> >(std::forward<MF>(mf));
    }

    //! Returns the contained type.
    const std::type_info& Type () const {
        return m_ptr->Type();
    }

    //! Returns a reference to the contained object.
    template <typename MF>
    MF& get () { return dynamic_cast<innards<MF>&>(*m_ptr).m_mf; }

    //! Returns a const reference to the contained object.
    template <typename MF>
    MF const& get () const { return dynamic_cast<innards<MF> const&>(*m_ptr).m_mf; }

    template <typename MF>
    bool is () const { return m_ptr->Type() == typeid(MF); }

private:
    struct innards_base {
        virtual const std::type_info& Type () const = 0;
    };

    template <typename MF>
    struct innards : innards_base
    {
        innards(MF && mf)
            : m_mf(std::forward<MF>(mf))
            {}

        virtual const std::type_info& Type () const override {
            return typeid(MF);
        }

        MF m_mf;
    };

    std::unique_ptr<innards_base> m_ptr;
};

}

#endif
