#ifndef AMREX_ANY_H_
#define AMREX_ANY_H_
#include <AMReX_Config.H>

#include <memory>
#include <typeinfo>
#include <type_traits>

namespace amrex {

/**
 * This class is similar to std::any.  However, amrex::Any can store
 * move-only types (e.g., amrex::MultiFab), whereas std::any is for copy
 * constructible types only.  An object of this class is non-copyable.
 */
class Any
{
public:

    Any () = default;
    ~Any () = default;

    Any (Any const& rhs) = delete;
    Any& operator= (Any const& rhs) = delete;

    Any (Any && rhs) = default;
    Any& operator= (Any && rhs) = default;

    //! Constructs by moving the given object.
    template <typename MF,
              std::enable_if_t<!std::is_same_v<Any,std::decay_t<MF>>,int> = 0>
    Any (MF && mf) // NOLINT(bugprone-forwarding-reference-overload)
        : m_ptr(std::make_unique<innards<MF> >(std::forward<MF>(mf)))
        {}

    //! Assigns by moving the given object.
    template <typename MF,
              std::enable_if_t<!std::is_same_v<Any,std::decay_t<MF>>,int> = 0>
    Any& operator= (MF && mf) {
        m_ptr = std::make_unique<innards<MF> >(std::forward<MF>(mf));
        return *this;
    }

    //! Returns the contained type.
    [[nodiscard]] const std::type_info& Type () const {
        if (m_ptr) {
            return m_ptr->Type();
        } else {
            return typeid(void);
        }
    }

    //! Returns a reference to the contained object.
    template <typename MF>
    [[nodiscard]] MF& get () {
        if (auto p0 = dynamic_cast<innards<MF>*>(m_ptr.get())) {
            return p0->m_mf;
        } else {
            return dynamic_cast<innards<MF&>&>(*m_ptr).m_mf;
        }
    }

    //! Returns a const reference to the contained object.
    template <typename MF>
    [[nodiscard]] MF const& get () const {
        if (auto p0 = dynamic_cast<innards<MF>*>(m_ptr.get())) {
            return p0->m_mf;
        } else if (auto p1 = dynamic_cast<innards<MF&>*>(m_ptr.get())) {
            return p1->m_mf;
        } else {
            return dynamic_cast<innards<MF const&> const&>(*m_ptr).m_mf;
        }
    }

    template <typename MF>
    [[nodiscard]] bool is () const { return Type() == typeid(MF); }

    [[nodiscard]] bool hasValue () const { return m_ptr != nullptr; }

private:
    struct innards_base // NOLINT(cppcoreguidelines-special-member-functions)
    {
        [[nodiscard]] virtual const std::type_info& Type () const = 0;
        virtual ~innards_base () = default;
    };

    template <typename MF>
    struct innards final : innards_base // NOLINT(cppcoreguidelines-special-member-functions)
    {
        innards (MF && mf)
            : m_mf(std::move(mf))
            {}

        ~innards () final = default;

        [[nodiscard]] const std::type_info& Type () const final {
            return typeid(MF);
        }

        MF m_mf;
    };

    std::unique_ptr<innards_base> m_ptr;
};

}

#endif
