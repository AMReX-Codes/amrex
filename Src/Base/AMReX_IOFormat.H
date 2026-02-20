#ifndef AMREX_IO_FORMAT_H_
#define AMREX_IO_FORMAT_H_

/**
 * \file AMReX_IOFormat.H
 * \brief RAII saver/restorer for \c std::basic_ios format state.
 */

#include <ios>

namespace amrex {

/**
 * \class amrex::IOFormatSaver
 * \ingroup amrex_utilities
 * \brief RAII saver and restorer for \c std::basic_ios format state.
 *
 * Saves the fill character, format flags, precision, and width of the
 * stream on construction and restores them on destruction.
 */
template <class CharT, class Traits>
class IOFormatSaver
{
public:
    using BasicIos = std::basic_ios<CharT, Traits>;

    //! Save the format state of \p ios.
    explicit IOFormatSaver (BasicIos& ios)
        : m_ios(&ios),
          m_fill(ios.fill()),
          m_flags(ios.flags()),
          m_precision(ios.precision()),
          m_width(ios.width())
    {}

    ~IOFormatSaver ()
    {
        m_ios->fill(m_fill);
        m_ios->flags(m_flags);
        m_ios->precision(m_precision);
        m_ios->width(m_width);
    }

    IOFormatSaver (IOFormatSaver const&) = delete;
    IOFormatSaver (IOFormatSaver &&) noexcept = delete;
    IOFormatSaver& operator= (IOFormatSaver const&) = delete;
    IOFormatSaver& operator= (IOFormatSaver &&) noexcept = delete;

private:
    BasicIos*                   m_ios;
    CharT                       m_fill;
    typename BasicIos::fmtflags m_flags;
    std::streamsize             m_precision;
    std::streamsize             m_width;
};

}

#endif
