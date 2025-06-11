#ifndef AMREX_IO_FORMAT_H_
#define AMREX_IO_FORMAT_H_

#include <ios>

namespace amrex {

/*
 * \brief I/O stream format saver
 *
 * This class can be used to save and restore I/O stream format in a RAII
 * way. It handles fill, fmtflag, precision, and width.
 */
template <class CharT, class Traits>
class IOFormatSaver
{
public:
    using BasicIos = std::basic_ios<CharT, Traits>;

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
