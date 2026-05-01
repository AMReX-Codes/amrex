#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstring>

#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"
#include "AMReX_Vector.H"
#include "AMReX_VectorIO.H"
#include "AMReX_Utility.H"

using namespace amrex;

namespace {

[[noreturn]] void ThrowingErrorHandler (const char* msg)
{
    throw std::runtime_error(msg ? msg : "amrex::Error");
}

class ErrorHandlerScope
{
public:
    explicit ErrorHandlerScope (amrex::ErrorHandler handler) noexcept
        : m_old_handler(amrex::system::error_handler)
    {
        amrex::SetErrorHandler(handler);
    }

    ~ErrorHandlerScope () noexcept
    {
        amrex::SetErrorHandler(m_old_handler);
    }

private:
    amrex::ErrorHandler m_old_handler;
};

}

void testIntIO(const IntDescriptor& id_out) {

    std::string data_file_name   = "int_data.dat";
    std::string header_file_name = "int_header_H";

    amrex::Vector<int> idata_out;
    for (int i = -99; i <= 100; ++i) {
        idata_out.push_back(i);
    }

    std::ofstream ofs;
    ofs.open(data_file_name.c_str(), std::ios::out|std::ios::binary);
    writeIntData(idata_out.data(), idata_out.size(), ofs, id_out);
    ofs.close();

    ofs.open(header_file_name.c_str(), std::ios::out);
    ofs << id_out << "\n";
    ofs.close();

    IntDescriptor id_in;
    std::ifstream ifs;
    ifs.open(header_file_name.c_str(), std::ios::in);
    ifs >> id_in;
    ifs.close();

    AMREX_ALWAYS_ASSERT(id_out == id_in);

    amrex::Vector<int> idata_in(idata_out.size());
    ifs.open(data_file_name.c_str(), std::ios::in|std::ios::binary);
    readIntData(idata_in.data(), idata_in.size(), ifs, id_in);
    ifs.close();

    for (int i = 0; i < static_cast<int>(idata_in.size()); ++i) {
        AMREX_ALWAYS_ASSERT(idata_in[i] == idata_out[i]);
    }
}

void testLongIO(const IntDescriptor& id_out) {

    std::string data_file_name   = "long_data.dat";
    std::string header_file_name = "long_header_H";

    amrex::Vector<Long> idata_out;
    for (int i = -99; i <= 100; ++i) {
        idata_out.push_back(i);
    }

    std::ofstream ofs;
    ofs.open(data_file_name.c_str(), std::ios::out|std::ios::binary);
    writeLongData(idata_out.data(), idata_out.size(), ofs, id_out);
    ofs.close();

    ofs.open(header_file_name.c_str(), std::ios::out);
    ofs << id_out << "\n";
    ofs.close();

    IntDescriptor id_in;
    std::ifstream ifs;
    ifs.open(header_file_name.c_str(), std::ios::in);
    ifs >> id_in;
    ifs.close();

    AMREX_ALWAYS_ASSERT(id_out == id_in);
    amrex::Vector<Long> idata_in(idata_out.size());

    ifs.open(data_file_name.c_str(), std::ios::in|std::ios::binary);
    readLongData(idata_in.data(), idata_in.size(), ifs, id_in);
    ifs.close();

    for (int i = 0; i < static_cast<int>(idata_in.size()); ++i) {
        AMREX_ALWAYS_ASSERT(idata_in[i] == idata_out[i]);
    }
}

void testTruncatedIntegerIO () {

    ErrorHandlerScope scope(ThrowingErrorHandler);

    bool native_int_failed = false;
    try {
        amrex::Vector<int> idata_in(2);
        std::istringstream is(std::string(sizeof(int), '\0'));
        readIntData(idata_in.data(), idata_in.size(), is, FPC::NativeIntDescriptor());
    } catch (std::runtime_error const&) {
        native_int_failed = true;
    }
    AMREX_ALWAYS_ASSERT(native_int_failed);

    bool native_long_failed = false;
    try {
        amrex::Vector<Long> idata_in(2);
        std::istringstream is(std::string(sizeof(Long), '\0'));
        readLongData(idata_in.data(), idata_in.size(), is, FPC::NativeLongDescriptor());
    } catch (std::runtime_error const&) {
        native_long_failed = true;
    }
    AMREX_ALWAYS_ASSERT(native_long_failed);

    bool converted_int_failed = false;
    try {
        auto const native_order = FPC::NativeIntDescriptor().order();
        auto const converted_order = (native_order == IntDescriptor::NormalOrder)
            ? IntDescriptor::ReverseOrder : IntDescriptor::NormalOrder;
        IntDescriptor id_in(FPC::NativeIntDescriptor().numBytes(), converted_order);
        amrex::Vector<int> idata_in(2);
        std::istringstream is(std::string(static_cast<std::size_t>(id_in.numBytes()), '\0'));
        readIntData(idata_in.data(), idata_in.size(), is, id_in);
    } catch (std::runtime_error const&) {
        converted_int_failed = true;
    }
    AMREX_ALWAYS_ASSERT(converted_int_failed);
}

void testRealIO(const RealDescriptor& rd_out) {

    std::string data_file_name   = "real_data.dat";
    std::string header_file_name = "real_header_H";

    amrex::Vector<Real> rdata_out;
    for (int i = -99; i <= 100; ++i) {
        rdata_out.push_back(amrex::Random());
    }

    std::ofstream ofs;
    ofs.open(data_file_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(rdata_out.data(), rdata_out.size(), ofs, rd_out);
    ofs.close();

    ofs.open(header_file_name.c_str(), std::ios::out);
    ofs << rd_out << "\n";
    ofs.close();

    RealDescriptor rd_in;
    std::ifstream ifs;
    ifs.open(header_file_name.c_str(), std::ios::in);
    ifs >> rd_in;
    ifs.close();

    AMREX_ALWAYS_ASSERT(rd_out == rd_in);

    amrex::Vector<Real> rdata_in(rdata_out.size());
    ifs.open(data_file_name.c_str(), std::ios::in|std::ios::binary);
    readRealData(rdata_in.data(), rdata_in.size(), ifs, rd_in);
    ifs.close();

    for (int i = 0; i < static_cast<int>(rdata_in.size()); ++i) {
        if (rd_in == FPC::Native32RealDescriptor() ||
            rd_in == FPC::Ieee32NormalRealDescriptor()) {
            AMREX_ALWAYS_ASSERT(std::abs(rdata_in[i] - rdata_out[i]) <= 1e-7);
        } else{
            AMREX_ALWAYS_ASSERT(rdata_in[i] == rdata_out[i]);
        }
    }
}

void testFloatIO(const RealDescriptor& rd_out) {

    std::string data_file_name   = "float_data.dat";
    std::string header_file_name = "float_header_H";

    amrex::Vector<float> rdata_out;
    for (int i = -99; i <= 100; ++i) {
        rdata_out.push_back(static_cast<float>(amrex::Random()));
    }

    std::ofstream ofs;
    ofs.open(data_file_name.c_str(), std::ios::out|std::ios::binary);
    writeFloatData(rdata_out.data(), rdata_out.size(), ofs, rd_out);
    ofs.close();

    ofs.open(header_file_name.c_str(), std::ios::out);
    ofs << rd_out << "\n";
    ofs.close();

    RealDescriptor rd_in;
    std::ifstream ifs;
    ifs.open(header_file_name.c_str(), std::ios::in);
    ifs >> rd_in;
    ifs.close();

    AMREX_ALWAYS_ASSERT(rd_out == rd_in);

    amrex::Vector<float> rdata_in(rdata_out.size());
    ifs.open(data_file_name.c_str(), std::ios::in|std::ios::binary);
    readFloatData(rdata_in.data(), rdata_in.size(), ifs, rd_in);
    ifs.close();

    for (int i = 0; i < static_cast<int>(rdata_in.size()); ++i) {
        if (rd_in == FPC::Native32RealDescriptor() ||
            rd_in == FPC::Ieee32NormalRealDescriptor()) {
            AMREX_ALWAYS_ASSERT(std::abs(rdata_in[i] - rdata_out[i]) <= 1e-7);
        } else{
            AMREX_ALWAYS_ASSERT(rdata_in[i] == rdata_out[i]);
        }
    }
}

void testDoubleIO(const RealDescriptor& rd_out) {

    std::string data_file_name   = "double_data.dat";
    std::string header_file_name = "double_header_H";

    amrex::Vector<double> rdata_out;
    for (int i = -99; i <= 100; ++i) {
        rdata_out.push_back(amrex::Random());
    }

    std::ofstream ofs;
    ofs.open(data_file_name.c_str(), std::ios::out|std::ios::binary);
    writeDoubleData(rdata_out.data(), rdata_out.size(), ofs, rd_out);
    ofs.close();

    ofs.open(header_file_name.c_str(), std::ios::out);
    ofs << rd_out << "\n";
    ofs.close();

    RealDescriptor rd_in;
    std::ifstream ifs;
    ifs.open(header_file_name.c_str(), std::ios::in);
    ifs >> rd_in;
    ifs.close();

    AMREX_ALWAYS_ASSERT(rd_out == rd_in);

    amrex::Vector<double> rdata_in(rdata_out.size());
    ifs.open(data_file_name.c_str(), std::ios::in|std::ios::binary);
    readDoubleData(rdata_in.data(), rdata_in.size(), ifs, rd_in);
    ifs.close();

    for (int i = 0; i < static_cast<int>(rdata_in.size()); ++i) {
        if (rd_in == FPC::Native32RealDescriptor() ||
            rd_in == FPC::Ieee32NormalRealDescriptor()) {
            AMREX_ALWAYS_ASSERT(std::abs(rdata_in[i] - rdata_out[i]) <= 1e-7);
        } else{
            AMREX_ALWAYS_ASSERT(rdata_in[i] == rdata_out[i]);
        }
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    IntDescriptor little16(2, IntDescriptor::ReverseOrder);
    IntDescriptor little32(4, IntDescriptor::ReverseOrder);
    IntDescriptor little64(8, IntDescriptor::ReverseOrder);

    IntDescriptor big16(2, IntDescriptor::NormalOrder);
    IntDescriptor big32(4, IntDescriptor::NormalOrder);
    IntDescriptor big64(8, IntDescriptor::NormalOrder);

    testIntIO(little16);
    testIntIO(little32);
    testIntIO(little64);

    testIntIO(big16);
    testIntIO(big32);
    testIntIO(big64);

    testLongIO(little16);
    testLongIO(little32);
    testLongIO(little64);

    testLongIO(big16);
    testLongIO(big32);
    testLongIO(big64);

    testTruncatedIntegerIO();

    testRealIO(FPC::NativeRealDescriptor());
    testRealIO(FPC::Native32RealDescriptor());
    testRealIO(FPC::Ieee32NormalRealDescriptor());
    testRealIO(FPC::Ieee64NormalRealDescriptor());

    testFloatIO(FPC::Native32RealDescriptor());
    testFloatIO(FPC::NativeRealDescriptor());
    testFloatIO(FPC::Ieee32NormalRealDescriptor());
    testFloatIO(FPC::Ieee64NormalRealDescriptor());

    testDoubleIO(FPC::Native32RealDescriptor());
    testDoubleIO(FPC::NativeRealDescriptor());
    testDoubleIO(FPC::Ieee32NormalRealDescriptor());
    testDoubleIO(FPC::Ieee64NormalRealDescriptor());

    amrex::Print() << "passed!" << '\n';

    amrex::Finalize();
}
