#include "AMReX_VectorIO.H"

using namespace amrex;

void amrex::writeIntData(const int* data, std::size_t size, std::ostream& os,
                         const IntDescriptor& id)
{
    if (id == FPC::NativeIntDescriptor())
    {
        os.write((char*) data, size*sizeof(int));
    }
    else if (id.numBytes() == 2)
    {
        writeIntData<std::int16_t, int>(data, size, os, id);
    }
    else if (id.numBytes() == 4)
    {
        writeIntData<std::int32_t, int>(data, size, os, id);
    }
    else if (id.numBytes() == 8)
    {
        writeIntData<std::int64_t, int>(data, size, os, id);
    }
    else {
        amrex::Error("Don't know how to work with this integer type.");
    }
}

void amrex::readIntData(int* data, std::size_t size, std::istream& is,
                        const IntDescriptor& id)
{
    if (id == FPC::NativeIntDescriptor()) 
    {
        is.read((char*) data, size * id.numBytes());
    } 
    else if (id.numBytes() == 2)
    {
        readIntData<int, std::int16_t>(data, size, is, id);
    }
    else if (id.numBytes() == 4)
    {
        readIntData<int, std::int32_t>(data, size, is, id);
    }
    else if (id.numBytes() == 8)
    {
        readIntData<int, std::int64_t>(data, size, is, id);
    }
    else {
        amrex::Error("Don't know how to work with this integer type.");
    }
}

void amrex::writeLongData(const Long* data, std::size_t size, std::ostream& os,
                          const IntDescriptor& id)
{
    if (id == FPC::NativeLongDescriptor())
    {
        os.write((char*) data, size*sizeof(Long));
    }
    else if (id.numBytes() == 2)
    {
        writeIntData<std::int16_t, Long>(data, size, os, id);
    }
    else if (id.numBytes() == 4)
    {
        writeIntData<std::int32_t, Long>(data, size, os, id);
    }
    else if (id.numBytes() == 8)
    {
        writeIntData<std::int64_t, Long>(data, size, os, id);
    }
    else {
        amrex::Error("Don't know how to work with this long type.");
    }
}

void amrex::readLongData(Long* data, std::size_t size, std::istream& is,
                         const IntDescriptor& id)
{
    if (id == FPC::NativeLongDescriptor()) 
    {
        is.read((char*) data, size * id.numBytes());
    } 
    else if (id.numBytes() == 2)
    {
        readIntData<Long, std::int16_t>(data, size, is, id);
    }
    else if (id.numBytes() == 4)
    {
        readIntData<Long, std::int32_t>(data, size, is, id);
    }
    else if (id.numBytes() == 8)
    {
        readIntData<Long, std::int64_t>(data, size, is, id);
    }
    else {
        amrex::Error("Don't know how to work with this long type.");
    }
}

void amrex::writeRealData(const Real* data, std::size_t size, std::ostream& os,
                          const RealDescriptor& rd)
{
    RealDescriptor::convertFromNativeFormat(os, static_cast<Long>(size), data, rd);
}

void amrex::readRealData(Real* data, std::size_t size, std::istream& is,
                         const RealDescriptor& rd)
{
    RealDescriptor::convertToNativeFormat(data, static_cast<Long>(size), is, rd);
}

void amrex::writeFloatData(const float* data, std::size_t size, std::ostream& os,
                           const RealDescriptor& rd)
{
    RealDescriptor::convertFromNativeFloatFormat(os, static_cast<Long>(size), data, rd);
}

void amrex::readFloatData(float* data, std::size_t size, std::istream& is,
                          const RealDescriptor& rd)
{
    RealDescriptor::convertToNativeFloatFormat(data, static_cast<Long>(size), is, rd);
}

void amrex::writeDoubleData(const double* data, std::size_t size, std::ostream& os,
                            const RealDescriptor& rd)
{
    RealDescriptor::convertFromNativeDoubleFormat(os, static_cast<Long>(size), data, rd);
}

void amrex::readDoubleData(double* data, std::size_t size, std::istream& is,
                           const RealDescriptor& rd)
{
    RealDescriptor::convertToNativeDoubleFormat(data, static_cast<Long>(size), is, rd);
}
