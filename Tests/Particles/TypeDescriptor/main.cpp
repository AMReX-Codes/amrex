#include <iostream>
#include <string>
#include <cstring>

#include "AMReX.H"
#include "AMReX_Print.H"
#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"
#include "AMReX_Vector.H"

#include "Convert.H"

using namespace amrex;

void writeRealData(const Real* data, std::size_t size, std::ostream& os)
{
    return;
}

void readRealData(Real* data, std::size_t size, std::istream& is,
                  const RealDescriptor& rd)
{
    return;
}

void writeIntData(const int* data, std::size_t size, std::ostream& os,
                  const IntDescriptor& id = FPC::NativeIntDescriptor())
{
    if (id == FPC::NativeIntDescriptor())
    {
        os.write((char*) data, size*sizeof(int));
    }
    else if (id.numBytes() == 4)
    {
        writeIntData<std::int32_t>(data, size, os, id);
    }
    else if (id.numBytes() == 8)
    {
        writeIntData<std::int64_t>(data, size, os, id);
    }
    else {
        amrex::Error("Don't know how to work with this integer type.");
    }
}

void readIntData(int* data, std::size_t size, std::istream& is,
                 const IntDescriptor& id)
{
    if (id == FPC::NativeIntDescriptor()) 
    {
        is.read((char*) data, size * id.numBytes());
    } 
    else if (id.numBytes() == 4)
    {
        readIntData<std::int32_t>(data, size, is, id);
    }
    else if (id.numBytes() == 8)
    {
        readIntData<std::int64_t>(data, size, is, id);
    }
    else {
        amrex::Error("Don't know how to work with this integer type.");
    }
}

void testIntIO(const IntDescriptor& id) {
  std::string file_name = "int_data.dat";
 
  amrex::Vector<int> idata_out;
  for (int i = 0; i < 100; ++i) {
      idata_out.push_back(i);
  }

  std::ofstream ofs;
  ofs.open(file_name.c_str(), std::ios::out|std::ios::binary);
  writeIntData(idata_out.data(), idata_out.size(), ofs, id);
  ofs.close();

  amrex::Vector<int> idata_in(100);

  std::ifstream ifs;
  ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
  readIntData(idata_in.data(), idata_in.size(), ifs, id);

  for (int i = 0; i < 100; ++i) {
      AMREX_ALWAYS_ASSERT(idata_in[i] == idata_out[i]);
  }
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  IntDescriptor little32(4, IntDescriptor::ReverseOrder);
  IntDescriptor little64(8, IntDescriptor::ReverseOrder);
  IntDescriptor big32(4, IntDescriptor::NormalOrder);
  IntDescriptor big64(8, IntDescriptor::NormalOrder);

  testIntIO(little32);
  testIntIO(little64);
  testIntIO(big32);
  testIntIO(big64);

  amrex::Finalize();
}
