#include <iostream>
#include <cstring>
//#include <unistd.h>

#include "AMReX_Print.H"
#include "AMReX_SPMD.H"


namespace amrex
{

/// these should work independent of MPI's existence
  template < >
  void linearIn<float>(float& a_outputT, const void* const a_inBuf)
  {
    //this fandango is to avoid unaligned accesses
    char realBuf[sizeof(float)];
    memcpy(realBuf, a_inBuf, sizeof(float));
    float* buffer = (float*)realBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearOut<float>(void* const a_outBuf, const float& a_inputT)
  {
    //this fandango is to avoid unaligned accesses
    char realBuf[sizeof(float)];
    float* realPtr = (float*)realBuf;
    *realPtr = a_inputT;
    memcpy(a_outBuf, realBuf, sizeof(float));
  }

  template < >
  int linearSize<float>(const float& inputfloat)
  {
    return sizeof(float);
  }

  template < >
  void linearIn<double>(double& a_outputT, const void* const a_inBuf)
  {
    //this fandango is to avoid unaligned accesses
    char realBuf[sizeof(double)];
    memcpy(realBuf, a_inBuf, sizeof(double));
    double* buffer = (double*)realBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearOut<double>(void* const a_outBuf, const double& a_inputT)
  {
    //this fandango is to avoid unaligned accesses
    char realBuf[sizeof(double)];
    double* realPtr = (double*)realBuf;
    *realPtr = a_inputT;
    memcpy(a_outBuf, realBuf, sizeof(double));
  }

  template < >
  int linearSize<double>(const double& inputdouble)
  {
    return sizeof(double);
  }

  template < >
  void linearIn<int>(int& a_outputT, const void* const a_inBuf)
  {
    int* buffer = (int*)a_inBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearIn<long long>(long long& a_outputT, const void* const a_inBuf)
  {
    long long* buffer = (long long*)a_inBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearIn<unsigned long long>(unsigned long long& a_outputT, const void* const a_inBuf)
  {
    unsigned long long* buffer = (unsigned long long*)a_inBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearOut<int>(void* const a_outBuf, const int& a_inputT)
  {
    int* buffer = (int*)a_outBuf;
    *buffer = a_inputT;
  }

  template < >
  void linearOut<long long>(void* const a_outBuf, const long long& a_inputT)
  {
    long long* buffer = (long long*)a_outBuf;
    *buffer = a_inputT;
  }

  template < >
  void linearOut<unsigned long long>(void* const a_outBuf, const unsigned long long& a_inputT)
  {
    unsigned long long* buffer = (unsigned long long*)a_outBuf;
    *buffer = a_inputT;
  }

  template < >
  int linearSize<int>(const int& a_input)
  {
    return sizeof(int);
  }

  template < >
  int linearSize<unsigned long long>(const unsigned long long& a_input)
  {
    return sizeof(unsigned long long);
  }

  template < >
  int linearSize<long long>(const long long& a_input)
  {
    return sizeof(long long);
  }

  template < >
  void linearIn<long>(long& a_outputT, const void* const a_inBuf)
  {
    long* buffer = (long*)a_inBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearOut<long>(void* const a_outBuf, const long& a_inputT)
  {
    long* buffer = (long*)a_outBuf;
    *buffer = a_inputT;
  }

  template < >
  int linearSize<long>(const long& a_input)
  {
    return sizeof(long);
  }

  template < >
  void linearIn<unsigned long>(unsigned long& a_outputT, const void* const a_inBuf)
  {
    unsigned long* buffer = (unsigned long*)a_inBuf;
    a_outputT = *buffer;
  }

  template < >
  void linearOut<unsigned long>(void* const a_outBuf, const unsigned long& a_inputT)
  {
    unsigned long* buffer = (unsigned long*)a_outBuf;
    *buffer = a_inputT;
  }

  template < >
  int linearSize<unsigned long>(const unsigned long& a_input)
  {
    return sizeof(unsigned long);
  }

// std::string specialization.
  template <>
  void linearIn<std::string>(std::string& a_outputT, const void* const a_inBuf)
  {
    int* intBuffer = (int*)a_inBuf;
    int length = intBuffer[0];
    if (length > 0)
    {
      const char* charBuffer = (const char*)(&intBuffer[1]);
      a_outputT.assign(charBuffer, length);
    }
    else a_outputT = "";
  }
  template <>
  void linearOut<std::string>(void* const a_outBuf, const std::string& a_inputT)
  {
    int* intBuffer = (int*)a_outBuf;
    intBuffer[0] = a_inputT.length();
    if (a_inputT.length() > 0)
    {
      char* charBuffer = (char*)(&intBuffer[1]);
      std::copy(a_inputT.begin(), a_inputT.end(), charBuffer);
    }
  }
  template <>
  int linearSize<std::string>(const std::string& a_input)
  {
    // A string is stored as its length + its data.
    return sizeof(int) + a_input.length() * sizeof(char);
  }

//std::vector<int>  specialization
  template < > int linearSize(const std::vector<int>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<int>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<int>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

//std::vector<unsigned long long>  specialization
  template < > int linearSize(const std::vector<unsigned long long>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<unsigned long long>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<unsigned long long>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

//std::vector<long>  specialization
  template < > int linearSize(const std::vector<long>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<long>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<long>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

//std::vector<Real>  specialization
  template < > int linearSize(const std::vector<float>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<float>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<float>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

  template < > int linearSize(const std::vector<double>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<double>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<double>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

//std::vector<std::string>  specialization
  template < > int linearSize(const std::vector<std::string>& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<std::string>& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<std::string>& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

//std::vector<std::vector<int> >  specialization
  template < > int linearSize(const std::vector<std::vector<int> >& a_input)
  {
    return linearListSize(a_input);
  }
  template < > void linearIn(std::vector<std::vector<int> >& a_outputT, const void* const inBuf)
  {
    linearListIn(a_outputT, inBuf);
  }
  template < > void linearOut(void* const a_outBuf, const std::vector<std::vector<int> >& a_inputT)
  {
    linearListOut(a_outBuf, a_inputT);
  }

}
