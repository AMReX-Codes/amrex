// ------------------------------------------------------------
// A test program for ReadAndBcastFile().
// ------------------------------------------------------------
#include <ParallelDescriptor.H>


// ------------------------------------------------------------
int main(int argc, char **argv) {
  BoxLib::Initialize(argc, argv);

  BL_PROFILE_VAR("main()", pmain);


  Array<double> dmSTimes(nStrategies, 0.0);

    ParallelDescriptor::Barrier();

    // Open the checkpoint header file for reading.
    std::string File(filename);
    File += "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    Array<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::stringbuf *pbuf = iss.rdbuf();

    is >> value;



  if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < nStrategies; ++i) {
        std::cout << std::endl << "Total times:" << std::endl;
	std::cout << dmSNames[i] << " time = " << dmSTimes[i] << std::endl;
      }
      std::cout << std::endl << std::endl;
  }




    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();

    return 0;
}
