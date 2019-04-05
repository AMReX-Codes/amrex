#include <algorithm> 
#include <string> 
#include <iostream>
#include <iomanip>
#include <fstream>

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_DataServices.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex; 

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

static
void
PrintUsage (char* progName)
{
    std::cout << "\nUsage:\n"
         << progName
         << "\n\tinfile = inputFileName"
         << "\n\t[-help]"
         << "\n\n";
    exit(1);
}

static
void
WriteFab (std::ostream&         os,
          const FArrayBox& fab,
          const char*      name)

{
    int nx = fab.box().length(0);
    int ny = fab.box().length(1);
    int dim = BL_SPACEDIM;


#if (BL_SPACEDIM == 2)
    for (int i = 0; i < nx; i++)
     for (int j = 0; j < ny; j++)
     {
       int index = j*nx + i;
       const Real * ptr = fab.dataPtr();
       os.write((char*)(ptr+index),sizeof(Real));
     }

#elif (BL_SPACEDIM == 3)
    int nz = fab.box().length(2);
    for (int i = 0; i < nx; i++)
     for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
      {
        int index = k*(nx*ny) + j*nx + i;
        const Real * ptr = fab.dataPtr();
        os.write((char*)(ptr+index),sizeof(Real));
      }
#endif
}


int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if(argc == 1)
      PrintUsage(argv[0]);

    if (ParallelDescriptor::NProcs() > 1){
      amrex::Error("This is an inherently serial program!");
    }
    ParmParse pp;

    std::string name; 
    pp.get("infile", name); 

    //
    // MatLab expects native floating-point format. 
    //
    FArrayBox::setFormat(FABio::FAB_NATIVE); 
      
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size); 
    char buf[128]; 
    std::string file = name; 
    file += ".mat"; 
    std::ofstream os; 
    
    os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size()); 
    
    os.open(file.c_str(), std::ios::out|std::ios::binary); 
    
    if(os.fail()) {
      amrex::FileOpenFailed(file); 
    }


    MultiFab in; 
    VisMF::Read(in, name); 
    const BoxArray ba = in.boxArray(); 

    //Fake a xlo to xhi
    for(int i = 0; i < ba.size(); ++i)
    {
      const Box& b = ba[i]; 
      Real xlo[BL_SPACEDIM], xhi[BL_SPACEDIM];
      for(int d = 0; d <BL_SPACEDIM; ++d)
      {
        xlo[d] = b.loVect()[d]; 
        xhi[d] = b.hiVect()[d]; 
        os.write((char*)&xlo,sizeof(Real)); 
        os.write((char*)&xhi,sizeof(Real)); 
      }         
    }
    //Write the Fab Data
    for(int i = 0; i < ba.size(); ++i)
    {
      WriteFab(os, in[i], buf);   
    }
    os.close();
  }
  Finalize(); 
}
