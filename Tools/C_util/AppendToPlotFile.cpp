// --------------------------------------------------------------------
// AppendToPlotFile.cpp
// --------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>

#include <AppendToPlotFile.H>

using namespace amrex;

void
AppendToPlotFile(AmrData&                   amrData,
                 const Vector<MultiFab*>&   mfout,
                 const std::string&         oFile,
                 const Vector<std::string>& nnames,
                 const std::string&         mfBaseName,
                 const std::string&         newHeaderName,
                 bool                       verbose)
{
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    std::string nFileHeader(oFile);
    nFileHeader += '/';
    nFileHeader += newHeaderName;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    std::ifstream is;

    os.precision(17);
    int finestLevel = mfout.size() - 1;
    
    int n_var = -1;

    if (ParallelDescriptor::IOProcessor())
    {
        if (verbose)
            std::cout << "Opening files = " << oFileHeader << " and " << nFileHeader << '\n';

        is.open(oFileHeader.c_str(), std::ios::in|std::ios::binary);
        os.open(nFileHeader.c_str(), std::ios::out|std::ios::binary);
        
        if (os.fail())
            FileOpenFailed(oFileHeader);
        if (is.fail())
            FileOpenFailed(nFileHeader);
        //
        // Start writing new plotfile header.
        //
//         std::cout << "reaching here" << std::endl;
        std::string version;
        is >> version;
        os << version << '\n';
        is >> n_var;
        os << n_var+nnames.size() << '\n';
        Vector<std::string> inames(n_var);
        for (int n = 0; n < n_var; n++) is >> inames[n];
        for (int n = 0; n < n_var; n++) os << inames[n] << '\n';
        for (int n = 0; n < nnames.size(); n++) os << nnames[n] << '\n';
        
        int sdim;
        is >> sdim;
        os << sdim << '\n';
        
        Real time;
        is >> time;
        os << time << '\n';
        
        int oFinestLevel;
        is >> oFinestLevel;
        
        BL_ASSERT(oFinestLevel>=finestLevel);
        os << finestLevel << '\n';
        
        Vector<Real> probLo(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM; i++) is >> probLo[i];
        for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
        os << '\n';
        
        Vector<Real> probHi(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM; i++) is >> probHi[i];
        for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
        os << '\n';
        
        Vector<int> refRatio(oFinestLevel);
        for (int i = 0; i < oFinestLevel; i++) is >> refRatio[i];
        for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
        os << '\n';
        
        Vector<Box> probDomain(oFinestLevel+1);
        for (int i = 0; i <= oFinestLevel; i++) is >> probDomain[i];
        for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
        os << '\n';
        
        int tmpI;
        for (int i = 0; i <= oFinestLevel; i++) is >> tmpI;
        for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
        os << '\n';
        
        Real dx[BL_SPACEDIM];
        for (int i = 0; i <= oFinestLevel; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                is >> dx[k];
            }
            if (i<=finestLevel)
            {
                for (int k = 0; k < BL_SPACEDIM; k++)
                {
                    os << dx[k] << ' ';
                }
                os << '\n';
            }
        }
        
        int coordSys;
        is >> coordSys;
        os << coordSys << '\n';
        
        int bndry;
        is >> bndry;
        os << bndry << '\n'; // The bndry data width.
    }
    //
    // Write out level by level.
    //
    std::string mfBaseName_Unique = mfBaseName; // Possibly modified if same name already used
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        int nGrids = amrData.boxArray(iLevel).size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);

        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        
        if (ParallelDescriptor::IOProcessor())
        {
            int ilev,ngrd;
            Real time;
            is >> ilev >> ngrd >> time;
            os << ilev << ' ' << ngrd << ' ' << time << '\n';

            int tmpI;
            is >> tmpI;
            os << tmpI << '\n';
            
            Real glocl,gloch;
            for (int i = 0; i < nGrids; ++i)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    is >> glocl >> gloch;
                    os << glocl
                       << ' '
                       << gloch
                       << '\n';
                }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            // NOTE: should already exist!
            //
            if (!UtilCreateDirectory(PathName, 0755))
                CreateDirectoryFailed(PathName);
        }
        //
        // Force other processors to wait until directory is "built".
        //
  //       std::cout << "Above Barrier" << std::endl;       
        ParallelDescriptor::Barrier();
  //       std::cout << "Below Barrir" << std::endl;

        std::string RelativePathNameNEW(buf);
        RelativePathNameNEW += '/';
        RelativePathNameNEW += mfBaseName_Unique;
        
        if (ParallelDescriptor::IOProcessor())
        {
            // account for multiple multifabs 
            int currentIndexComp(0);
            while(currentIndexComp < n_var) {
                
                std::string RelativePathNameOLD;
                is >> RelativePathNameOLD;

                // Avoid name clash
                if (RelativePathNameOLD == RelativePathNameNEW) {
                    mfBaseName_Unique += (const std::string)("x");

                    RelativePathNameNEW = buf;
                    RelativePathNameNEW += '/';
                    RelativePathNameNEW += mfBaseName_Unique;
                }

                string mfName(oFile);
                string corrected ;
                mfName += '/';
                mfName += RelativePathNameOLD;
                corrected += RelativePathNameOLD;
                VisMF visMF(mfName);
                os << corrected << '\n';
                std::cout << corrected << std::endl;
                currentIndexComp += visMF.nComp();
            }

            // Reform name in case we had to rename the mf
            //  and add new name after exisiting ones
            RelativePathNameNEW = buf;
            RelativePathNameNEW += '/';
            RelativePathNameNEW += mfBaseName_Unique;
            os << RelativePathNameNEW << '\n';

        }

        //
        // Now build the full relative pathname of the MultiFab.
        //
        PathName += '/';
        PathName += mfBaseName;
//        std::cout << "Aboe Vismf" <<std::endl;
       
        // Write the new multifab
        VisMF::Write(*mfout[iLevel], PathName, VisMF::OneFilePerCPU);
    }
//       std::cout <<"just above file close" << std::endl;   
    os.close();
    is.close();
}
