
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_MPI
#define MY_USE_MPI 1 
#else 
#define MY_USE_MPI 0
#endif 

using namespace amrex;

void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile1=f1 infile2=f2, ngrow=ngrow" << std::endl;
    exit(1);
}

int main(int argc, char* argv[])
{
    static_assert(MY_USE_MPI!=1, "cannot work with MPI" ); 
    amrex::Initialize(argc,argv);
    {
        if (argc < 3)
            print_usage(argc,argv);

        std::string name1, name2;
        int ngrow = -1;
        {
        	ParmParse pp;
        	pp.get("infile1", name1);
          	pp.get("infile2", name2);
        	pp.query("ngrow", ngrow);
        }

        MultiFab mf1, mf2;

        Print() << "Reading " << name1 << std::endl;
        VisMF::Read(mf1, name1);

        Print() << "Reading " << name2 << std::endl;
        VisMF::Read(mf2, name2);

        if (ngrow < 0) ngrow = std::min(mf1.nGrow(), mf2.nGrow());
        if (ngrow > mf1.nGrow())
            Abort(" ngrow bigger than infile1's ngrow! ");
         if (ngrow > mf2.nGrow())
            Abort(" ngrow bigger than infile2's ngrow! ");
   
        if (mf1.boxArray() != mf2.boxArray()) {
        	Abort("The two multifabs have different BoxArray");
        }

        const int ncomp = mf1.nComp();

/*        std::vector<Real> mf1_min(ncomp);
        std::vector<Real> mf1_max(ncomp);
        std::vector<Real> mf2_min(ncomp);
        std::vector<Real> mf2_max(ncomp);

        for (int icomp = 0; icomp < ncomp; ++icomp) {
            mf1_min[icomp] = mf1.min(icomp,ngrow);
            mf1_max[icomp] = mf1.max(icomp,ngrow);
            mf2_min[icomp] = mf2.min(icomp,ngrow);
            mf2_max[icomp] = mf2.max(icomp,ngrow);
        }
*/ 

        MultiFab mfdiff(mf1.boxArray(), mf1.DistributionMap(), ncomp, ngrow);

        MultiFab::Copy(mfdiff, mf1, 0, 0, ncomp, ngrow);
        MultiFab::Subtract(mfdiff, mf2, 0, 0, ncomp, ngrow);

        for (int icomp = 0; icomp < ncomp; ++icomp) {
            Print() << "Component " << icomp << std::endl; 
            Print() << "diff Min,max: " << mfdiff.min(icomp,ngrow) 
    		  << " , " << mfdiff.max(icomp,ngrow) << std::endl;
        }
/*    	if (ncomp > 1) {
    	    Print() << " for component " << icomp;
    	}
        	Print() << " 1st mf min,max: " << mf1_min[icomp]
    		  << ", " << mf1_max[icomp]
                      << ", 2nd mf min,max: ";
        	Print() << mf2_min[icomp]
    	    	  << ", " << mf2_max[icomp] << "\n";
         } */

         Print() << "Writing mfdiff" << std::endl;
         VisMF::Write(mfdiff, "mfdiff");
    }
    amrex::Finalize();
}

