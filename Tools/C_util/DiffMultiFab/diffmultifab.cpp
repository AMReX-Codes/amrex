
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
print_usage (int, char* argv[])
{
    Print()<<"\n"
           <<"This program differences two MultiFabs.\n\n"
           << "usage:\n"
           << argv[0] << " infile1=mf1 infile2=mf2 ngrow=ngrow\n"
           << std::endl;
    exit(1);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        if (argc < 3) {
            print_usage(argc,argv);
        }

        const std::string farg = amrex::get_command_argument(1);
        if (farg == "-h" || farg == "--help")
        {
            print_usage(argc, argv);
        }

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

        if (ngrow < 0) {
            ngrow = std::min(mf1.nGrow(), mf2.nGrow());
        }
        if (ngrow > mf1.nGrow()) {
            Abort(" ngrow bigger than infile1's ngrow! ");
        }
        if (ngrow > mf2.nGrow()) {
            Abort(" ngrow bigger than infile2's ngrow! ");
        }
        if (mf1.boxArray() != mf2.boxArray()) {
            Abort("The two multifabs have different BoxArray");
        }

        const int ncomp = mf1.nComp();

        std::vector<Real> mf1_min(ncomp);
        std::vector<Real> mf1_max(ncomp);
        std::vector<Real> mf2_min(ncomp);
        std::vector<Real> mf2_max(ncomp);

        for (int icomp = 0; icomp < ncomp; ++icomp) {
            mf1_min[icomp] = mf1.min(icomp,ngrow);
            mf1_max[icomp] = mf1.max(icomp,ngrow);
            mf2_min[icomp] = mf2.min(icomp,ngrow);
            mf2_max[icomp] = mf2.max(icomp,ngrow);
        }

        MultiFab mfdiff(mf1.boxArray(), mf2.DistributionMap(), ncomp, ngrow);
#ifdef AMREX_USE_MPI
        {
            MultiFab tmp(mf1.boxArray(), mf1.DistributionMap(), ncomp, ngrow);
            MultiFab::Copy(tmp, mf1, 0, 0, ncomp, ngrow);
            mfdiff.Redistribute(tmp, 0, 0, ncomp, IntVect(ngrow));
        }
#else
        MultiFab::Copy(mfdiff, mf1, 0, 0, ncomp, ngrow);
#endif
        MultiFab::Subtract(mfdiff, mf2, 0, 0, ncomp, ngrow);

        for (int icomp = 0; icomp < ncomp; ++icomp) {
            Real mn = mfdiff.min(icomp,ngrow);
            Real mx = mfdiff.max(icomp,ngrow);
            if (ncomp > 1) {
                Print() << "Component " << icomp << "\n";
            }
            Print() << "    Min and max of the diff are " << mn << " and " << mx << "\n";
            if (mn != 0.0) {
                Print() << "    Min Index: " << mfdiff.minIndex(icomp,ngrow) << "\n";
            }
            if (mx != 0.0) {
                Print() << "    Max Index: " << mfdiff.maxIndex(icomp,ngrow) << "\n";
            }
            Print() << "    Min and max of 1st mf are " << mf1_min[icomp]
                    << " and " << mf1_max[icomp] << "\n";
            Print() << "    Min and max of 2nd mf are " << mf2_min[icomp]
                    << " and " << mf2_max[icomp] << "\n";
        }

        Print() << "Writing mfdiff" << "\n";
        VisMF::Write(mfdiff, "mfdiff");
    }
    amrex::Finalize();
}
