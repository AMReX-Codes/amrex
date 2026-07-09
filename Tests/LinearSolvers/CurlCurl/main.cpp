#include "MyTest.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        {
            ParmParse pp;
            Vector<int> variable_alpha;
            Vector<int> variable_beta;
            if (pp.contains("variable_alpha")) {
                int t;
                pp.get("variable_alpha", t);
                variable_alpha.push_back(t);
            } else {
                variable_alpha.push_back(0);
                variable_alpha.push_back(1);
            }
            if (pp.contains("variable_beta")) {
                int t;
                pp.get("variable_beta", t);
                variable_beta.push_back(t);
            } else {
                variable_beta.push_back(0);
                variable_beta.push_back(1);
            }
            for (auto alpha : variable_alpha) {
                for (auto beta : variable_beta) {
                    pp.add("variable_alpha", alpha);
                    pp.add("variable_beta", beta);
                    amrex::Print() << "\nTesting variable_alpha = " << alpha << " variable_beta = " << beta << "\n";
                    MyTest mytest;
                    mytest.solve();
                }
            }
        }
    }

    amrex::Finalize();
}
