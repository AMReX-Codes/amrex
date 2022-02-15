
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BoxList.H>
#include <AMReX_BoxArray.H>
#include <fstream>

using namespace amrex;

void test ();
BoxArray readBoxList (const std::string& file, Box& domain);

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    test();
    amrex::Finalize();
}

void test ()
{
    BL_PROFILE("test");
    const int ngrids = 10;
    Vector<Box> domains(ngrids);
    Vector<BoxArray> grids(ngrids);

    for (int igrid=0; igrid < ngrids; ++igrid)
    {
        grids[igrid] = readBoxList("grids/grids_"+std::to_string(igrid+1), domains[igrid]);
    }

    for (int igrid=0; igrid < ngrids; ++igrid)
    {
        amrex::Print() << "grid # " << igrid+1 << ": size " << grids[igrid].size()
                       << " min box " << grids[igrid].minimalBox()
                       << "\n                "
                       << " domain box " << domains[igrid] << "\n";
        BL_PROFILE("BoxList::"+std::to_string(igrid+1));
        BoxList bl;
        bl.complementIn(domains[igrid], grids[igrid]);
        amrex::Print() << "BoxList # " << igrid+1 << ": size " << bl.size() << "\n\n";
    }
}

BoxArray
readBoxList (const std::string& file, Box& domain)
{
    BoxArray retval;

    std::ifstream ifs;
    ifs.open(file.c_str(), std::ios::in);
    ifs >> domain;
    ifs.ignore(1000,'\n');
    retval.readFrom(ifs);

    return retval;
}
