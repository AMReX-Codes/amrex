
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <functional>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

    {
        Vector<int> vec({1});//,2,3,4,5,6,7,7,7});
        amrex::RemoveDuplicates<int,std::hash<int> >(vec);
        for (int i = 0; i < vec.size(); ++i) {
            amrex::Print() << " " << i << ": " << vec[i] << std::endl;
        }
    }

    amrex::Finalize();
}

