#include <iostream>
#include <stdint.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    uint64_t num_cells = 512;
    int max_grid_size = 128;
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(num_cells - 1,
                                   num_cells - 1,
                                   num_cells - 1));
    const Box domain(domain_lo, domain_hi);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dmap(ba);
    
    MultiFab zhi(ba, dmap, 1, 0);

    std::string file_name = "zhi.bin";
    std::ifstream ifs;
    ifs.open(file_name.c_str(), std::ios::in|std::ios::binary);
    if ( not ifs ) {
        amrex::Print() << "Failed to open file " << file_name << " for reading. \n";
        amrex::Abort();
    }
    Vector<float> values(num_cells*num_cells*num_cells);
    ifs.read((char*) &values[0], num_cells*num_cells*num_cells*sizeof(float));    

    for (MFIter mfi(zhi); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        for (unsigned x = box.smallEnd(0); x <= box.bigEnd(0); ++x) {
            for (unsigned y = box.smallEnd(1); y <= box.bigEnd(1); ++y) {
                for (unsigned z = box.smallEnd(2); z <= box.bigEnd(2); ++z) {
                    IntVect iv(x, y, z);
                    uint64_t index = x*num_cells*num_cells + y*num_cells + z;
                    zhi[mfi](iv) = values[index];
                }
            }
        }
    }

    amrex::VisMF::Write(zhi, "zhi/zhi");
    
    amrex::Finalize();
}
