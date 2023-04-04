#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_iMultiFab.H"

int MyMain();

int main(int argc, char** argv) {
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#else
    amrex::ignore_unused(argc,argv);
#endif
    // Let me throw exceptions for triggering my debugger
    amrex::Initialize(MPI_COMM_WORLD, std::cout, std::cerr, [](const char* msg) { throw std::runtime_error(msg); });
    int ret = MyMain();
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
    return ret;
}

amrex::iMultiFab InitializeMultiFab(const amrex::Box& domain)
{
    amrex::BoxArray ba(domain);
    // do some work across MPI ranks.
    ba.maxSize(8);
    amrex::DistributionMapping dm(ba);
    amrex::iMultiFab mf(ba, dm, 1, 0);
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto array = mf.array(mfi);
        ParallelFor(mfi.tilebox(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            array(i,j,k) = i + j*nx + k*nx*ny;
        });
    }
    return mf;
}

bool ParallelCopyWithItselfIsCorrect(amrex::iMultiFab& mf, const amrex::Box& domain) {
    const amrex::IntVect e_x = amrex::IntVect::TheDimensionVector(0);
    const amrex::Box source_box{domain.smallEnd(), domain.bigEnd() - domain.bigEnd(0) * e_x, domain.ixType()};
    const amrex::Box dest_box{domain.smallEnd() + domain.bigEnd(0) * e_x, domain.bigEnd(), domain.ixType()};
    amrex::NonLocalBC::MultiBlockIndexMapping dtos;
    dtos.offset = dest_box.smallEnd() - source_box.smallEnd();
    amrex::NonLocalBC::ParallelCopy(mf, dest_box, mf, 0, 0, 1, amrex::IntVect{0}, dtos);
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    int fails = 0;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box section = dest_box & mfi.tilebox();
        if (section.isEmpty()) continue;
        auto array = mf.const_array(mfi);
        amrex::LoopOnCpu(section, [&](int i, int j, int k)
        {
            amrex::Dim3 si = dtos(amrex::Dim3{i,j,k});
            int value = si.x + si.y*nx + si.z*nx*ny;
            fails += (array(i,j,k) != value);

            AMREX_ASSERT(fails==0);  // If DEBUG, crash on first error.
        });
    }
    return fails == 0;
}

int GetFaceDir(amrex::IndexType iv)
{
    static_assert(sizeof(amrex::IndexType) == sizeof(unsigned int), "IndexType is not punnable to unsigned int");
    unsigned int value{0};
    std::memcpy(&value, &iv, sizeof(unsigned int));
    return int(((value & 0b001) + (value & 0b010) + (value & 0b100)) >> 1);
}

amrex::Box GetFaceBoundary(const amrex::Box& domain, amrex::Orientation::Side side)
{
    int dir = GetFaceDir(domain.ixType());
    const amrex::IntVect face_normal = amrex::IntVect::TheDimensionVector(dir);
    amrex::Box box{};
    AMREX_ASSERT(side == amrex::Orientation::Side::low || side == amrex::Orientation::Side::high);
    if (side == amrex::Orientation::Side::low) {
        box = amrex::Box{domain.smallEnd(), domain.bigEnd() - domain.bigEnd(dir) * face_normal, domain.ixType()};
    } else {
        box = amrex::Box{domain.smallEnd() + domain.bigEnd(dir) * face_normal, domain.bigEnd(), domain.ixType()};
    }
    return box;
}

bool ParallelCopyFaceToFace(amrex::iMultiFab& dest, const amrex::Box& domain_dest, amrex::Orientation::Side dest_side,
                            const amrex::iMultiFab& src, const amrex::Box& domain_src, amrex::Orientation::Side src_side)
{
    int sdir = GetFaceDir(domain_src.ixType());
    int ddir = GetFaceDir(domain_dest.ixType());
    const amrex::Box src_box = GetFaceBoundary(domain_src, src_side);
    const amrex::Box dest_box = GetFaceBoundary(domain_dest, dest_side);

    // Default construction is identity
    amrex::NonLocalBC::MultiBlockIndexMapping dtos{};
    // Change permutation to get a correct mapping between index types
    std::swap(dtos.permutation[sdir], dtos.permutation[ddir]);
    // Map smallest destination box index as an index in the source space
    const amrex::IntVect dest_smallEnd_in_src = amrex::NonLocalBC::Apply(dtos, dest_box.smallEnd());
    // Compute the offset to get the proper shift in the source space
    dtos.offset =  dest_smallEnd_in_src - src_box.smallEnd();
    // Sanity-Check that the correct box is being mapped
    // Note, this checks index types, too!
    AMREX_ASSERT(amrex::NonLocalBC::Image(dtos, dest_box) == src_box);

    amrex::NonLocalBC::ParallelCopy(dest, dest_box, src, 0, 0, 1, amrex::IntVect{0}, dtos);
    int fails = 0;
    const int nx = domain_src.length(0);
    const int ny = domain_src.length(1);
    for (amrex::MFIter mfi(dest); mfi.isValid(); ++mfi) {
        const amrex::Box section = dest_box & mfi.tilebox();
        if (section.isEmpty()) continue;
        auto darray = dest.const_array(mfi);
        amrex::LoopOnCpu(section, [&](int i, int j, int k)
        {
            amrex::Dim3 si = dtos(amrex::Dim3{i,j,k});
            int value = si.x + si.y*nx + si.z*nx*ny;
            fails += (darray(i,j,k) != value);

            AMREX_ASSERT(fails==0); // If in debug, crash on first error.
        });
    }
    return fails == 0;
}

int MyMain()
{
    using namespace amrex;
    Box domain{IntVect{}, IntVect{AMREX_D_DECL(31, 31, 31)}};
    // Loop over all index types
    for (int i = 0; i < AMREX_D_TERM(2,*2,*2); ++i) {
        AMREX_D_TERM(const auto ix = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b001));,
                     const auto iy = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b010));,
                     const auto iz = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b100));)
        IndexType itype{AMREX_D_DECL(ix, iy, iz)};
        Box converted_domain = convert(domain, itype);
        iMultiFab mf = InitializeMultiFab(converted_domain);
        ParallelCopyWithItselfIsCorrect(mf, converted_domain);
    }
    // Loop over all combinations of face orientations
    const IntVect AMREX_D_DECL(e_x = IntVect::TheDimensionVector(0), e_y = IntVect::TheDimensionVector(1), e_z = IntVect::TheDimensionVector(2));
    IndexType dirs[AMREX_SPACEDIM] = {AMREX_D_DECL(IndexType(e_x), IndexType(e_y), IndexType(e_z))};
    Orientation::Side sides[2] = {Orientation::low, Orientation::high};
    for (int ii = 0; ii < 2*AMREX_SPACEDIM; ++ii) {
        int i = ii / 2;
        Orientation::Side iside = sides[ii % 2];
        for (int jj = 0; jj < 2*AMREX_SPACEDIM; ++jj) {
            int j = jj / 2;
            Orientation::Side jside = sides[j % 2];
            Box converted_domain_x = convert(domain, dirs[i]);
            iMultiFab dest = InitializeMultiFab(converted_domain_x);
            Box converted_domain_y = convert(domain, dirs[j]);
            const iMultiFab src = InitializeMultiFab(converted_domain_y);
            ParallelCopyFaceToFace(dest, converted_domain_x, iside, src, converted_domain_y, jside);
        }
    }
    return 0;
}
