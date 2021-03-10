#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_iMultiFab.H"

int MyMain();

int main(int argc, char** argv) {
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
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
    amrex::DistributionMapping dm(ba);
    amrex::iMultiFab mf(ba, dm, 1, 0);
    const int nx = domain.bigEnd(0) + 1;
    const int ny = domain.bigEnd(1) + 1;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto array = mf.array(mfi);
        ParallelFor(mfi.tilebox(), [=](int i, int j, int k) 
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
    int fails = 0;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box section = dest_box & mfi.tilebox();
        if (section.isEmpty()) continue; 
        auto array = mf.const_array(mfi);
        amrex::LoopOnCpu(section, [&](int i, int j, int k)
        {
            amrex::Dim3 si = dtos(amrex::Dim3{i,j,k});
            fails += (array(i,j,k) != array(si.x, si.y, si.z));
        });
    }
    return fails == 0;
}

int GetFaceDir(amrex::IndexType iv)
{
    static_assert(sizeof(amrex::IndexType) == sizeof(unsigned int), "IndexType is not punnable to unsigned int");
    unsigned int value{0};
    std::memcpy(&value, &iv, sizeof(unsigned int));
    return ((value & 0b001) + (value & 0b010) + (value & 0b100)) >> 1;
}

bool ParallelCopyFaceToFace(amrex::iMultiFab& dest, const amrex::Box& domain_dest, const amrex::iMultiFab& src, const amrex::Box& domain_src) {
    int ds = GetFaceDir(domain_src.ixType());
    int dd = GetFaceDir(domain_dest.ixType());
    const amrex::IntVect face_normal_src = amrex::IntVect::TheDimensionVector(ds);
    const amrex::IntVect face_normal_dest = amrex::IntVect::TheDimensionVector(dd);
    const amrex::Box src_box{domain_src.smallEnd(), domain_src.bigEnd() - domain_src.bigEnd(ds) * face_normal_src, domain_src.ixType()};
    const amrex::Box dest_box{domain_dest.smallEnd(), domain_dest.bigEnd() - domain_dest.bigEnd(dd) * face_normal_dest, domain_dest.ixType()};
    amrex::NonLocalBC::MultiBlockIndexMapping dtos;
    std::swap(dtos.permutation[ds], dtos.permutation[dd]);
    amrex::Print() << dest_box << '\n';
    amrex::Print() << src_box << '\n';
    amrex::Print() << amrex::NonLocalBC::Image(dtos, dest_box) << '\n';
    AMREX_ASSERT(amrex::NonLocalBC::Image(dtos, dest_box) == src_box);
    if (amrex::NonLocalBC::Image(dtos, dest_box) != src_box) {
        return false;
    }
    amrex::NonLocalBC::ParallelCopy(dest, dest_box, src, 0, 0, 1, amrex::IntVect{0}, dtos);
    int fails = 0;
    for (amrex::MFIter mfi(dest); mfi.isValid(); ++mfi) {
        const amrex::Box section = dest_box & mfi.tilebox();
        if (section.isEmpty()) continue; 
        auto darray = dest.const_array(mfi);
        auto sarray = src.const_array(mfi);
        amrex::LoopOnCpu(section, [&](int i, int j, int k)
        {
            amrex::Dim3 si = dtos(amrex::Dim3{i,j,k});
            fails += (darray(i,j,k) != sarray(si.x, si.y, si.z));
        });
    }
    return fails == 0;
}


int MyMain()
{
    using namespace amrex;
    Box domain{IntVect{}, IntVect{AMREX_D_DECL(1, 1, 1)}};
    // Loop over all index types
    for (int i = 0; i < AMREX_D_TERM(2,*2,*2); ++i) {
        const auto ix = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b001));
        const auto iy = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b010));
        const auto iz = static_cast<IndexType::CellIndex>(static_cast<bool>(i & 0b100));
        IndexType itype{AMREX_D_DECL(ix, iy, iz)};
        Box converted_domain = convert(domain, itype);
        iMultiFab mf = InitializeMultiFab(converted_domain);
        if (!ParallelCopyWithItselfIsCorrect(mf, converted_domain)) {
            return 1;
        }
    }
    // Loop over face directions
    const IntVect AMREX_D_DECL(e_x = IntVect::TheDimensionVector(0), e_y = IntVect::TheDimensionVector(1), e_z = IntVect::TheDimensionVector(2));
    IndexType dirs[AMREX_SPACEDIM] = {AMREX_D_DECL(IndexType(e_x), IndexType(e_y), IndexType(e_z))};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            Box converted_domain_x = convert(domain, dirs[i]);
            iMultiFab dest = InitializeMultiFab(converted_domain_x);
            Box converted_domain_y = convert(domain, dirs[j]);
            const iMultiFab src = InitializeMultiFab(converted_domain_y);
            if (!ParallelCopyFaceToFace(dest, converted_domain_x, src, converted_domain_y)) {
                return 1;
            }
        }
    }
    return 0;
}