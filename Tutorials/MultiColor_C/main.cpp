
#include <MultiFab.H>
#include <ParmParse.H>
#include <VisMF.H>

namespace {
    int ncomp = 16;
}

void f(MultiFab& mf, const MultiFab& mf2);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BoxArray ba, ba2;
    {
	ParmParse pp;
	
	int n_cell, max_grid_size;
	pp.get("n_cell", n_cell);
	pp.get("max_grid_size", max_grid_size);

	Box domain(IntVect(D_DECL(       0,       0,       0)),
		   IntVect(D_DECL(n_cell-1,n_cell-1,n_cell-1)));

	ba.define(domain);
	ba.maxSize(max_grid_size);

	ba2.define(BoxLib::grow(domain, -n_cell/4));
	ba2.maxSize(max_grid_size/2);
    }

    MultiFab mf(ba, ncomp, 0);
    mf.setVal(1.0);

    MultiFab mf2(ba2, ncomp, 0);
    mf2.setVal(2.0);

    int ncomm = ParallelDescriptor::NColors();
    // Note that the number of sub-communicators is controlled by runtime parameter boxlib.ncolors.

    int ncomp_sub = ncomp / ncomm;
    BL_ASSERT(ncomm*ncomp_sub == ncomp);

    //  Build MultiFabs in sub-communicators and copy data into them
    PArray<MultiFab> mf_sub(ncomm, PArrayManage);
    PArray<MultiFab> mf2_sub(ncomm, PArrayManage);
    for (int i = 0; i < ncomm; ++i) {
	mf_sub.set(i, new MultiFab(ba, ncomp_sub, 0, ParallelDescriptor::Color(i)));
	mf2_sub.set(i, new MultiFab(ba2, ncomp_sub, 0, ParallelDescriptor::Color(i)));
	mf_sub[i].copy(mf, ncomp_sub*i, 0, ncomp_sub);
	mf2_sub[i].copy(mf2, ncomp_sub*i, 0, ncomp_sub);
    }

    // Do some work in sub-communicators
    int icolor = ParallelDescriptor::SubCommColor().to_int();
    f(mf_sub[icolor], mf2_sub[icolor]);

    // Copy data back into global MultiFab
    for (int i = 0; i < ncomm; ++i) {
	mf.copy(mf_sub[i], 0, ncomp_sub*i, ncomp_sub);
    }    

    VisMF::Write(mf, "mf");

    BoxLib::Finalize();
}

void f(MultiFab& mf, const MultiFab& mf2)
{
    BL_ASSERT(mf.color() == mf2.color());

    int icolor = mf.color().to_int();

    // Let's build a new mf within the subcommunicator
    MultiFab mftmp(mf.boxArray(),  mf.nComp(), mf.nGrow(), mf.DistributionMap());
    mftmp.setVal(icolor);

    // Let's do some parallel copy
    if (icolor == 0) {
	mftmp.copy(mf2);
	mftmp.copy(mf2);  // For testing SeqNum(), we do parallel copy twice on Color 0
    } else {
	mftmp.copy(mf2);
    }

    MultiFab::Add(mf, mftmp, 0, 0, mf.nComp(), 0);
}

