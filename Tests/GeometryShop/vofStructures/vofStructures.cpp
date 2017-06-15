#include <cmath>
#include <type_traits>

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_FlatPlateGeom.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBISBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BLProfiler.H"
#include <AMReX_BLFort.H>
#include <AMReX_EBDebugDump.H>
using namespace amrex;

#include <Node.H>

extern "C"
{
    void do_eb_work(const int* lo, const int* hi,
                    CNode* cnodes, const int* cnum,
                    D_DECL(FNode* fnodes0,
                           FNode* fnodes1,
                           FNode* fnodes2), const int* fnum,
                    const BL_FORT_IFAB_ARG(mask));
}

static void
Copy(intDIM& out, const IntVect& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

static void
Copy(IntVect& out, const intDIM& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

#if 0
static void
dumpInts(const intDIM& iv)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        std::cout << iv[d] << " ";
    }
    std::cout << '\n';
}
#endif

void
get_EGLG(EBLevelGrid& eblg)
{
    Real radius = 0.5;
    Real domlen = 1;
    std::vector<Real> centervec(SpaceDim);
    std::vector<int>  ncellsvec(SpaceDim);

    ParmParse pp;

    std::string geom_type;

    pp.get("domain_length", domlen);                     
    pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
    IntVect ivlo = IntVect::TheZeroVector();
    IntVect ivhi;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
        ivhi[idir] = ncellsvec[idir] - 1;
    }
    Box domain(ivlo, ivhi);
    RealVect origin = RealVect::Zero;
    Real dx = domlen/ncellsvec[0];

    int verbosity = 0;
    pp.get("verbosity", verbosity);

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    pp.get("geom_type",geom_type);
    if (geom_type == "sphere")
    {
        pp.get(   "sphere_radius", radius);
        pp.getarr("sphere_center", centervec, 0, SpaceDim);
        RealVect center;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            center[idir] = centervec[idir];
        }

        //bool insideRegular = false;
        bool insideRegular = true;
        SphereIF sphere(radius, center, insideRegular);

        GeometryShop gshop(sphere, verbosity);
        ebisPtr->define(domain, origin, dx, gshop);
    }
    else if (geom_type == "flat_plate") 
    {
        std::vector<Real>  platelovec(SpaceDim);
        std::vector<Real>  platehivec(SpaceDim);

        int normalDir;
        Real plateLoc;
    
        pp.getarr("plate_lo", platelovec, 0, SpaceDim);
        pp.getarr("plate_hi", platehivec, 0, SpaceDim);
        pp.get("plate_location", plateLoc);
        pp.get("plate_normal", normalDir);

        RealVect plateLo, plateHi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            plateLo[idir] = platelovec[idir];
            plateHi[idir] = platehivec[idir];
        }
        FlatPlateGeom flat_plate(normalDir, plateLoc, plateLo, plateHi);
        ebisPtr->define(domain, origin, dx, flat_plate);
    }
    else
    {
        amrex::Abort("Unknown geom_type");
    }

    BoxArray ba(domain);
    int maxboxsize = 16;
    pp.query("maxboxsize",maxboxsize);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    int nGrowEBLG = 2;
    eblg.define(ba, dm, domain, nGrowEBLG);
}

struct tface
{
    /*
      A simple class to simplify uniquifying faces
     */
    tface() : mL(-1), mR(-1) {}
    tface(int L, int R, const SideIterator& s) : mL(L), mR(R) {if (s()==Side::Lo) flip();}
    bool operator< (const tface& rhs) const {
	if (mL == rhs.mL)
	{
	    return mR < rhs.mR;
	}
	return mL < rhs.mL;
    }
    void flip ()
	{
	    int t=mL;
	    mL = mR;
	    mR = t;
	}
    int mL, mR;
};

static void
AssignFaces(CNode&              cnode,
	    const FNode&        fnode,
	    int                 idir,
	    const SideIterator& sit)
{
    int iside = sit() == Side::Lo ? 0 : 1;
    for (int L=0; L<cnode.nCells; ++L) {
	for (int nc=0; nc<cnode.cells[L].Nnbr[idir][iside]; ++nc) {
	    
	    cnode.cells[L].faceID[idir][iside][nc] = REGULAR_FACE;

	    int c1=L;
	    int c2=cnode.cells[L].nbr[idir][iside][nc];

	    for (int nf=0; nf<fnode.nFaces; ++nf) {
		FNode::CutFace cf = fnode.faces[nf];
		int f1 = cf.cellLo;
		int f2 = cf.cellHi;		
		if ( (f1==c1 && f2==c2) || (f1==c2 && f2==c1) ) {
		    //cnode.cells[L].faceID[idir][iside][nc] = nf; // Use local numbering
		    cnode.cells[L].faceID[idir][iside][nc] = cf.ebFaceID; // Use box-wide numbering
		}
	    }
	}
    }
}

static void
BuildFortranGraphNodes(std::map<int,std::vector<CNode> >& graphCNodes,
                       std::map<int,std::array<std::vector<FNode>, BL_SPACEDIM> >& graphFNodes,
                       iMultiFab& ebmask,
                       const EBLevelGrid& eblg)
{
    for (MFIter mfi(ebmask); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        EBISBox ebis_box = eblg.getEBISL()[mfi];
        const Box gbox = ebmask[mfi].box() & ebis_box.getDomain();

	int ebCellID = 0; // Reset to zero for each box
        if (!ebis_box.isAllRegular())
        {
	    std::array<std::map<IntVect, std::set<tface> >, BL_SPACEDIM> tf;
	    
            const EBGraph& ebgr = ebis_box.getEBGraph();
            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                const IntVect iv = bit();
                if (ebis_box.isIrregular(iv))
                {
		    // Set up cnode data structure
                    std::vector<VolIndex> gbox_vofs = ebis_box.getVoFs(iv);
                    if (gbox_vofs.size() > 0)
                    {
                        graphCNodes[gid].push_back(CNode());
                        CNode& gn = graphCNodes[gid].back();
                        Copy(gn.iv,iv);
                        gn.nCells = gbox_vofs.size();

                        for (int icc=0; icc<gn.nCells; ++icc)
                        {
                            gn.cells[icc].ebCellID = ebCellID++;
                            const VolIndex& vof = gbox_vofs[icc];
                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int iside = sit() == Side::Lo ? 0 : 1;
                                    std::vector<FaceIndex> faces = ebgr.getFaces(vof,idir,sit());

				    IntVect iv_nbr = iv + sign(sit())*BASISV(idir);
				    IntVect iv_face = iv + iside*BASISV(idir);
				    
                                    int nValid = 0;
                                    for (int iface=0; iface<faces.size(); ++iface)
                                    {
                                        if (!ebis_box.isCovered(iv_nbr))
                                        {
                                            if (ebis_box.isRegular(iv_nbr))
                                            {
                                                gn.cells[icc].nbr[idir][iside][nValid] = REGULAR_CELL;
                                            }
                                            else 
                                            {
                                                gn.cells[icc].nbr[idir][iside][nValid] = faces[iface].cellIndex(sit());
                                            }
                                            // Add this face to the set of faces (and including faces between cut and regular cells)
					    tf[idir][iv_face].insert(tface(icc,gn.cells[icc].nbr[idir][iside][nValid],sit));
                                            nValid++;
					    BL_ASSERT(nValid<NCELLMAX);
                                        }
                                    }
                                    gn.cells[icc].Nnbr[idir][iside] = nValid;
                                }
                            }			    
                        }
                    }
                }
            }

	    const Box& vbox = mfi.validbox();
	    std::array<std::map<IntVect,int>,BL_SPACEDIM> fnodeAtIV;
	    for (int idir = 0; idir < SpaceDim; idir++)
	    {
		int ebFaceID = 0;
		const Box fbox = amrex::surroundingNodes(vbox,idir);
		for (BoxIterator bit(fbox); bit.ok(); ++bit)
		{
		    const IntVect& iv_face = bit();
		    std::map<IntVect, std::set<tface> >::const_iterator it=tf[idir].find(iv_face);
		    if (it!=tf[idir].end()) {
			graphFNodes[gid][idir].push_back(FNode());
			FNode& fn = graphFNodes[gid][idir].back();
			const std::set<tface>& tfiv=it->second;
			fn.nFaces = tfiv.size();
			Copy(fn.iv,iv_face);
			int cnt = 0;
			for (std::set<tface>::iterator fit=tfiv.begin(); fit!=tfiv.end(); ++fit)
			{
			    fn.faces[cnt].cellLo = fit->mL;
			    fn.faces[cnt].cellHi = fit->mR;
                            bool is_face_between_reg_and_cut_cells
                                = fn.faces[cnt].cellLo==REGULAR_CELL
                                || fn.faces[cnt].cellHi==REGULAR_CELL;

                            if ( !is_face_between_reg_and_cut_cells )
                            {
                                fn.faces[cnt].ebFaceID = ebFaceID++;
                                cnt++;
                            }
                            else
                            {
                                fn.faces[cnt].ebFaceID = REGULAR_FACE;
                            }
			}
			if (cnt!=0)
			{
                            // Number the FNode (vs. CutFaces, this is one per IntVect)
			    fnodeAtIV[idir][iv_face] = graphFNodes[gid][idir].size() - 1;
			}
		    }
		}
	    }

	    for (int ic=0; ic<graphCNodes[gid].size(); ++ic)
	    {
		CNode& gn = graphCNodes[gid][ic];
		IntVect iv;
		Copy(iv,gn.iv);
		for (int icc=0; icc<gn.nCells; ++icc)
		{
		    for (int idir = 0; idir < SpaceDim; idir++)
		    {
			for (SideIterator sit; sit.ok(); ++sit)
			{
			    int iside = sit() == Side::Lo ? 0 : 1;
			    IntVect iv_face = iv + iside*BASISV(idir);
			    std::map<IntVect,int>::const_iterator it = fnodeAtIV[idir].find(iv_face);
			    if (it!=fnodeAtIV[idir].end());
			    {
				for (int iface=0; iface<gn.cells[icc].Nnbr[idir][iside]; ++iface)
				{
				    const FNode& fn = graphFNodes[gid][idir][it->second];
				    AssignFaces(gn,fn,idir,sit);
				}
			    }
			}
		    }
		}
	    }
	}

        BaseFab<int>& mask_fab = ebmask[mfi];
        mask_fab.setVal(REGULAR_CELL);
        for (int inode=0; inode<graphCNodes[gid].size(); ++inode)
        {
            IntVect iv;
            Copy(iv,graphCNodes[gid][inode].iv);
            mask_fab(iv,0) = inode;
        }
        if (ebis_box.isAllCovered()) {
            mask_fab.setVal(COVERED_CELL);
        }
        else if (!ebis_box.isAllRegular()) {
            // Not all covered, or all regular, there might be covered cells...
            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                if (ebis_box.isCovered(bit())) {
                    mask_fab(bit(),0) = COVERED_CELL;
                }
            }
        }
    }
}


int myTest()
{

    std::cout << "CNode is POD: " << std::is_pod<CNode>::value << std::endl;
    std::cout << "FNode is POD: " << std::is_pod<FNode>::value << std::endl;
    std::cout << "Size of CNode: " << sizeof(CNode)/sizeof(int) << " ints" << std::endl;
    std::cout << "Size of FNode: " << sizeof(FNode)/sizeof(int) << " ints" << std::endl;

    EBLevelGrid eblg;
    get_EGLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();

    IntVect tilesize(AMREX_D_DECL(10240,8,32));
    std::map<int,std::vector<CNode> > graphCNodes;
    std::map<int,std::array<std::vector<FNode>, BL_SPACEDIM> > graphFNodes;
    int nGrow = 1;
    iMultiFab ebmask(ba,dm,1,nGrow);

    BuildFortranGraphNodes(graphCNodes,graphFNodes,ebmask,eblg);
    
    MultiFab state(ba,dm,1,nGrow);

    for (MFIter mfi(ebmask); mfi.isValid(); ++mfi)
    {
        const Box& vbox = mfi.validbox();
        BaseFab<int>& mask_fab = ebmask[mfi];
        int cnum = graphCNodes[mfi.index()].size();
        int fnum[BL_SPACEDIM];
        for (int idir=0; idir<BL_SPACEDIM; ++idir) 
        {
            fnum[idir] = graphFNodes[idir][mfi.index()].size();
        }
        do_eb_work(vbox.loVect(), vbox.hiVect(),
                   graphCNodes[mfi.index()].data(), &cnum,
                   D_DECL(graphFNodes[0][mfi.index()].data(),
                          graphFNodes[1][mfi.index()].data(),
                          graphFNodes[2][mfi.index()].data()),fnum,
                   BL_TO_FORTRAN_N(mask_fab,0));
    }

    return 0;
}

int
main(int argc,char **argv)
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE_VAR("main()", pmain);

        int retFlag = myTest();

        if (retFlag != 0)
        {
            amrex::Print() << "non zero return detected = " << retFlag << '\n';
            amrex::Print() << "sphere test failed" << '\n';
        }
        else
        {
            amrex::Print() << "sphere test passed" << '\n';
        }

        BL_PROFILE_VAR_STOP(pmain);
    }
    amrex::Finalize();
    return 0;
}
