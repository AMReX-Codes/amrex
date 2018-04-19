

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>

using std::cout;
using std::endl;

#include <unistd.h>

#ifdef USE_TEC_BIN_IO
#include "TECIO.h"
#endif

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Geometry.H"

struct Node
{
    enum typeEnum{INIT=0, COVERED=1, VALID=2};
    Node()
        : level(-1), iv(IntVect(AMREX_D_DECL(-1,-1,-1))), grid(-1), type(Node::INIT) {}
    Node(const IntVect& idx, int lev, int grd, typeEnum typ = INIT)
      : level(lev), iv(idx), grid(grd), type(typ) {}
    inline bool operator< (const Node& rhs) const
        {
            if (level < rhs.level) return true;
            if ((level == rhs.level) && iv < rhs.iv) return true;
            return false;
        }
    inline bool operator!= (const Node& rhs) const
        {
            return ((*this) < rhs || rhs < (*this));
        }
    int level;
    IntVect iv;
    int grid;
    typeEnum type;
};

std::ostream& operator<< (std::ostream&  os, const Node& node)
{
    os << "Node: IntVect=" << node.iv << ", level=" << node.level << ", grid=" << node.grid << ", type="; 
    if (node.type==Node::INIT)
        os << "INIT";
    else if (node.type==Node::COVERED)
        os << "COVERED";
    else
        os << "VALID";
    if (os.fail())
        amrex::Error("operator<<(std::ostream&,Node&) failed");
    return os;
}

struct Element
{
#if (BL_SPACEDIM==2)
#define MYLEN 4
    Element(const Node& a, const Node& b, const Node& c, const Node& d)
        {  n[0]=&a; n[1]=&b; n[2]=&c; n[3]=&d; }
#else
#define MYLEN 8
    Element(const Node& a, const Node& b, const Node& c, const Node& d,
            const Node& e, const Node& f, const Node& g, const Node& h)
        {  n[0]=&a; n[1]=&b; n[2]=&c; n[3]=&d; n[4]=&e; n[5]=&f; n[6]=&g; n[7]=&h; }
#endif
    const Node* n[MYLEN];
    inline bool operator< (const Element& rhs) const
        {
            for (int i=0; i<MYLEN; ++i)
                if (*n[i] != *rhs.n[i])
                    return *n[i] < *rhs.n[i];
            return false;
        }
};

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr 
	<< endl
	<< " usage:" << endl
	<< endl
	<< "    " << argv[0] << " [ParmParse input file] [additional keyword=input]" << endl
	<< endl
	<< "       The first argument is a file name unless it contains an =." << endl
	<< "       With no arguments prints this help." << endl
        << endl
	<< " keywords:" << endl
	<< endl
	<< "    box = int list of two subbox coords, LL and UR  (default all)" << endl
	<< "    comps = integer comp list  (overrides sComp and nComp)" << endl
	<< "    finestLevel = <#> finest level to use  (default all)" << endl
	<< "    help = <anything> prints this help" << endl
	<< "    infile = <plotfile>  (required)" << endl
	<< "    nComp = <#> number of comps  (default all)" << endl
	<< "    nGrowPer = <#> number of lev-0 cells by which to" << endl
	<< "               extend periodic boundaries  (default 0)" << endl
	<< "    sComp = start comp  (default 0)" << endl
	<< "    connect_cc = Generate flattened structure by connecting cells centers,"
        << "                 otherwise, generate node at all cell corners and copy cc"
        << "                 value out (default 1)" << endl
	<< endl
	<< " if nGrowPer > 1 then these additional keywords are needed:" << endl
	<< endl	
	<< "    geometry.coord_sys = <0 Cartesion, 1 rz>" << endl
	<< "    geometry.is_periodic = <0, false or 1, true for each axis>" << endl
	<< "    geometry.prob_lo = <lower limits for the axes>" << endl
	<< "    geometry.prob_hi = <upper limits for the axes>" << endl
	<< endl;
    amrex::Finalize();
    exit(1);
}

static
BoxArray
GetBndryCells (const BoxArray& ba,
               int             ngrow,
               const Geometry& geom)
{
    //
    // First get list of all ghost cells.
    //
    BoxList gcells, bcells;

    for (int i = 0; i < ba.size(); ++i)
	gcells.join(amrex::boxDiff(amrex::grow(ba[i],ngrow),ba[i]));
    //
    // Now strip out intersections with original BoxArray.
    //
    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(*it);

        if (isects.empty())
            bcells.push_back(*it);
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces;
            for (int i = 0; i < isects.size(); i++)
                pieces.push_back(isects[i].second);
            BoxList leftover = amrex::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }
    //
    // Now strip out overlaps.
    //
    gcells.clear();
    gcells = amrex::removeOverlap(bcells);
    bcells.clear();

    if (geom.isAnyPeriodic())
    {
        Vector<IntVect> pshifts(27);

        const Box& domain = geom.Domain();

        for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
        {
            if (!domain.contains(*it))
            {
                //
                // Add in periodic ghost cells shifted to valid region.
                //
                geom.periodicShift(domain, *it, pshifts);

                for (int i = 0; i < pshifts.size(); i++)
                {
                    const Box& shftbox = *it + pshifts[i];

                    const Box& ovlp = domain & shftbox;
                    BoxList bl = amrex::complementIn(ovlp,BoxList(ba));
                    bcells.catenate(bl);
                }
            }
        }

        gcells.catenate(bcells);
    }

    return BoxArray(gcells);
}

class FABdata
{
public:
    FABdata(int i, int n)
        {fab.resize(Box(IntVect::TheZeroVector(),
                        IntVect(AMREX_D_DECL(i-1,0,0))),n);}
    Real& operator[](int i) {return fab.dataPtr()[i];}
    Real* dataPtr (int N = 0) {return fab.dataPtr(N);}
    const Real* dataPtr (int N = 0) const {return fab.dataPtr(N);}
    FArrayBox fab;
    int boxSize;
};

void
Collate(Vector<Real>& NodeRaw,
        Vector<int>&  EltRaw,
        int          nCompPerNode)
{

#if BL_USE_MPI
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    BL_ASSERT(IOProc==0);
    const int nProcs = ParallelDescriptor::NProcs(); 
    Vector<int> nmdataR(nProcs,0);
    Vector<int> offsetR(nProcs,0);
    //
    // Tell root CPU how many Real data elements each CPU will be sending.
    //
    int countR = NodeRaw.size();
    MPI_Gather(&countR,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmdataR.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
            offsetR[i] = offsetR[i-1] + nmdataR[i-1];

        NodeRaw.resize(offsetR[nProcs-1] + nmdataR[nProcs-1]);
    }
    //
    // Gather all the Real data to IOProc into NodeRaw.
    //
    MPI_Gatherv(NodeRaw.dataPtr(),
                countR,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                NodeRaw.dataPtr(),
                nmdataR.dataPtr(),
                offsetR.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    // Now communicate the element info
    Vector<int> nmdataI(nProcs,0);
    Vector<int> offsetI(nProcs,0);
    int countI = EltRaw.size();
    MPI_Gather(&countI,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmdataI.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
            offsetI[i] = offsetI[i-1] + nmdataI[i-1];

        EltRaw.resize(offsetI[nProcs-1] + nmdataI[nProcs-1]);
    }
    //
    // Gather all the data to IOProc into EltRaw
    //
    MPI_Gatherv(EltRaw.dataPtr(),
                countI,
                ParallelDescriptor::Mpi_typemap<int>::type(),
                EltRaw.dataPtr(),
                nmdataI.dataPtr(),
                offsetI.dataPtr(),
                ParallelDescriptor::Mpi_typemap<int>::type(),
                IOProc,
                ParallelDescriptor::Communicator());
    //
    // Shift nodeIDs in element definitions
    //
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
        {
            const int nodeOffset = offsetR[i]/nCompPerNode;
            for (int j = 0; j < nmdataI[i]; j++)                
                EltRaw[offsetI[i]+j] += nodeOffset;
        }
    }
#endif
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool verbose = false;
    pp.query("verbose",verbose);

    if (verbose>1)
        AmrData::SetVerbose(true);

    std::string infile; pp.get("infile",infile);
    std::string outfile_DEF;

    std::string outType = "tec";
#ifdef USE_TEC_BIN_IO
    bool doBin = true;
    pp.query("doBin",doBin);
    outfile_DEF = infile+(doBin ? ".plt" : ".dat" );
#else
    bool doBin=false;
    outfile_DEF = infile+".dat";
#endif

    // 
    bool connect_cc = true; pp.query("connect_cc",connect_cc);

    std::string outfile(outfile_DEF); pp.query("outfile",outfile);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);

    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();

    const Vector<std::string>& names = amrData.PlotVarNames();

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = amrData.NComp();
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    Box subbox;
    if (int nx=pp.countval("box"))
    {
        Vector<int> barr;
        pp.getarr("box",barr,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(barr.size()==2*d);
        subbox=Box(IntVect(AMREX_D_DECL(barr[0],barr[1],barr[2])),
                   IntVect(AMREX_D_DECL(barr[d],barr[d+1],barr[d+2]))) & amrData.ProbDomain()[0];

    } else {

        subbox = amrData.ProbDomain()[0];
    }

    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Vector<BoxArray> gridArray(Nlev);
    Vector<Box> subboxArray(Nlev);

    int nGrowPer = 0; pp.query("nGrowPer",nGrowPer);
    PArray<Geometry> geom(Nlev);

    Vector<Real> LO(BL_SPACEDIM,0);
    Vector<Real> HI(BL_SPACEDIM,1);
    RealBox rb(LO.dataPtr(),HI.dataPtr());
    int coordSys = 0;
    Vector<int> isPer(BL_SPACEDIM,0);

    for (int lev=0; lev<Nlev; ++lev)
    {
        subboxArray[lev]
            = (lev==0 ? subbox
               : Box(subboxArray[lev-1]).refine(amrData.RefRatio()[lev-1]));

        geom.set(lev,new Geometry(amrData.ProbDomain()[lev],&rb,coordSys,const_cast<int*>(isPer.dataPtr())));

        if (nGrowPer>0 && lev==0)
        {
            for (int i=0; i<BL_SPACEDIM; ++i)
            {
                if (geom[lev].isPeriodic(i))
                {
                    if (subboxArray[lev].smallEnd()[i] == amrData.ProbDomain()[lev].smallEnd()[i])
                        subboxArray[lev].growLo(i,nGrowPer);
                    if (subboxArray[lev].bigEnd()[i] == amrData.ProbDomain()[lev].bigEnd()[i])
                        subboxArray[lev].growHi(i,nGrowPer);
                }
            }
        }

        gridArray[lev] = amrex::intersect(amrData.boxArray(lev), subboxArray[lev]);

        if (nGrowPer>0 && geom[lev].isAnyPeriodic() && gridArray[lev].size()>0)
        {
	  //const Box& probDomain = amrData.ProbDomain()[lev];
            const BoxArray& ba = amrData.boxArray(lev);
            BoxList bladd;
            Vector<IntVect> shifts;
            for (int i=0; i<ba.size(); ++i)
            {
                geom[lev].periodicShift(subboxArray[lev],ba[i],shifts);
                for (int j=0; j<shifts.size(); ++j)
                {
                    Box ovlp = Box(ba[i]).shift(shifts[j]) & subboxArray[lev];
                    if (ovlp.ok())
                        bladd.push_back(ovlp);
                }
            }
            bladd.simplify();
            BoxList blnew(gridArray[lev]);
            blnew.join(bladd);
            gridArray[lev] = BoxArray(blnew);
        }

        if (!gridArray[lev].size())
        {
            Nlev = lev;
            finestLevel = Nlev-1;
            gridArray.resize(Nlev);
            subboxArray.resize(Nlev);
        }
    }

    const int nGrow = 1;
    typedef BaseFab<Node> NodeFab;
    typedef FabArray<NodeFab> MultiNodeFab; 
    PArray<MultiNodeFab> nodes(Nlev,PArrayManage);

    std::cerr << "Before nodes allocated" << endl;
    for (int lev=0; lev<Nlev; ++lev)
        nodes.set(lev,new MultiNodeFab(gridArray[lev],1,nGrow));
    std::cerr << "After nodes allocated" << endl;

    int cnt = 0;
    typedef std::map<Node,int> NodeMap;
    NodeMap nodeMap;
    for (int lev=0; lev<Nlev; ++lev)
    {
        for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodes[lev][fai];
            const Box& box = ifab.box() & subboxArray[lev];
            for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv))
                ifab(iv,0) = Node(iv,lev,fai.index(),Node::VALID);
        }
            
        if (lev != 0)
        {
            const int ref = amrData.RefRatio()[lev-1];
            const Box& rangeBox = Box(IntVect::TheZeroVector(),
                                     (ref-1)*IntVect::TheUnitVector());

            BoxArray bndryCells = GetBndryCells(nodes[lev].boxArray(),ref,geom[lev]);

            for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
            {
                const Box& box = amrex::grow(fai.validbox(),ref) & subboxArray[lev];
                NodeFab& ifab = nodes[lev][fai];
                std::vector< std::pair<int,Box> > isects = bndryCells.intersections(box);
                for (int i = 0; i < isects.size(); i++)
                {
                    Box co = isects[i].second & fai.validbox();
                    if (co.ok())
                        std::cout << "BAD ISECTS: " << co << std::endl;

                    const Box& dstBox = isects[i].second;
                    const Box& srcBox = amrex::coarsen(dstBox,ref);

                    NodeFab dst(dstBox,1);
                    for (IntVect iv(srcBox.smallEnd());
                         iv<=srcBox.bigEnd();
                         srcBox.next(iv))
                    {
                        const IntVect& baseIV = ref*iv;
                        for (IntVect ivt(rangeBox.smallEnd());ivt<=rangeBox.bigEnd();rangeBox.next(ivt))
                            dst(baseIV + ivt,0) = Node(iv,lev-1,-1,Node::VALID);
                    }
                    const Box& ovlp = dstBox & ifab.box();

                    Box mo = ovlp & fai.validbox();
                    if (mo.ok())
                    {
                        std::cout << "BAD OVERLAP: " << mo << std::endl;
                        std::cout << "         vb: " << fai.validbox() << std::endl;
                    }
                    if (ovlp.ok())
                        ifab.copy(dst,ovlp,0,ovlp,0,1);
                }
            }
        }

        // Block out cells covered by finer grid
        if (lev < finestLevel)
        {
            const BoxArray coarsenedFineBoxes =
                BoxArray(gridArray[lev+1]).coarsen(amrData.RefRatio()[lev]);
                
            for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
            {
                NodeFab& ifab = nodes[lev][fai];
                const Box& box = ifab.box();
                std::vector< std::pair<int,Box> > isects = coarsenedFineBoxes.intersections(box);
                for (int i = 0; i < isects.size(); i++)
                {
                    const Box& ovlp = isects[i].second;
                    for (IntVect iv=ovlp.smallEnd(); iv<=ovlp.bigEnd(); ovlp.next(iv))
                        ifab(iv,0) = Node(iv,lev,fai.index(),Node::COVERED);
                }
            }
        }

        // Add the unique nodes from this level to the list
        for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodes[lev][fai];
            const Box& box = fai.validbox() & subboxArray[lev];
            for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
            {
                if (ifab(iv,0).type == Node::VALID)
                {
                    if (ifab(iv,0).level != lev) 
                        std::cout << "bad level: " << ifab(iv,0) << std::endl;
                    nodeMap[ifab(iv,0)] = cnt++;
                }
            }
        }
    }

    std::cerr << "After nodeMap built, size=" << nodeMap.size() << endl;

    typedef std::set<Element> EltSet;
    EltSet elements;
    for (int lev=0; lev<Nlev; ++lev)
    {
        for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodes[lev][fai];                        
            Box box = ifab.box() & subboxArray[lev];

            for (int dir=0; dir<BL_SPACEDIM; ++dir)
                box.growHi(dir,-1);

            for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
            {
#if (BL_SPACEDIM == 2)
                const Node& n1 = ifab(iv,0);
                const Node& n2 = ifab(IntVect(iv).shift(amrex::BASISV(0)),0);
                const Node& n3 = ifab(IntVect(iv).shift(IntVect::TheUnitVector()),0);
                const Node& n4 = ifab(IntVect(iv).shift(amrex::BASISV(1)),0);
                if (n1.type==Node::VALID && n2.type==Node::VALID &&
                    n3.type==Node::VALID && n4.type==Node::VALID )
                    elements.insert(Element(n1,n2,n3,n4));
#else
                const IntVect& ivu = IntVect(iv).shift(amrex::BASISV(2));
                const Node& n1 = ifab(iv ,0);
                const Node& n2 = ifab(IntVect(iv ).shift(amrex::BASISV(0)),0);
                const Node& n3 = ifab(IntVect(iv ).shift(amrex::BASISV(0)).shift(amrex::BASISV(1)),0);
                const Node& n4 = ifab(IntVect(iv ).shift(amrex::BASISV(1)),0);
                const Node& n5 = ifab(ivu,0);
                const Node& n6 = ifab(IntVect(ivu).shift(amrex::BASISV(0)),0);
                const Node& n7 = ifab(IntVect(ivu).shift(amrex::BASISV(0)).shift(amrex::BASISV(1)),0);
                const Node& n8 = ifab(IntVect(ivu).shift(amrex::BASISV(1)),0);
                if (n1.type==Node::VALID && n2.type==Node::VALID &&
                    n3.type==Node::VALID && n4.type==Node::VALID &&
                    n5.type==Node::VALID && n6.type==Node::VALID &&
                    n7.type==Node::VALID && n8.type==Node::VALID )
                    elements.insert(Element(n1,n2,n3,n4,n5,n6,n7,n8));
#endif
            }
        }
    }

    int nElts = (connect_cc ? elements.size() : nodeMap.size());
    std::cerr << "Before connData allocated " << elements.size() << " elements"  << endl;
    Vector<int> connData(MYLEN*nElts);
    std::cerr << "After connData allocated " << elements.size() << " elements" << endl;

    if (connect_cc) {
        cnt = 0;
        for (EltSet::const_iterator it = elements.begin(); it!=elements.end(); ++it)
        {
            for (int j=0; j<MYLEN; ++j)
            {
                const NodeMap::const_iterator noit = nodeMap.find(*((*it).n[j]));
                if (noit == nodeMap.end())
                {
                    std::cout << "Node not found in node map" << std::endl;
                    std::cout << *((*it).n[j]) << std::endl;
                }
                else
                {
                    connData[cnt++] = noit->second+1;
                }
            }
        }
    } else {
        cnt = 1;
        for (int i=0; i<nElts; ++i) {
            for (int j=0; j<MYLEN; ++j) {
                connData[i*MYLEN+j] = cnt++;
            }
        }
    }

    std::cerr << "Final elements built" << endl;

    // Invert the map
    std::vector<Node> nodeVect(nodeMap.size());
    for (NodeMap::const_iterator it=nodeMap.begin(); it!=nodeMap.end(); ++it)
    {
        if (it->second>=nodeVect.size() || it->second<0)
            std::cout << "Bad id: " << it->second << "  bad node: " << it->first << std::endl;
        BL_ASSERT(it->second>=0);
        BL_ASSERT(it->second<nodeVect.size());
        nodeVect[it->second] = (*it).first;
    }
    std::cerr << "Final nodeVect built (" << nodeVect.size() << " nodes)" << endl;
		
    nodeMap.clear();
    elements.clear();
    nodes.clear();

    std::cerr << "Temp nodes, elements cleared" << endl;

    PArray<MultiFab> fileData(Nlev);
    int ng = nGrowPer;
    for (int lev=0; lev<Nlev; ++lev)
    {
        if (lev!=0)
            ng *= amrData.RefRatio()[lev-1];
        const BoxArray& ba = gridArray[lev];
        fileData.set(lev,new MultiFab(ba,comps.size(),0));
        fileData[lev].setVal(1.e30);
        std::cerr << "My data set alloc'd at lev=" << lev << endl;

        MultiFab pData, pDataNG;
        if (ng>0 && geom[lev].isAnyPeriodic())
        {
            const Box& pd = amrData.ProbDomain()[lev];
            //const BoxArray& vba = amrData.boxArray(lev);
            Box shrunkenDomain = pd;
            for (int i=0; i<BL_SPACEDIM; ++i)
                if (geom[lev].isPeriodic(i))
                    shrunkenDomain.grow(i,-ng);

            const BoxArray edgeVBoxes = amrex::boxComplement(pd,shrunkenDomain);
            pData.define(edgeVBoxes,1,ng,Fab_allocate);
            pDataNG.define(BoxArray(edgeVBoxes).grow(ng),1,0,Fab_allocate);
        }

        for (int i=0; i<comps.size(); ++i)
        {
            BoxArray tmpBA = amrex::intersect(fileData[lev].boxArray(),amrData.ProbDomain()[lev]);
            MultiFab tmpMF(tmpBA,1,0);
            tmpMF.setVal(2.e30);
            amrData.FillVar(tmpMF,lev,names[comps[i]],0);
            fileData[lev].copy(tmpMF,0,i,1);
            if (ng>0 && geom[lev].isAnyPeriodic())
            {
                pData.setVal(3.e30);
                pDataNG.copy(tmpMF);
                for (MFIter mfi(pData); mfi.isValid(); ++mfi)
                    pData[mfi].copy(pDataNG[mfi]);
                amrData.FillVar(pData,lev,names[comps[i]],0);
		pData.EnforcePeriodicity(geom[lev]);
                for (MFIter mfi(pData); mfi.isValid(); ++mfi)
                    pDataNG[mfi].copy(pData[mfi]);
                fileData[lev].copy(pDataNG,0,i,1);
            }            
        }

        if (fileData[lev].max(0) > 1.e29)
        {
            std::cerr << "Bad mf data" << std::endl;
            VisMF::Write(fileData[lev],"out.mfab");
            amrex::Abort();
        }
    }
    std::cerr << "File data loaded" << endl;

    int nNodesFINAL = (connect_cc  ?  nodeVect.size() : nElts*MYLEN );
    int nCompsFINAL = BL_SPACEDIM+comps.size();
    FABdata tmpData(nNodesFINAL,nCompsFINAL);
    int tmpDatLen = nCompsFINAL*nNodesFINAL;
    std::cerr << "Final node data allocated (size=" << tmpDatLen << ")" << endl;

    //const Vector<Vector<Real> >& dxLevel = amrData.DxLevel();  // Do not trust dx from file...compute our own
    Vector<Vector<Real> > dxLevel(Nlev);
    for (int i=0; i<Nlev; ++i)
    {
        dxLevel[i].resize(BL_SPACEDIM);
        for (int j=0; j<BL_SPACEDIM; ++j)
            dxLevel[i][j] = amrData.ProbSize()[j]/amrData.ProbDomain()[i].length(j);
    }
    const Vector<Real>& plo = amrData.ProbLo();

    // Handy structures for loading data in usual fab ordering (transpose of mef/flt ordering)
#define BIN_POINT
#undef BIN_POINT
#ifdef BIN_POINT
    Real* data = tmpData.dataPtr();
#else /* BLOCK ordering */
    Vector<Real*> fdat(comps.size()+BL_SPACEDIM);
    for (int i=0; i<fdat.size(); ++i)
        fdat[i] = tmpData.dataPtr(i);
#endif

    cnt = 0;
    int Nnodes = nodeVect.size();
    int jGridPrev = 0;
    int levPrev = -1;
    for (int i=0; i<Nnodes; ++i)
    {
        const Node& node = nodeVect[i];

        const Vector<Real>& dx = dxLevel[node.level];
        const IntVect& iv = node.iv;

        const BoxArray& grids = fileData[node.level].boxArray();
        int jGrid = node.grid;
        if (jGrid<0)
        {
            bool found_it = false;
            
            // Try same grid as last time, otherwise search list
            if (node.level==levPrev && grids[jGridPrev].contains(iv))
            {
                found_it = true;
                jGrid = jGridPrev;
            }
            for (int j=0; j<grids.size() && !found_it; ++j)
            {
                if (grids[j].contains(iv))
                {
                    found_it = true;
                    jGrid = j;
                }
            }
            BL_ASSERT(found_it);
        }
	// Remember these for next time
	levPrev = node.level;
	jGridPrev = jGrid;

	Vector<IntVect> ivt;
        if (connect_cc) {
            ivt.resize(1,iv);
        } else {
            ivt.resize(AMREX_D_PICK(1,4,8),iv);
            ivt[1] += amrex::BASISV(0);
            ivt[2] = ivt[1] + amrex::BASISV(1);
            ivt[3] += amrex::BASISV(1);
#if BL_SPACEDIM==3
            for (int n=0; n<4; ++n) {
                ivt[4+n] = iv[n] + amrex::BASISV(2);
            }
#endif
        }

	for (int j=0; j<ivt.size(); ++j) {

            Real offset = (connect_cc ? 0.5 : 0);

#ifdef BIN_POINT
	  for (int dir=0; dir<BL_SPACEDIM; ++dir)
            data[cnt++] = plo[dir] + (ivt[j][dir] + offset) * dx[dir];
#else /* BLOCK ordering */
	  for (int dir=0; dir<BL_SPACEDIM; ++dir)
            fdat[dir][cnt] = plo[dir] + (ivt[j][dir] + offset) * dx[dir];
#endif /* BIN_POINT */


#ifdef BIN_POINT
	  for (int n=0; n<comps.size(); ++n) {
            data[cnt++] = fileData[node.level][jGrid](iv,n);
          }
#else /* BLOCK ordering */
	  for (int n=0; n<comps.size(); ++n) {
              fdat[n+BL_SPACEDIM][cnt] = fileData[node.level][jGrid](iv,n);
          }
	  cnt++;
#endif /* BIN_POINT */
	}
    }

    //Collate(tmpData,connData,MYLEN);

    //
    // Write output
    //
    const int nState = BL_SPACEDIM + comps.size();
    std::string vars = AMREX_D_TERM("X"," Y"," Z");
    for (int j=0; j<comps.size(); ++j)
        vars += " " + amrData.PlotVarNames()[comps[j]];

    if (outType == "tec")
    {

#ifdef BIN_POINT
        string block_or_point = "FEPOINT";
#else /* BLOCK ordering */
        string block_or_point = "FEBLOCK";
#endif


        if (doBin)
        {
#ifdef USE_TEC_BIN_IO
            INTEGER4 Debug = 0;
            INTEGER4 VIsDouble = 1;
            INTEGER4 EltID = AMREX_D_PICK(0,1,3);
            TECINI((char*)"Pltfile data", (char*)vars.c_str(), (char*)outfile.c_str(), (char*)".", &Debug, &VIsDouble);
            INTEGER4 nPts = nNodesFINAL;
            TECZNE((char*)infile.c_str(), &nPts, &nElts, &EltID, (char*)block_or_point.c_str(), NULL);
            TECDAT(&tmpDatLen,tmpData.fab.dataPtr(),&VIsDouble);
            TECNOD(connData.dataPtr());
            TECEND();
#else
            amrex::Abort("Need to recompile with USE_TEC_BIN_IO defined");
#endif
        }
        else
        {
            std::ofstream osf(outfile.c_str(),std::ios::out);
            osf << AMREX_D_TERM("VARIABLES= \"X\"", " \"Y\"", " \"Z\"");
            for (int j=0; j<comps.size(); ++j)
                osf << " \""  << amrData.PlotVarNames()[comps[j]] << "\"";
            char buf[100];
            sprintf(buf,"%g",amrData.Time());
            osf << endl << "ZONE T=\"" << infile << " time = " << buf
                << "\", N=" << nNodesFINAL << ", E=" << nElts << ", F=" << "FEPOINT" << " ET="
                //<< "\", N=" << nPts << ", E=" << nElts << ", F=" << block_or_point << " ET="
                << AMREX_D_PICK("POINT","QUADRILATERAL","BRICK") << endl;
            
            for (int i=0; i<nNodesFINAL; ++i)
            {
                for (int j=0; j<nState; ++j)
                    osf << tmpData.dataPtr(j)[i] << " ";
                osf << endl;
            }
            for (int i=0; i<nElts; ++i)
            {
                for (int j=0; j<MYLEN; ++j)
                    osf << connData[i*MYLEN+j] << " ";
                osf << endl;
            }
            osf << endl;
            osf.close();
        }
    }
    else
    {
        std::ofstream ofs;
        std::ostream* os =
            (outfile=="-" ? (std::ostream*)(&std::cout) : (std::ostream*)(&ofs) );
        if (outfile!="-")
            ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
        (*os) << infile << " time = " << amrData.Time() << endl;
        (*os) << vars << endl;
        (*os) << nElts << " " << MYLEN << endl;
        tmpData.fab.writeOn(*os);
        (*os).write((char*)connData.dataPtr(),sizeof(int)*connData.size());
        if (outfile!="-")
            ofs.close();
    }

    amrex::Finalize();
    return 0;
}
