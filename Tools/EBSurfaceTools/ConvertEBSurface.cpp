#include <AMReX.H>
#include <AMReX_GeomIntersectUtils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VectorIO.H>
#include <AMReX_Utility.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  ParmParse pp;

  std::string infile;
  pp.get("infile",infile);

  std::string FullHeaderPath = infile;
  if (!FullHeaderPath.empty() && FullHeaderPath[FullHeaderPath.size()-1] != '/')
    FullHeaderPath += '/';
  FullHeaderPath += "Header";

  if (ParallelDescriptor::IOProcessor())
  {
    std::ifstream ifs(FullHeaderPath.c_str());
    std::string versionStr;
    ifs >> versionStr;

    int nDataComp;
    ifs >> nDataComp;

    int nNodesPerElt;
    ifs >> nNodesPerElt;

    int nFiles;
    ifs >> nFiles;

    Vector<long> fileNum(nFiles);
    Vector<long> nNodes(nFiles);
    Vector<long> nElts(nFiles);
    for (int i=0; i<nFiles; ++i)
    {
      ifs >> fileNum[i];
    }
    for (int i=0; i<nFiles; ++i)
    {
      ifs >> nNodes[i];
    }
    for (int i=0; i<nFiles; ++i)
    {
      ifs >> nElts[i];
    }

    IntDescriptor id;
    ifs >> id;
    RealDescriptor rd;
    ifs >> rd;
    IntDescriptor ld;
    ifs >> ld;

    int SizeOfEdgeData = 2 * AMREX_SPACEDIM;
    int SizeOfNodeData = nDataComp;
    int SizeOfFragData = nNodesPerElt;

    NodeMap nodeMap;
    Vector<NodeMapIt> sortedNodes;
    Vector<Triangle> triangles;

    int nNodesRedundant = 0;
    for (int i=0; i<nFiles; ++i)
    {
      auto FullDataPath = infile;
      if (!FullDataPath.empty() && FullDataPath[FullDataPath.size()-1] != '/')
        FullDataPath += '/';
      FullDataPath += "Data";

      auto nDigits = std::log10(fileNum[i]) + 1;
      FullDataPath += Concatenate("_",fileNum[i],nDigits);

      std::ifstream ifsd(FullDataPath.c_str(), std::ios::binary);
      if( ! ifsd.good()) { amrex::FileOpenFailed(FullDataPath); }

        
      Vector<int>  flattenedEdges(SizeOfEdgeData * nNodes[i]);
      Vector<Real> flattenedNodes(SizeOfNodeData * nNodes[i]);
      Vector<long> flattenedFrags(SizeOfFragData * nElts[i]);

      readIntData( flattenedEdges.dataPtr(), flattenedEdges.size(), ifsd, id);
      readRealData(flattenedNodes.dataPtr(), flattenedNodes.size(), ifsd, rd);
      readLongData(flattenedFrags.dataPtr(), flattenedFrags.size(), ifsd, ld);

      ifsd.close();

      IntVect L,R;
      Vector<Real> data(SizeOfNodeData);
      Vector<int> globalNum(nNodes[i]);

      for (int j=0; j<nNodes[i]; ++j)
      {
        auto eoffset = j*SizeOfEdgeData;
        auto noffset = j*SizeOfNodeData;
        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
          L[d] = flattenedEdges[eoffset+d               ];
          R[d] = flattenedEdges[eoffset+d+AMREX_SPACEDIM];
        }
        for (int d=0; d<SizeOfNodeData; ++d)
        {
          data[d] = flattenedNodes[noffset+d];
        }
        Edge e(L,R);
        e.ID = j;
        auto p = make_pair(e,data);
        AMREX_ASSERT(ret.second);

        auto res = nodeMap.find(e);
        if (res == nodeMap.end()) {
          p.first.ID = nodeMap.size();
          globalNum[j] = p.first.ID;
          res = nodeMap.insert(p).first;
          sortedNodes.push_back(res);
        }
        else {
          globalNum[j] = res->first.ID;
        }
      }

      for (int j=0; j<nElts[i]; ++j)
      {
        auto foffset = j*SizeOfFragData;
        triangles.push_back(Triangle());
        auto& frag = triangles.back();
        AMREX_ASSERT(frag.size() == SizeOfFragData);
        for (int d=0; d<SizeOfFragData; ++d)
        {
          frag[d] = sortedNodes[globalNum[flattenedFrags[foffset+d]]];
        }
      }
      nNodesRedundant += nNodes[i];
    }

    Print() << sortedNodes.size() << " unique nodes and " << triangles.size() << " triangles found";
    if (nNodesRedundant != sortedNodes.size())
    {
      Print() << " (" << nNodesRedundant - sortedNodes.size() << " redundant nodes removed)"; 
    }
    Print() << '\n';

    std::string outfile;
    pp.get("outfile",outfile);
    std::ofstream ofs(outfile.c_str());
    ofs << sortedNodes.size() << " " << triangles.size() << std::endl;

    for (int j=0; j<sortedNodes.size(); ++j)
    {
      const auto& vec = sortedNodes[j]->second;
      ofs << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
    }
    for (int j=0; j<triangles.size(); ++j)
    {
      const Triangle& tri = triangles[j];
      for (int k=0; k<3; ++k)
      {
        ofs << tri[k]->first.ID + 1 << " ";
      }
      ofs << std::endl;
    }
    ofs.close();

  }

  amrex::Finalize();

  return 0;
}

