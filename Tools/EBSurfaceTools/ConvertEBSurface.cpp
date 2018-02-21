#include <AMReX.H>
#include <AMReX_GeomIntersectUtils.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VectorIO.H>
#include <AMReX_Utility.H>

using namespace amrex;
using std::list;

list<list<Segment>> MakePolyLines(Vector<Segment>& segVec);

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
#if AMREX_SPACEDIM==2
    Vector<Segment> surfaceFragmentsG;
#elif AMREX_SPACEDIM==3
    Vector<Triangle> surfaceFragmentsG;
#else
    amrex::Abort("This code not applicable to 1D");
#endif

    auto fileNumberMAX = fileNum[0];
    for (int i=1; i<nFiles; ++i)
    {
      fileNumberMAX = std::max(fileNumberMAX,fileNum[i]);
    }
    int nDigits = std::log10(fileNumberMAX) + 1;

    int nNodesRedundant = 0;
    for (int i=0; i<nFiles; ++i)
    {
      auto FullDataPath = infile;
      if (!FullDataPath.empty() && FullDataPath[FullDataPath.size()-1] != '/')
        FullDataPath += '/';
      FullDataPath += "Data";

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
#if AMREX_SPACEDIM==2
        surfaceFragmentsG.push_back(Segment());
#elif AMREX_SPACEDIM==3
        surfaceFragmentsG.push_back(Triangle());
#endif
        auto& frag = surfaceFragmentsG.back();
        AMREX_ASSERT(frag.size() == SizeOfFragData);
        for (int d=0; d<SizeOfFragData; ++d)
        {
          frag[d] = sortedNodes[globalNum[flattenedFrags[foffset+d]]];
        }
      }
      nNodesRedundant += nNodes[i];
    }

    Print() << sortedNodes.size() << " unique nodes and " << surfaceFragmentsG.size() << " surfaceFragmentsG found";
    if (nNodesRedundant != sortedNodes.size())
    {
      Print() << " (" << nNodesRedundant - sortedNodes.size() << " redundant nodes removed)"; 
    }
    Print() << '\n';

    std::string outfile;
    pp.get("outfile",outfile);
    std::ofstream ofs(outfile.c_str());

#if AMREX_SPACEDIM==2
    list<list<Segment>> contours = MakePolyLines(surfaceFragmentsG);
    if (contours.size() > 0)
    {
      ofs << contours.size() << '\n';
      for (std::list<list<Segment>>::iterator it = contours.begin(); it!=contours.end(); ++it)
      {
        ofs << it->size() + 1 << '\n';
        auto slit=it->begin();
        const auto& pt0 = (*slit)[0]->second;
        ofs << pt0[0] << " " << pt0[1] << '\n';
        for (  ; slit!=it->end(); ++slit)
        {
          const auto& pt = (*slit)[1]->second;
          ofs << pt[0] << " " << pt[1] << '\n';
        }
      }
    }
#else
    ofs << sortedNodes.size() << " " << surfaceFragmentsG.size() << std::endl;

    for (int j=0; j<sortedNodes.size(); ++j)
    {
      const auto& vec = sortedNodes[j]->second;
      ofs << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
    }
    for (int j=0; j<surfaceFragmentsG.size(); ++j)
    {
      const auto& tri = surfaceFragmentsG[j];
      for (int k=0; k<tri.size(); ++k)
      {
        ofs << tri[k]->first.ID + 1 << " ";
      }
      ofs << std::endl;
    }
#endif

    ofs.close();

  }

  amrex::Finalize();

  return 0;
}

static
list<Segment>::iterator FindMySeg(list<Segment>& segs, const NodeMapIt& idx)
{
  for (list<Segment>::iterator it=segs.begin(); it!=segs.end(); ++it)
  {
    AMREX_ASSERT(it->size() == 2);
    if ( ((*it)[0] == idx) || ((*it)[1] == idx) )
      return it;
  }
  return segs.end();
}

list<list<Segment>> MakePolyLines(Vector<Segment>& segVec)
{
  list<Segment> segments;
  for (int i=0; i<segVec.size(); ++i) segments.push_back(segVec[i]);


  list<list<Segment>> contourLines;
  contourLines.push_back(list<Segment>());

  if (segments.size() > 0)
  {
    contourLines.push_back(list<Segment>());
    contourLines.back().push_back(segments.front());
    segments.pop_front();

    auto idx = contourLines.back().back()[1];
    while (segments.begin() != segments.end())
    {
      auto segIt = FindMySeg(segments,idx);
      if (segIt != segments.end())
      {
        const auto& idx_l = (*segIt)[0];
        const auto& idx_r = (*segIt)[1];

        if ( idx_l == idx )
        {
          idx = idx_r;
          contourLines.back().push_back(*segIt);
        }
        else
        {
          idx = idx_l;
          contourLines.back().push_back(Segment(idx_r,idx_l));
        }

        segments.erase(segIt);
      }
      else
      {
        contourLines.push_back(list<Segment>());
        contourLines.back().push_back(segments.front());
        segments.pop_front();
        
        idx = contourLines.back().back()[1];
      }
    }
  }

  // Connect up the line segments as much as possible
  bool changed;
  do
  {
    changed = false;
    for (std::list<list<Segment>>::iterator it = contourLines.begin(); it!=contourLines.end(); ++it)
    {
      if (!it->empty())
      {
        const auto& idx_l = it->front()[0];
        const auto& idx_r = it->back()[1];
        for (std::list<list<Segment>>::iterator it1 = contourLines.begin(); it1!=contourLines.end(); ++it1)
        {
          if (!it1->empty() && it!=it1)
          {
            if (idx_r == it1->front()[0])
            {
              it->splice(it->end(),*it1);
              changed = true;
            }
            else if (idx_r == it1->back()[1])
            {
              it1->reverse();
              for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
              {
                const auto tmp = (*it2)[0];
                (*it2)[0] = (*it2)[1];
                (*it2)[1] = tmp;
              }
              it->splice(it->end(),*it1);
              changed = true;
            }
            else if (idx_l == it1->front()[0])
            {
              it1->reverse();
              for (list<Segment>::iterator it2=it1->begin(); it2!=it1->end(); ++it2)
              {
                const auto tmp = (*it2)[0];
                (*it2)[0] = (*it2)[1];
                (*it2)[1] = tmp;
              }
              it->splice(it->begin(),*it1);
              changed = true;
            }
          }
        }
      }
    }
  } while(changed);

  // Clear out empty placeholders for lines we connected up to others.
  for (std::list<list<Segment>>::iterator it = contourLines.begin(); it!=contourLines.end();)
  {
    if (it->empty())
      contourLines.erase(it++);
    else
      it++;
  }

  return contourLines;
}
