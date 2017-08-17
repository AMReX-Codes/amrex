#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include "kdtree_F.H"

using namespace amrex;

struct KDNode {
    
    Box box;
    KDNode* left;
    KDNode* right;
    Real cost;
    int num_procs_left;    

    KDNode(const Box& box_in, Real cost_in, int num_procs_in)
        : box(box_in), cost(cost_in), num_procs_left(num_procs_in)
    {
        left = NULL;
        right = NULL;
    }
};

class KDTree {
    
public:
    
    KDTree(const Box& domain, const FArrayBox& cost, int num_procs) {
        Real total_cost = cost.sum(0);
        root = new KDNode(domain, total_cost, num_procs);
        buildKDTree(root, cost);
    }
    
    ~KDTree() {
        freeKDTree(root);
    }
    
    void GetBoxes(BoxList& bl, Array<Real>& costs) {        
        walkKDTree(root, bl, costs);
    }
    
private:
    
    void buildKDTree(KDNode* node, const FArrayBox& cost) {
        if (node->num_procs_left == 1) return;
        
        partitionNode(node, cost);
        
        buildKDTree(node->left,  cost);
        buildKDTree(node->right, cost);
    }

    void freeKDTree(KDNode* node) {
        
        if (node != NULL) {
            
            KDNode* left  = node->left;
            KDNode* right = node->right;
            
            delete node;
            
            freeKDTree(left);
            freeKDTree(right);
        }
    }

    void walkKDTree(KDNode* node, BoxList& bl, Array<Real>& costs) {
        
        if (node->left == NULL && node->right == NULL) {
            costs.push_back(node->cost);
            bl.push_back(node->box);
            return;
        }
        
        walkKDTree(node->left,  bl, costs);
        walkKDTree(node->right, bl, costs);
    }
    
    void partitionNode(KDNode* node, const FArrayBox& cost) {
        
        const Box& box = node->box;
        BL_ASSERT(cost.box().contains(box));
        
        int split;
        Real cost_left, cost_right;
        int dir = getLongestDir(box);
        compute_best_partition(cost.dataPtr(), cost.loVect(), cost.hiVect(),
                               box.loVect(), box.hiVect(), node->cost, dir,
                               &cost_left, &cost_right, &split);
        
        Box left, right;
        splitBox(split, dir, box, left, right);
        
        node->left  = new KDNode(left,  cost_left,  node->num_procs_left/2);
        node->right = new KDNode(right, cost_right, node->num_procs_left/2);
    }

    int getLongestDir(const Box& box) {
        IntVect size = box.size();
        int argmax = 0;
        int max = size[0];
        for (int i = 1; i < BL_SPACEDIM; ++i) {
            if (size[i] > max) {
                max = size[i];
                argmax = i;
            }
        }
        return argmax;
    }
    
    void splitBox(int split, int dir,
                  const Box& box, Box& left, Box& right) {
        left = box;
        right = box;
        
        left.setBig(dir, split);
        right.setSmall(dir, split+1);        
    }

    KDNode* root;
    
};

class MyParticleContainer 
    : public ParticleContainer<0>
{
 public:
    
    MyParticleContainer (const Geometry& geom, 
                         const DistributionMapping& dmap,
                         const BoxArray& ba)
        : ParticleContainer<0> (geom, dmap, ba) {}
};

void computeCost(MultiFab& local_cost, 
                 MultiFab& global_cost,
                 MyParticleContainer& myPC) {
    
    const int lev = 0;
    const Geometry& geom = myPC.Geom(lev);
    const Box& domain = geom.Domain();
    const BoxArray& ba = myPC.ParticleBoxArray(lev);
    const DistributionMapping& dm = myPC.ParticleDistributionMap(lev);

    BoxList global_bl;
    Array<int> procs_map;
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
        global_bl.push_back(domain);
        procs_map.push_back(i);
    }

    BoxArray global_ba(global_bl);    
    DistributionMapping global_dm(procs_map);    
    
    local_cost.define(ba, dm, 1, 0);
    local_cost.setVal(0.0);
    myPC.Increment(local_cost, lev);

    global_cost.define(global_ba, global_dm, 1, 0);
    global_cost.copy(local_cost, 0, 0, 1);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int num_cells, max_grid_size, num_procs;
   
    ParmParse pp;    
    pp.get("num_cells", num_cells);
    pp.get("num_procs", num_procs);
    pp.get("max_grid_size", max_grid_size);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(num_cells - 1, num_cells - 1, num_cells - 1));
    const Box domain(domain_lo, domain_hi);
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0; 
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dmap(ba);

    MyParticleContainer myPC(geom, dmap, ba);    
    MyParticleContainer::ParticleInitData pdata = {};
    myPC.InitFromBinaryFile("binary_particle_file.dat", 0);

    MultiFab local_cost;
    MultiFab global_cost;    
    computeCost(local_cost, global_cost, myPC);

    FArrayBox *cost;
    for ( MFIter mfi(global_cost); mfi.isValid(); ++mfi ) {
        cost = &global_cost[mfi];
    }

    KDTree tree = KDTree(domain, *cost, num_procs);

    BoxList new_bl;
    Array<Real> box_costs;
    tree.GetBoxes(new_bl, box_costs);
    BoxArray new_ba(new_bl);
    
    for (int i = 0; i < box_costs.size(); ++i) {
        std::cout << box_costs[i] << " " << new_ba[i] << std::endl;
    }
    
    Array<int> new_pmap;
    for (int i = 0; i < new_ba.size(); ++i) {
        new_pmap.push_back(0);
    }

    DistributionMapping new_dm(new_pmap);

    myPC.SetParticleBoxArray(0, new_ba);
    myPC.SetParticleDistributionMap(0, new_dm);
    
    myPC.Redistribute();
    MultiFab new_local_cost;
    MultiFab new_global_cost;    
    computeCost(new_local_cost, new_global_cost, myPC);

    WriteSingleLevelPlotfile("plt00000", new_local_cost, {"cost"},
                             geom, 0.0, 0);
    myPC.Checkpoint("plt00000", "particle0", true);
    
    amrex::Finalize();
}
