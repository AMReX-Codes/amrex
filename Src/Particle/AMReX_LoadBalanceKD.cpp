#include "AMReX_LoadBalanceKD.H"

#include "AMReX_KDTree_F.H"

using namespace amrex;

KDTree::KDTree(const Box& domain, const FArrayBox& cost, int num_procs) {
    Real total_cost = cost.sum(0);
    root = new KDNode(domain, total_cost, num_procs);
    buildKDTree(root, cost);
}

KDTree::~KDTree() {
    freeKDTree(root);
}
    
void KDTree::GetBoxes(BoxList& bl, Array<Real>& costs) {        
    walkKDTree(root, bl, costs);
}
    
void KDTree::buildKDTree(KDNode* node, const FArrayBox& cost) {
    if (node->num_procs_left == 1) return;
    
    partitionNode(node, cost);
    
    buildKDTree(node->left,  cost);
    buildKDTree(node->right, cost);
}

void KDTree::freeKDTree(KDNode* node) {
    
    if (node != NULL) {
        
        KDNode* left  = node->left;
        KDNode* right = node->right;
        
        delete node;
        
        freeKDTree(left);
        freeKDTree(right);
    }
}

void KDTree::walkKDTree(KDNode* node, BoxList& bl, Array<Real>& costs) {
    
    if (node->left == NULL && node->right == NULL) {
        costs.push_back(node->cost);
        bl.push_back(node->box);
        return;
    }
    
    walkKDTree(node->left,  bl, costs);
    walkKDTree(node->right, bl, costs);
}

void KDTree::partitionNode(KDNode* node, const FArrayBox& cost) {
    
    const Box& box = node->box;
    BL_ASSERT(cost.box().contains(box));
    
    int split;
    Real cost_left, cost_right;
    int dir = getLongestDir(box);
    amrex_compute_best_partition(cost.dataPtr(), cost.loVect(), cost.hiVect(),
                                 box.loVect(), box.hiVect(), node->cost, dir,
                                 &cost_left, &cost_right, &split);
    
    Box left, right;
    splitBox(split, dir, box, left, right);
    
    node->left  = new KDNode(left,  cost_left,  node->num_procs_left/2);
    node->right = new KDNode(right, cost_right, node->num_procs_left/2);
}

int KDTree::getLongestDir(const Box& box) {
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

void KDTree::splitBox(int split, int dir,
                      const Box& box, Box& left, Box& right) {
    left = box;
    right = box;
    
    left.setBig(dir, split);
    right.setSmall(dir, split+1);        
}
