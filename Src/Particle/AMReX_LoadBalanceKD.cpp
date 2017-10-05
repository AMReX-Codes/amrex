#include "AMReX_LoadBalanceKD.H"

using namespace amrex;

int KDTree::min_box_size = 4;

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
    
    bool success = partitionNode(node, cost);
    
    if (success) {
        buildKDTree(node->left,  cost);
        buildKDTree(node->right, cost);
    }
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

bool KDTree::partitionNode(KDNode* node, const FArrayBox& cost) {
    
    const Box& box = node->box;
    BL_ASSERT(cost.box().contains(box));
    
    int split;
    Real cost_left, cost_right;
    int dir = getLongestDir(box);
    amrex_compute_best_partition(cost.dataPtr(), cost.loVect(), cost.hiVect(),
                                 box.loVect(), box.hiVect(), node->cost, dir,
                                 &cost_left, &cost_right, &split);
    
    Box left, right;
    bool success = splitBox(split, dir, box, left, right);

    BL_ASSERT(left.numPts()  > 0);
    BL_ASSERT(right.numPts() > 0);
    
    if (success) {
        node->left  = new KDNode(left,  cost_left,  node->num_procs_left/2);
        node->right = new KDNode(right, cost_right, node->num_procs_left/2);
    }

    return success;
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

bool KDTree::splitBox(int split, int dir,
                      const Box& box, Box& left, Box& right) {

    int lo = box.smallEnd(dir);
    int hi = box.bigEnd(dir);
    int length = hi - lo + 1;

    if (length < 2*min_box_size) return false;

    int new_length_left  = split - lo + 1;
    int new_length_right = hi - split; 
    if (new_length_left < min_box_size) {
        split = min_box_size + lo - 1;
    } else if (new_length_right < min_box_size) {
        split = hi - min_box_size;
    }

    left = box;
    right = box;
    
    left.setBig(dir, split);
    right.setSmall(dir, split+1);

    return true;
}
