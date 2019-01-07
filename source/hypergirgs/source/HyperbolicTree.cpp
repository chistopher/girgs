
#include <hypergirgs/HyperbolicTree.h>


namespace hypergirgs {


HyperbolicTree::HyperbolicTree(std::vector<double> &radii, std::vector<double> &angles, double T, double R) {

}

std::vector<std::pair<int, int>> HyperbolicTree::generate(int seed) {
    return std::vector<std::pair<int, int>>();
}

void HyperbolicTree::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {

    auto touching = m_helper.touching(cellA, cellB, level);
    if(cellA == cellB || touching) {
        // sample all type 1 occurrences with this cell pair
        for(auto& layer_pair : m_layer_pairs[level]){
            if(cellA != cellB || layer_pair.first <= layer_pair.second)
                sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
        }

    } else { // not touching
        if (m_T == 0)
            return;
        // sample all type 2 occurrences with this cell pair
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    if(touching) {
        // recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they don't
        auto fA = AngleHelper::firstChild(cellA);
        auto fB = AngleHelper::firstChild(cellA);
        visitCellPair(fA + 0, fB + 0, level+1);
        visitCellPair(fA + 0, fB + 1, level+1);
        visitCellPair(fA + 1, fB + 1, level+1);
        if(cellA != cellB)
            visitCellPair(fA + 1, fB + 0, level+1); // if they are the same we already did this call 3 lines above
    }
}

void HyperbolicTree::sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j) {

}

void HyperbolicTree::sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j) {

}


} // namespace hypergirgs
