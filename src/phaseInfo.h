//
// Created by caronkey on 10/23/21.
//

#ifndef PEDHAP_PHASEINFO_H
#define PEDHAP_PHASEINFO_H
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "results.h"
extern float T1;
extern float T2;
extern int P_ENSURE_BLOCK;
extern int P_ENSURE_SIDE;
class PInfo{
public:
//    std::unordered_map<int,int> covered_blocks;
    std::unordered_map<int, int> blocks; //b_id : confilict
    std::unordered_map<int,bool> block_reverses;
    std::vector<int> uncertain_blocks;
    std::unordered_map<int, std::vector<int>*> side0_support;
    std::unordered_map<int, std::vector<int>*> side1_support;
    int infoId;
    int firstBlock;
    explicit PInfo(int infoId);
    void set_covered_call(int ps, int side, int pos);
    void init_blocks();
};

class InfoSet {
public:
    std::vector<int> rank;
    std::unordered_map<int, int> parent;
    std::unordered_map<int, bool> blocks_reverse_info;
    std::vector<int> uncertain_blocks;
    std::vector<int> confilict_poses;
public:
    // Constructor to create and
    // initialize sets of n items
    InfoSet() = default;

    bool n_reverse_info(int block_id) {
        if(blocks_reverse_info.find(block_id) == blocks_reverse_info.end()) return false;
        return blocks_reverse_info[block_id];
    }

    void add_read(PInfo* pinfo){
        pinfo->init_blocks();
//        ensure block is better
        auto first_block = pinfo->blocks[pinfo->firstBlock];
        auto first_reverse = pinfo->block_reverses[pinfo->firstBlock];
//        if first block found and conflict exists.
        if(this->blocks_reverse_info.find(first_block) != blocks_reverse_info.end()
                    && first_reverse != this->blocks_reverse_info[first_block]) {
            for(int i = 0; i < pinfo->blocks.size(); i++) {
                pinfo->block_reverses[i] = !pinfo->block_reverses[i];
            }
        }

        for(auto blk_id : pinfo->blocks) {
            auto b_id = blk_id.first;
            auto r = pinfo->block_reverses[b_id];
            if (parent.find(b_id) != parent.end()) {
                this->parent[b_id] = b_id;
            }
            this->blocks_reverse_info.emplace(b_id, r);
            Union(b_id, first_block);
        }
    }

    // Finds set of given item x
    int find(int x)
    {
        // Finds the representative of the set
        // that x is an element of
        if(parent.find(x) == parent.end()) return -1;
        if (parent[x] != x) {

            // if x is not the parent of itself
            // Then x is not the representative of
            // his set,
            parent[x] = find(parent[x]);

            // so we recursively call Find on its parent
            // and move i's node directly under the
            // representative of this set
        }

        return parent[x];
    }

    // Do union of two sets represented
    // by x and y.
    void Union(int x, int y)
    {
        // Find current sets of x and y
        int xset = find(x);
        int yset = find(y);

        // If they are already in same set
        if (xset == yset)
            return;

        // Put smaller ranked item under
        // bigger ranked item if ranks are
        // different
        if (rank[xset] < rank[yset]) {
            parent[xset] = yset;
        }
        else if (rank[xset] > rank[yset]) {
            parent[yset] = xset;
        }

            // If ranks are same, then increment
            // rank.
        else {
            parent[yset] = xset;
            rank[xset] = rank[xset] + 1;
        }
    }
};


#endif //PEDHAP_PHASEINFO_H