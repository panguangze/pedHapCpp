//
// Created by caronkey on 10/23/21.
//

#ifndef PEDHAP_PHASEINFO_H
#define PEDHAP_PHASEINFO_H
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "results.h"
#include <map>
extern float T1;
extern float T2;
extern int P_ENSURE_BLOCK;
extern int P_ENSURE_SIDE;
extern int CORRECT_SCORE;
extern int ERROR_SCORE;
class PInfo{
public:
//    std::unordered_map<int,int> covered_blocks;
    std::map<int, int> blocks; //b_id : confilict
    std::unordered_map<int,bool> block_reverses;
    std::vector<int> uncertain_blocks;
    std::vector<int> certain_blocks;
    std::map<int, std::vector<int>*> side0_support;
    std::map<int, std::vector<int>*> side1_support;
    int infoId;
    explicit PInfo(int infoId);
    void set_covered_call(int ps, int side, int pos, bool is_homo_support);
    void init_blocks(std::vector<int>& confilict_poses);
};

class InfoSet {
public:
    std::unordered_map<int, int> rank;
    std::unordered_map<int, int> parent;
    std::unordered_map<int, bool> blocks_reverse_info;
    std::vector<int> uncertain_blocks;
    std::vector<int> confilict_poses;
    std::unordered_map<int, std::vector<int>*> side0_support;
    std::unordered_map<int, std::vector<int>*> side1_support;
    std::vector<int> prev_link_heads;
//    std::unordered_map<int, std::vector<int>*> root_leaves;
public:
    // Constructor to create and
    // initialize sets of n items
    InfoSet() = default;

    bool n_reverse_info(int block_id) {
        if(blocks_reverse_info.find(block_id) == blocks_reverse_info.end()) return false;
        return blocks_reverse_info[block_id];
    }

    void flip(int x) {
        auto root = parent[x];
        for (auto item : parent) {
            if(item.second == root) {
                blocks_reverse_info[item.first] = !blocks_reverse_info[item.first];
            }
        }
    }

    void interact(InfoSet* other) {

    }

    void add_read(PInfo* pinfo, bool is_homo){
        pinfo->init_blocks(confilict_poses);
//        ensure block is better
        auto first_block = pinfo->certain_blocks[0];
        prev_link_heads.push_back(first_block) ;
        auto first_reverse = pinfo->block_reverses[first_block];
//        if first block found and conflict exists.
        if (!is_homo) {
            if(this->blocks_reverse_info.find(first_block) != blocks_reverse_info.end()
               && first_reverse != this->blocks_reverse_info[first_block]) {
                for(auto it: pinfo->blocks) {
                    pinfo->block_reverses[it.first] = !pinfo->block_reverses[it.first];
                }
            }
        }
        auto prev_block_id = first_block;
        for(auto blk_id : pinfo->certain_blocks) {
            if (side0_support.find(blk_id) == side0_support.end()) {
//                blocks.emplace(b_ps, 1);
                side0_support[blk_id] = new std::vector<int>();
                side1_support[blk_id] = new std::vector<int>();
                side0_support[blk_id] = pinfo->side0_support[blk_id];
                side1_support[blk_id] = pinfo->side1_support[blk_id];
            }
            auto r = pinfo->block_reverses[blk_id];
//            if (parent.find(b_id) == parent.end()) {
//                this->parent[b_id] = b_id;
//                this->rank[b_id] = 1;
//            }
            if (blk_id > 0) {
                if(find(blk_id) != -1) {
                    // if current max and two flip not equal, flip all origin
                    if (blocks_reverse_info.find(blk_id) != blocks_reverse_info.end() && blocks_reverse_info[blk_id] != r) {
                        auto origin_side0_size = side0_support[blk_id]->size();
                        auto origin_side1_size = side1_support[blk_id]->size();
                        auto current_side0_size = pinfo->side0_support[blk_id]->size();
                        auto current_side1_size = pinfo->side1_support[blk_id]->size();

                        auto origin_max = origin_side0_size > origin_side1_size ? origin_side0_size : origin_side1_size;
                        auto origin_min = origin_side0_size < origin_side1_size ? origin_side0_size : origin_side1_size;
                        auto current_max = current_side0_size > current_side1_size ? current_side0_size : current_side1_size;
                        auto current_min = current_side0_size < current_side1_size ? current_side0_size : current_side1_size;
                        if (std::find(prev_link_heads.begin(), prev_link_heads.end(), blk_id) != prev_link_heads.end()) {
//                            blocks_reverse_info.emplace(blk_id, r);
                            flip(blk_id);
                        } else {
                            if (CORRECT_SCORE * origin_max - ERROR_SCORE * origin_min <= CORRECT_SCORE * current_max - ERROR_SCORE * current_min) {
//                            flip(blk_id);
                                blocks_reverse_info.emplace(blk_id, r);
                            }
                        }
                        Union(blk_id, prev_block_id);
                        prev_block_id = blk_id;
                    } else {
                        blocks_reverse_info.emplace(blk_id, r);
                        Union(blk_id, prev_block_id);
                        prev_block_id = blk_id;
                    }
                } else {
                    this->blocks_reverse_info.emplace(blk_id, r);
                    Union(blk_id, prev_block_id);
                    prev_block_id = blk_id;
                }
            }
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
        if(xset == -1) {
            this->parent[x] = x;
            this->rank[x] = 1;
            xset = x;
        }
        if(yset == -1) {
            this->parent[y] = y;
            this->rank[y] = 1;
            yset = y;
        }

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
