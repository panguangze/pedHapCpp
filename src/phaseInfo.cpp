//
// Created by caronkey on 10/23/21.
//

#include "phaseInfo.h"

PInfo::PInfo(int infoId) {
    this->infoId = -infoId;
    this->firstBlock = -1;
}

void PInfo::set_covered_call(int b_ps, int side, int pos) {
    if(firstBlock == -1) firstBlock = b_ps;
    if (blocks.find(b_ps) != blocks.end()) {
        if (side == 0) {
            side0_support[b_ps]->push_back(pos);
        }
        if (side == 1) {
            side1_support[b_ps]->push_back(pos);
        }
    } else {
        blocks.emplace(b_ps, 1);
        side0_support[b_ps] = new std::vector<int>();
        side1_support[b_ps] = new std::vector<int>();
    }
}

void PInfo::init_blocks() {
    for(auto it: blocks) {
        int block_id = it.first;
        int side = it.second;
        if (P_ENSURE_BLOCK !=0 && block_id == P_ENSURE_BLOCK && P_ENSURE_SIDE == 0) {
            block_reverses[block_id] = false;
            continue;
        }
        if (block_id == infoId) continue;
        auto v0 = (float)side0_support[block_id]->size();
        auto v1 = (float)side1_support[block_id]->size();
        if (v1 == v0 || (v0 != 0 && v1 != 0 && std::max(v0, v1)/std::min(v0, v1) <= T1) || (abs(v1 - v0) <= T2 && (v1 ==0 || v0 == 0))){
            uncertain_blocks.push_back(block_id);
        } else {
            auto n_r = false;
            if (v1 > v0) {
                n_r = true;
                it.second = 0;
            }
            block_reverses[block_id] = n_r;
        }
    }
    blocks.emplace(infoId, -1);
    block_reverses[infoId] = false;
}
