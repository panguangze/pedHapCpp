//
// Created by caronkey on 10/23/21.
//

#include "phaseInfo.h"

PInfo::PInfo(int infoId) {
    this->infoId = infoId;
}

void PInfo::set_covered_call(int b_ps, int side, int pos, bool is_homo_support) {
    if (blocks.find(b_ps)== blocks.end()) {
        blocks.emplace(b_ps, 1);
        side0_support[b_ps] = new std::vector<int>();
        side1_support[b_ps] = new std::vector<int>();
    }
    if (side == 0) {
        side0_support[b_ps]->push_back(pos);
        if (is_homo_support) {
            side1_support[b_ps]->push_back(pos);
        }
    }
    if (side == 1) {
        side1_support[b_ps]->push_back(pos);
        if (is_homo_support) {
            side0_support[b_ps]->push_back(pos);
        }
    }
}

void PInfo::init_blocks(std::vector<int>& confilict_poses) {
    for(auto it: blocks) {
        int block_id = it.first;
        int side = it.second;
//        if (P_ENSURE_BLOCK !=0 && block_id == P_ENSURE_BLOCK && P_ENSURE_SIDE == 0) {
//            block_reverses[block_id] = false;
//            continue;
//        }
        if (block_id == infoId) continue;
        auto v0 = (float)side0_support[block_id]->size();
        auto v1 = (float)side1_support[block_id]->size();
//        if (v1 == v0 || (v0 != 0 && v1 != 0 && std::max(v0, v1)/std::min(v0, v1) <= T1) || (abs(v1 - v0) <= T2 && (v1 ==0 || v0 == 0))){
        if (std::abs(v1 - v0) >= T2 ||((v1 == 0 || v0 == 0) && std::abs(v1 - v0) >= T1)){
            certain_blocks.push_back(block_id);
            auto n_r = false;
            if (v1 > v0) {
                n_r = true;
                it.second = 0;
                for (auto item : *side0_support[block_id]) {
                    confilict_poses.push_back(item);
                }
            } else {
                for (auto item : *side1_support[block_id]) {
                    confilict_poses.push_back(item);
                }
            }
            block_reverses[block_id] = n_r;
//if (v1 == v0 || (v0 != 0 && v1 != 0 && std::min(v1, v0) / std::max(v1, v0)  > T1) || (std::abs(v1 - v0) < T2 && (v1 ==0 || v0 == 0))){
        } else {
            uncertain_blocks.push_back(block_id);
        }
    }
    certain_blocks.push_back(infoId);
    block_reverses[infoId] = false;
}
