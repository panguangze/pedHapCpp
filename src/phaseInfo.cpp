//
// Created by caronkey on 10/23/21.
//

#include "phaseInfo.h"

PInfo::PInfo(int ps) {
    this->ps = ps;
}

void PInfo::set_covered_call(int b_ps, int side, int pos) {
    if (covered_blocks.find(b_ps) != covered_blocks.end()) {
        if (side == 0) {
            if(side0_support.find(b_ps) == side0_support.end()) {
                side0_support[b_ps] = new std::vector<int>();
            }
            side0_support[b_ps]->push_back(pos);
        }
        if (side == 1) {
            if(side1_support.find(b_ps) == side1_support.end()) {
                side1_support[b_ps] = new std::vector<int>();
            }
            side1_support[b_ps]->push_back(pos);
        }
    }
}

void PInfo::init_blocks(std::vector<int> &ensure_blocks) {

}
