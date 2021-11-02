//
// Created by yyh on 12/10/2018.
//
#include "type.h"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#define SPECIFIC_HOMO_BLOCK -123456789

Fragment::Fragment(const Fragment &rhs)
{
    this->read_qual = rhs.read_qual;
    this->snps = rhs.snps;
    this->insertion_size = rhs.insertion_size;
    this->rescued = rhs.rescued;
    this->dm = rhs.dm;
    this->start = rhs.start;
    this->end = rhs.end;
}

//------------------------------- PhasedBlock Class ------------------------------------//

PhasedBlock::PhasedBlock() { std::cout << "You shouldn't be here"; }
PhasedBlock::PhasedBlock(uint start_idx)
{
    this->start_variant_idx = start_idx;
    this->end_variant_idx = start_idx;
}
//PhasedBlock::PhasedBlock(uint variant_start_idx, std::shared_ptr<VcfRecord> result)
//{
//    this->start_variant_idx = variant_start_idx;
//    calls[variant_start_idx] = result;
//    variant_idxes.insert(start_variant_idx);
//}

//do not use this constructor
PhasedBlock::PhasedBlock(const PhasedBlock &rhs)
{
    start_variant_idx = rhs.start_variant_idx;
    calls = rhs.calls;                  //bug here, should make duplicate
    variant_idxes = rhs.variant_idxes;
}

//PhasedBlock::PhasedBlock(const PhasedBlock &rhs, std::unordered_map<uint, std::shared_ptr<VcfRecord>> &results_dup)
//{
//    start_variant_idx = rhs.start_variant_idx;
//    variant_idxes = rhs.variant_idxes;
//    for (uint var_idx : variant_idxes)
//        calls[var_idx] = results_dup[var_idx];
//}

PhasedBlock::~PhasedBlock()
{
    calls.clear();
    variant_idxes.clear();
}

//PhasedBlock::PhasedBlock(const PhasedBlock &rhs) {}

void PhasedBlock::flip(int idx)
{
    for (auto it : calls)
        it.second->flip();
}


ChromoPhaser::ChromoPhaser(const uint &chr_id, const std::string &chr_name, int nsmp) {
    this->chr_id = chr_id;
    this->chr_name = chr_name;
    this->sample_count = nsmp;
    for (int i = 0; i < sample_count; ++i) {
        auto tmp = new sample_PhasedBlocks();
        this->phased_blocks_info.push_back(tmp);
    }
}

ChromoPhaser::~ChromoPhaser() {

}

void ChromoPhaser::setSampleCount(int count) {
    this->sample_count = count;
    for(int i = 0; i< count; i++) {
        auto tmp = new std::unordered_map<uint, uint>();
        this->variant_to_block_id.push_back(tmp);
    }
}

void ChromoPhaser::add_result(std::shared_ptr<VcfRecord> result) {
    this->results_for_variant.push_back(result);
//    for(int i =0; i < this->sample_count; i++) {
//        auto s_blocks_map = this->phased_blocks_info[i];
//        auto* tcall = result->calls[i];
//        if (s_blocks_map->find(tcall->ps) != s_blocks_map->end()) {
//            (*s_blocks_map)[tcall->ps]->add_call(tcall);
////            (*s_blocks_map)[tcall->ps]->calls.emplace(tcall->pos, tcall);
//        } else {
//            (*s_blocks_map)[tcall->ps] = std::make_shared<PhasedBlock>(tcall->pos);
//        }
//    }
}

void ChromoPhaser::phase_with_hete(int idx1, int idx2, int side) {
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    std::unordered_map<uint, PInfo*> reads;
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        Call* s1_call = result->calls[idx1];
        Call* s2_call = result->calls[idx2];
        if( s1_call->isHomo() || s2_call->isHomo() || !s2_call->isPhased()) continue;
//        check mendel?
        if(s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2 &&
            s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2)
            continue;
        if(reads.find(s2_call->ps) == reads.end()) {
            auto read = new PInfo(s2_call->ps);
            reads[s2_call->ps] = read;
        }

        if(s1_call->isPhased()) {
            auto o_side = side;
            if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                o_side = abs(side -1);
            } else {
                o_side = abs(side);
            }
            reads[s2_call->ps]->set_covered_call(s1_call->ps, o_side, s1_call->pos);
        } else {
            if (s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            }
            s1_call->ps = -s2_call->ps;
        }
    }
    InfoSet hete_reads;
    for(auto it: reads) {
        hete_reads.add_read(it.second);
    }
    extend(idx1, hete_reads);
}

void ChromoPhaser::extend(int idx, InfoSet& infoSet) {
    std::unordered_map<int ,int> finalize_new_block_ids;
    int f_new_id;
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        Call* s1_call = result->calls[idx];
        if (s1_call->isHomo()) continue;
        if (s1_call->ps == 0) continue;
        auto new_block_id = infoSet.find(s1_call->ps);
        auto n_r = infoSet.blocks_reverse_info[s1_call->ps];
        if (new_block_id == -1) {
            if (finalize_new_block_ids.find(s1_call->ps) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[s1_call->ps];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                finalize_new_block_ids.emplace(s1_call->ps, f_new_id);
            }
        } else {
            if(finalize_new_block_ids.find(new_block_id) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[new_block_id];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                finalize_new_block_ids.emplace(new_block_id, f_new_id);
            }
        }
        s1_call->ps = f_new_id;
        if(n_r){
            auto t = s1_call->allele1;
            s1_call->allele1 = s1_call->allele2;
            s1_call->allele2 = t;
        }
    }
    for(auto it : finalize_new_block_ids) {
        if (it.first != it.second) {
            P_ENSURE_BLOCK = it.first;
            P_ENSURE_SIDE = it.second;
        }
    }
}

void ChromoPhaser::phase_with_homo(int idx1, int idx2, int side) {
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        Call *s1_call = result->calls[idx1];
        Call *s2_call = result->calls[idx2];
        if (s1_call->isHomo() || !s2_call->isHomo()) continue;
//        check mendel?
        if ((s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2) ||
                (s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2))
            continue;
//        if (reads.find(s2_call->ps) == reads.end()) {
//            auto read = new PInfo(s2_call->ps);
//            reads[s2_call->ps] = read;
//        }

        if (s1_call->isPhased()) {
            auto o_side = side;
            if (s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                o_side = abs(side - 1);
            } else {
                o_side = abs(side);
            }
            read->set_covered_call(s1_call->ps, o_side, s1_call->pos);

        } else {
            s1_call->ps = SPECIFIC_HOMO_BLOCK;
            if(side == 0 && s1_call->allele1 != s2_call->allele1) {
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            } else if(side == 0 && s1_call->allele1 == s2_call->allele2){
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            }
        }
    }
    InfoSet hete_reads;
    hete_reads.add_read(read);
    extend(idx1, hete_reads);
}
