//
// Created by yyh on 12/10/2018.
//
#include "type.h"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>

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
        if( s1_call->isHomo() || s2_call->isHomo() || !s2_call->phased) continue;
//        check mendel?
        if(s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2 &&
            s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2)
            continue;
        if(reads.find(s2_call->ps) == reads.end()) {
            auto read = new PInfo(s2_call->ps);
            reads[s2_call->ps] = read;
        }

        if(s1_call->phased) {
            auto o_side = side;
            if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                o_side = abs(side -1);
            } else {
                o_side = abs(side);
            }
            reads[s2_call->ps]->set_covered_call(s1_call, o_side, s1_call->pos);

        } else {
            s1_call->ps = -s2_call->ps;
            s1_call->allele1 = s2_call->allele1;
            s2_call->allele2 = s2_call->allele2;
        }
        auto hete_reads = new InfoSet();
        for(auto it: reads) {
            hete_reads.add_read(it.second)
        }
    }

}
