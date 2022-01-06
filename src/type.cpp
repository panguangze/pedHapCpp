//
// Created by yyh on 12/10/2018.
//
#include "type.h"
#include <iostream>
#include <cmath>
#include <set>
#include "util.h"
#define SPECIFIC_HOMO_BLOCK -123456789

ChromoPhaser::ChromoPhaser(const uint &chr_id, const std::string &chr_name, int nsmp) {
    this->chr_id = chr_id;
    this->chr_name = chr_name;
    this->sample_count = nsmp;
    this->variant_count = 0;
}


void ChromoPhaser::add_result(const std::shared_ptr<VcfRecord>& result) {
    this->results_for_variant.push_back(result);
}

void ChromoPhaser::phase_with_hete(int idx1, int idx2, int side) {
    logging(std::cerr,"phasing hete");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    std::unordered_map<uint, PInfo*> reads;
    for(int mendel_pas : this->mendel_pass) {
        auto result = results_for_variant[mendel_pas];
        if(result->bnd) continue;
        Call* s1_call = result->calls[idx1];
        Call* s2_call = result->calls[idx2];
        if( s1_call->isHomo() || s2_call->isHomo() || !s2_call->isPhased()) continue;
//        check mendel?
//        if(s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2 &&
//            s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2)
//            continue;
        if(reads.find(s2_call->block_id) == reads.end()) {
            auto read = new PInfo(-s2_call->block_id);
            reads[s2_call->block_id] = read;
        }

        if(s1_call->isPhased()) {
            auto o_side = side;
            if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                o_side = abs(side -1);
            } else {
                o_side = abs(side);
            }
            reads[s2_call->block_id]->set_covered_call(s1_call->block_id, o_side, s1_call->pos);
        } else {
            if ((s1_call->allele1 != s2_call->allele2 && s1_call->allele1 != s2_call->allele1) ||
                    (s1_call->allele2 != s2_call->allele2 && s1_call->allele2 != s2_call->allele1))
                continue;
            if (s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            }
            s1_call->block_id = -s2_call->block_id;
        }
    }
    InfoSet hete_reads;
    auto s = reads.size();
    for(auto it: reads) {
        hete_reads.add_read(it.second);
    }
    extend(idx1, hete_reads, side);
}

void ChromoPhaser::extend(int idx, InfoSet& infoSet, int side) {
    std::unordered_map<int ,int> finalize_new_block_ids;
    int f_new_id;
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        Call* s1_call = result->calls[idx];
        if (s1_call->isHomo()) continue;
        if (s1_call->block_id == 0) continue;
        auto new_block_id = infoSet.find(s1_call->block_id);
        auto n_r = infoSet.blocks_reverse_info[s1_call->block_id];
        if (new_block_id == -1) {
            if (finalize_new_block_ids.find(s1_call->block_id) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[s1_call->block_id];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                finalize_new_block_ids.emplace(s1_call->block_id, f_new_id);
            }
        } else {
            if(finalize_new_block_ids.find(new_block_id) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[new_block_id];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                finalize_new_block_ids.emplace(new_block_id, f_new_id);
            }
        }
        s1_call->block_id = f_new_id;
        if(n_r){
            auto t = s1_call->allele1;
            s1_call->allele1 = s1_call->allele2;
            s1_call->allele2 = t;
        }
    }
    for(auto it : finalize_new_block_ids) {
        if (it.first != it.second) {
            P_ENSURE_BLOCK = it.first;
            P_ENSURE_SIDE = side;
            break;
        }
    }
}

void ChromoPhaser::phase_with_homo(int idx1, int idx2, int side) {
    logging(std::cerr,"phasing homo");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    for(int mendel_pas : this->mendel_pass) {
        auto result = results_for_variant[mendel_pas];
        Call *s1_call = result->calls[idx1];
        Call *s2_call = result->calls[idx2];
        if (s1_call->isHomo() || !s2_call->isHomo()) continue;
//        check mendel?
        if ((s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2) &&
                (s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2))
            continue;
//        if (reads.find(s2_call->block_id) == reads.end()) {
//            auto read = new PInfo(s2_call->block_id);
//            reads[s2_call->block_id] = read;
//        }

        if (s1_call->isPhased()) {
            auto o_side = side;
            if (s1_call->allele2 == s2_call->allele1) {
                o_side = abs(side - 1);
            } else {
                o_side = abs(side);
            }
            read->set_covered_call(s1_call->block_id, o_side, s1_call->pos);
        } else {
            s1_call->block_id = SPECIFIC_HOMO_BLOCK;
            if(side == 0 && s1_call->allele1 != s2_call->allele1) {
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            } else if(side == 1 && s1_call->allele1 == s2_call->allele1){
                auto t = s1_call->allele1;
                s1_call->allele1 = s1_call->allele2;
                s1_call->allele2 = t;
            }
        }
    }
    InfoSet hete_reads;
    hete_reads.add_read(read);
    extend(idx1, hete_reads,side);
}

bool ChromoPhaser::check_mendel(int idx1, int idx2, int idx3) {
    for (int i = 0; i < results_for_variant.size(); ++i) {
        auto rec = results_for_variant[i];
        auto child = rec->calls[idx1];
        auto father = rec->calls[idx2];
        auto mother = rec->calls[idx3];
        auto c1 = child->allele1;
        auto c2 = child->allele2;
        auto f1 = father->allele1;
        auto f2 = father->allele2;
        auto m1 = mother->allele1;
        auto m2 = mother->allele2;
        if (((c1 == f1 || c1 == f2) && (c2 == m1 || c2 == m2)) || ((c2 == f1 || c2 == f2) && (c1 == m1 || c1 == m2))) mendel_pass.push_back(i);
        mendel_cs.push_back(i);
    }
}
