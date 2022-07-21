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

void ChromoPhaser::correct_conflict(int idx){
    std::set<int> intersect;

    set_intersection(conflicts1.begin(), conflicts1.end(), conflicts2.begin(), conflicts2.end(),
                     std::inserter(intersect, intersect.begin()));
    for (auto item : intersect) {
        if (item == 168385) {
            int ddd=89;
        }
        if (results_for_variant[item]->pos == 168385) {
            int tmp = 99;
        }
        results_for_variant[item]->calls[idx]->flip();
        auto tmp = results_for_variant[item]->pos;
        auto tmp2 = 0;
    }
    conflicts1.clear();
    conflicts2.clear();
}

void ChromoPhaser::phase_with_hete(int idx1, int idx2, int side, InfoSet* infoSet) {
    logging(std::cerr,"phasing hete");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    std::unordered_map<uint, PInfo*> reads;
    auto tmp = (idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0) ? this->mendel_passf : this->mendel_passm;
    int current_link_block_id = 2;
    int prev_s2_block_id = 0;
    std::unordered_map<int, int> s2_block_id2current_id;
    for(int mendel_pas : tmp) {
        auto result = results_for_variant[mendel_pas];
        if (result->pos == 168404) {
            int tmp3=1;
        }
//        if(result->bnd) continue;
        Call* s1_call = result->calls[idx1];
        Call* s2_call = result->calls[idx2];
        if( s1_call->isHomo() || (!s2_call->isHomo() and !s2_call->isPhased()) || s2_call->isHomo()) continue;
        if ((!s2_call->isHomo() && s2_call->block_id != prev_s2_block_id)) {
            if (s2_block_id2current_id.find(s2_call->block_id) != s2_block_id2current_id.end()) {
                current_link_block_id = s2_block_id2current_id[s2_call->block_id];
                prev_s2_block_id = s2_call->block_id;
            } else {
                if (prev_s2_block_id != 0) {
                    current_link_block_id = current_link_block_id + 1;
                }
                prev_s2_block_id = s2_call->block_id;
                s2_block_id2current_id.emplace(s2_call->block_id,current_link_block_id);
            }
//            if (prev_s2_block_id == 0) continue;
        }
//        if( s1_call->isHomo() || s2_call->isHomo() || !s2_call->isPhased()) continue;

//        check mendel?
//        if(s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2 &&
//            s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2)
//            continue;
//        if (s2_call->isHomo())
        if(reads.find(current_link_block_id) == reads.end()) {
            auto read = new PInfo(-current_link_block_id);
            reads[current_link_block_id] = read;
        }

        if(s1_call->isPhased()) {
            auto o_side = side;
            if (s2_call->isHomo()) {
                if(s1_call->allele1 == s2_call->allele1) {
                    o_side = abs(side);
                } else {
                    o_side = abs(side - 1);
                }
                reads[current_link_block_id]->set_covered_call(s1_call->block_id, o_side, mendel_pas);
            } else {
                if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
                    o_side = abs(side -1);
                } else {
                    o_side = abs(side);
                }
                reads[current_link_block_id]->set_covered_call(s1_call->block_id, o_side, mendel_pas);
            }
//            auto o_side = side;
        } else {
//            if ((s1_call->allele1 != s2_call->allele2 && s1_call->allele1 != s2_call->allele1) &&
//                    (s1_call->allele2 != s2_call->allele2 && s1_call->allele2 != s2_call->allele1))
//                continue;
//            if (s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
//                if (side % 2 == 0) {
//                    conflicts1.insert(mendel_pas);
//                } else {
//                    conflicts2.insert(mendel_pas);
//                }
//                auto t = s1_call->allele1;
//                s1_call->allele1 = s1_call->allele2;
//                s1_call->allele2 = t;
//            }
//            s1_call->block_id = -s2_call->block_id;
        }
    }
//    InfoSet hete_reads;
    auto s = reads.size();
    for(auto it: reads) {
        if (it.first == 2) {
            auto mmm = 33;
        }
        infoSet->add_read(it.second, false);
    }
//    extend(idx1, hete_reads, side);
    if (idx1 == 0) {
        for (auto item : infoSet->confilict_poses) {
            if (side % 2 == 0) {
                conflicts1.insert(item);
            } else {
                conflicts2.insert(item);
            }
        }
    }
}

void ChromoPhaser::extend(int idx, InfoSet* infoSet, int side) {
    std::unordered_map<int ,int> finalize_new_block_ids;
    int f_new_id;
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        if (result->pos == 168386) {
            int tmp=1;
        }
        Call* s1_call = result->calls[idx];
        if (s1_call->isHomo()) continue;
        if (s1_call->block_id == 0) continue;
        if (result->pos == 168386) {
            int tmp = 2;
            auto v = infoSet->blocks_reverse_info.find(69538) == infoSet->blocks_reverse_info.end();
            int tp = 2;
        }
        auto new_block_id = infoSet->find(s1_call->block_id);
        auto n_r = infoSet->blocks_reverse_info[s1_call->block_id];
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

void ChromoPhaser::phase_with_homo(int idx1, int idx2, int side, InfoSet* infoSet) {
    logging(std::cerr,"phasing homo");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    auto tmp = ((idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0)) ? this->mendel_passf : this->mendel_passm;

    for(int i = 0 ; i < tmp.size(); i++) {
        auto mendel_pas = tmp[i];
        auto result = results_for_variant[mendel_pas];
        if (result->pos == 97921) {
            int tmp3=1;
        }
        Call *s1_call = result->calls[idx1];
        Call *s2_call = result->calls[idx2];
        if(s1_call->block_id == 11) {
            int mmm = 33;
        }
        if (s1_call->isHomo() || (!s2_call->isHomo() )) continue;
//        check mendel?
//        if ((s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2) &&
//                (s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2))
//            continue;
//        if (reads.find(s2_call->block_id) == reads.end()) {
//            auto read = new PInfo(s2_call->block_id);
//            reads[s2_call->block_id] = read;
//        }

        if (s1_call->isPhased()) {
            auto o_side = side;
            if (o_side == 0 && this->conflicts1.find(mendel_pas) != this->conflicts1.end())  {
                if (s1_call->allele2 == s2_call->allele1) {
                    o_side = abs(side);
                } else {
                    o_side = abs(side - 1);
                }
            } else {
                if (o_side == 0) {
                    if (s1_call->allele1 == s2_call->allele1) {
                        o_side = abs(side);
                    } else {
                        o_side = abs(side - 1);
                    }
                } else {
                    if (s1_call->allele2 == s2_call->allele1) {
                        o_side = abs(side - 1);
                    } else {
                        o_side = abs(side);
                    }
                }
            }
            read->set_covered_call(s1_call->block_id, o_side, mendel_pas);
        } else {
//            s1_call->block_id = SPECIFIC_HOMO_BLOCK;
//            if(side == 0 && s1_call->allele1 != s2_call->allele1) {
//                auto t = s1_call->allele1;
//                s1_call->allele1 = s1_call->allele2;
//                s1_call->allele2 = t;
//            } else if(side == 1 && s1_call->allele1 == s2_call->allele1){
//                auto t = s1_call->allele1;
//                s1_call->allele1 = s1_call->allele2;
//                s1_call->allele2 = t;
//            }
        }
    }
//    InfoSet hete_reads;
    infoSet->add_read(read, true);
//    extend(idx1, hete_reads,side);
    if (idx1 == 0) {
        for (auto item : infoSet->confilict_poses) {
            if (side % 2 == 0) {
                conflicts1.insert(item);
            } else {
                conflicts2.insert(item);
            }
        }
    }
}


void ChromoPhaser::phase_with_homo2(int idx1, int idx2, int side, InfoSet* infoSet) {
    logging(std::cerr,"phasing homo");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    auto tmp = ((idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0)) ? this->mendel_passf : this->mendel_passm;

    for(int i = 0 ; i < tmp.size(); i++) {
        auto mendel_pas = tmp[i];
        auto result = results_for_variant[mendel_pas];
        if (result->pos == 1236104) {
            int tmp3=1;
        }
        Call *s1_call = result->calls[idx1];
        Call *s2_call = result->calls[idx2];
        if(s1_call->block_id == 11) {
            int mmm = 33;
        }
//        if (!) {
//
//        }
//        if (s1_call->isHomo() || (!s2_call->isHomo() )) continue;
        if ((!s2_call->isPhased() && !s2_call->isHomo()) || (s2_call->isPhased() && s2_call->block_id != 1)) continue;
//        check mendel?
//        if ((s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2) &&
//                (s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2))
//            continue;
//        if (reads.find(s2_call->block_id) == reads.end()) {
//            auto read = new PInfo(s2_call->block_id);
//            reads[s2_call->block_id] = read;
//        }

        if (s1_call->isPhased()) {
            auto o_side = side;
            if (o_side == 0) {
                if (s2_call->allele1 == s1_call->allele1) {
                    o_side = abs(side);
                } else {
                    o_side = abs(side - 1);
                }
            } else {
                if (s2_call->allele2 == s1_call->allele1) {
                    o_side = abs(side - 1);
                } else {
                    o_side = abs(side);
                }
            }
//            if (o_side == 0 && this->conflicts1.find(mendel_pas) != this->conflicts1.end())  {
//                if (s1_call->allele2 == s2_call->allele1) {
//                    o_side = abs(side);
//                } else {
//                    o_side = abs(side - 1);
//                }
//            } else {
////                if (o_side == 0) {
//                    if (s1_call->allele1 == s2_call->allele1) {
//                        o_side = 0;
//                    } else {
//                        o_side = 1;
//                    }
////                } else {
////                    if (s1_call->allele2 == s2_call->allele1) {
////                        o_side = abs(side - 1);
////                    } else {
////                        o_side = abs(side);
////                    }
////                }
//            }
            read->set_covered_call(s1_call->block_id, o_side, mendel_pas);
        } else {
//            s1_call->block_id = SPECIFIC_HOMO_BLOCK;
//            if(side == 0 && s1_call->allele1 != s2_call->allele1) {
//                auto t = s1_call->allele1;
//                s1_call->allele1 = s1_call->allele2;
//                s1_call->allele2 = t;
//            } else if(side == 1 && s1_call->allele1 == s2_call->allele1){
//                auto t = s1_call->allele1;
//                s1_call->allele1 = s1_call->allele2;
//                s1_call->allele2 = t;
//            }
        }
    }
//    InfoSet hete_reads;
    infoSet->add_read(read, true);
//    extend(idx1, hete_reads,side);
    if (idx1 == 0) {
        for (auto item : infoSet->confilict_poses) {
            if (side % 2 == 0) {
                conflicts1.insert(item);
            } else {
                conflicts2.insert(item);
            }
        }
    }
}


void ChromoPhaser::check_mendel(int idx1, int idx2, int idx3) {
    for (int i = 0; i < results_for_variant.size(); ++i) {
        auto rec = results_for_variant[i];
        if (rec->pos == 13849720) {
            int tmp = 0;
        }
        auto child = rec->calls[idx1];
        auto father = rec->calls[idx2];
        auto mother = rec->calls[idx3];
        auto c1 = child->allele1;
        auto c2 = child->allele2;
        auto f1 = father->allele1;
        auto f2 = father->allele2;
        auto m1 = mother->allele1;
        auto m2 = mother->allele2;
        if(c1 == -1) {
            mendel_cs.push_back(i);
            continue;
        } else if (f1 == -1 && m1 != -1) {
            if (c1 == m1 || c1 == m2 || c2 == m1 || c2 == m2){
                mendel_passm.push_back(i);
                continue;
            }
        } else if (f1 != -1 && m1 == -1) {
            if (c1 == f1 || c1 == f2 || c2 == f1 || c2 == f2){
                mendel_passf.push_back(i);
                continue;
            }
//            mendel_passf.push_back(i);
//            continue;
        } else if (f1 == -1 && m1 == -1){
            mendel_cs.push_back(i);
        } else {
            if (((c1 == f1 && c2 == m1) || (c1 == f1 && c2 == m2) || (c1 == f2 && c2 == m1) || (c1 == f2 && c2 == m2) ||
                 (c2 == f1 && c1 == m1) || (c2 == f1 && c1 == m2) || (c2 == f2 && c1 == m1) || (c2 == f2 && c1 == m2))){
                mendel_passm.push_back(i);
                mendel_passf.push_back(i);
            }else if (c1 == m1 || c1 == m2 || c2 == m1 || c2 == m2) {
                mendel_passm.push_back(i);
            } else if (c1 == f1 || c1 == f2 || c2 == f1 || c2 == f2) {
                mendel_passf.push_back(i);
            }else
                mendel_cs.push_back(i);
        }
    }
}
