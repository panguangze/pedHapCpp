//
// Created by yyh on 12/10/2018.
//
#include "type.h"
#include <iostream>
#include <cmath>
#include <set>
#include "util.h"
#define SPECIFIC_HOMO_BLOCK -123456789l

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
    if (idx2 == -1) return;
    logging(std::cerr,"phasing hete");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    std::unordered_map<uint, PInfo*> reads;
    auto tmp = (idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0) ? this->mendel_passf : this->mendel_passm;
    int current_link_block_id = 2;
//    int prev_s2_block_id = 0;
    int prev_hete_block_id = 0;
    std::unordered_map<int, int> s2_block_id2current_id;
    for(int mendel_pas : tmp) {
        auto result = results_for_variant[mendel_pas];
        if(result->ID == ".") {
//            auto tmp3 = 5;
            continue;
        }
//        if(result->bnd) continue;
        Call* s1_call = result->calls[idx1];
        Call* s2_call = result->calls[idx2];
        if (s1_call->pos == 878665 || s1_call->pos == 878664) {
            int tmp3=1;
        }
        if( s1_call->isHomo() || (!s2_call->isHomo() and !s2_call->isPhased())) continue;
        if (s2_call->isHomo()) {
            s2_call->block_id = prev_hete_block_id;
        } else {
            prev_hete_block_id = s2_call->block_id;
        }
        current_link_block_id = s2_call->block_id;
//        if (s2_call->block_id != prev_s2_block_id) {
//            if (s2_block_id2current_id.find(s2_call->block_id) != s2_block_id2current_id.end()) {
//                current_link_block_id = s2_block_id2current_id[s2_call->block_id];
//                prev_s2_block_id = s2_call->block_id;
//            } else {
//                if (prev_s2_block_id != 0) {
//                    current_link_block_id = current_link_block_id + 1;
//                }
//                prev_s2_block_id = s2_call->block_id;
//                s2_block_id2current_id.emplace(s2_call->block_id,current_link_block_id);
//            }
//        }
//        if (s2_call->isHomo()) {
//        } else {
//
//        }

//        if ((!s2_call->isHomo() && s2_call->block_id != prev_s2_block_id)) {
//            if (s2_block_id2current_id.find(s2_call->block_id) != s2_block_id2current_id.end()) {
//                current_link_block_id = s2_block_id2current_id[s2_call->block_id];
//                prev_s2_block_id = s2_call->block_id;
//            } else {
//                if (prev_s2_block_id != 0) {
//                    current_link_block_id = current_link_block_id + 1;
//                }
//                prev_s2_block_id = s2_call->block_id;
//                s2_block_id2current_id.emplace(s2_call->block_id,current_link_block_id);
//            }
////            if (prev_s2_block_id == 0) continue;
//        } else {
//            if (s2_call->isHomo() && s2_call->block_id != prev_s2_block_id) {
//                current_link_block_id = current_link_block_id + 1;
//                prev_s2_block_id = s2_call->block_id;
//            }
//        }
//        if( s1_call->isHomo() || s2_call->isHomo() || !s2_call->isPhased()) continue;

//        check mendel?
//        if(s1_call->allele1 != s2_call->allele1 && s1_call->allele1 != s2_call->allele2 &&
//            s1_call->allele2 != s2_call->allele1 && s1_call->allele2 != s2_call->allele2)
//            continue;
//        if (s2_call->isHomo())
        if (current_link_block_id == 5) {
            auto mm = 3;
        }
        if(reads.find(current_link_block_id) == reads.end()) {
            auto read = new PInfo(-current_link_block_id);
            reads[current_link_block_id] = read;
        }

        if(s1_call->isPhased()) {
            auto o_side = side;
            if (!s2_call->isHomo()) {
                if(s1_call->allele1 == s2_call->allele1) {
                    o_side = abs(side);
                } else {
                    o_side = abs(side - 1);
                }
                reads[current_link_block_id]->set_covered_call(s1_call->block_id, o_side, mendel_pas, true);
            } else {
//                if s2 is ok
//                if (s2_call->block_id == 1) {
//                    if (side == 0) {
//                        if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
//                            o_side = abs(side -1);
//                        } else {
//                            o_side = abs(side);
//                        }
//                    } else {
//                        if(s1_call->allele2 == s2_call->allele2 || s1_call->allele1 == s2_call->allele1) {
//                            o_side = abs(side -1);
//                        } else {
//                            o_side = abs(side);
//                        }
//                    }
//                } else {
//                    if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
//                        o_side = abs(side -1);
//                    } else {
//                        o_side = abs(side);
//                    }
//                }
//                if(s1_call->allele1 == s2_call->allele2 || s1_call->allele2 == s2_call->allele1) {
//                    o_side = abs(side -1);
//                } else {
//                    o_side = abs(side);
//                }
//                reads[current_link_block_id]->set_covered_call(s1_call->block_id, o_side, mendel_pas, false);
            }
//            auto o_side = side;
        } else {
//            if ((s1_call->allele1 != s2_call->allele2 && s1_call->allele1 != s2_call->allele1) &&
//                    (s1_call->allele2 != s2_call->allele2 && s1_call->allele2 != s2_call->allele1))
//                continue;
//            if (s2_call->isHomo()) continue;
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
        if (it.first == 1) {
            auto mmm = 33;
            auto s0 = it.second->side0_support[194418247];
            auto s1 = it.second->side1_support[194418247];
            auto mmmm = 44;
        }
        if (it.first == 5) {
            int lkk = 9;
        }
        if(it.second->blocks.size() <= 1) continue;
//        extract_lst(idx2, idx1, it.second, results_for_variant, this->prev_contig_variant_count);
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
    for (auto it: reads) {
        free(it.second);
    }
    reads.clear();
}

void ChromoPhaser::extend(int idx, InfoSet* infoSet, int side, int type) {
    std::unordered_map<int ,int> finalize_new_block_ids;
    int f_new_id;
    for(int i = 0; i < this->results_for_variant.size(); i++) {
        auto result = results_for_variant[i];
        if (result->pos == 23343) {
            int tmp=1;
        }
        Call* s1_call = result->calls[idx];
        Call* s2_call;
        Call* s3_call;
        if (idx == 0) {
            s2_call = result->calls[idx+1];
            s3_call = result->calls[idx+2];
        }
        if (s1_call->isHomo()) continue;
        if (s1_call->block_id == 0) continue;
        if (s1_call->block_id == SPECIFIC_HOMO_BLOCK) {
            s1_call->block_id = 1;
            continue;
        }
        if (result->pos == 23343) {
            int tmp = 2;
            auto v = infoSet->blocks_reverse_info.find(10) == infoSet->blocks_reverse_info.end();
            int tp = 2;
        }
        auto new_block_id = infoSet->find(s1_call->block_id);
        auto n_r = infoSet->blocks_reverse_info[s1_call->block_id];
        if (new_block_id == -1) {
            if (finalize_new_block_ids.find(s1_call->block_id) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[s1_call->block_id];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                auto const result = finalize_new_block_ids.insert(std::make_pair(s1_call->block_id, f_new_id));
                if (not result.second) { result.first->second = f_new_id; }
//                finalize_new_block_ids.emplace(s1_call->block_id, f_new_id);
            }
        } else {
            if(finalize_new_block_ids.find(new_block_id) != finalize_new_block_ids.end()) {
                f_new_id = finalize_new_block_ids[new_block_id];
            } else {
                f_new_id = finalize_new_block_ids.size() + 1;
                auto const result = finalize_new_block_ids.insert(std::make_pair(new_block_id, f_new_id));
                if (not result.second) { result.first->second = f_new_id; }
//                finalize_new_block_ids.emplace(new_block_id, f_new_id);
            }
        }
//        if (f_new_id == 1) {
//            if(type == 0) {
//                s1_call->block_id = f_new_id + 1;
//            } else if ((s2_call->isHomo() || s3_call->isHomo())){
//                s1_call->block_id = f_new_id;
//            }
//        } else {
//            s1_call->block_id = f_new_id + 1;
//        }
        if(type == 0) {
            s1_call->block_id = f_new_id + 1;

        } else if (new_block_id !=-1){
            s1_call->block_id = 1;
        }
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
    if (idx2 == -1) return;
    logging(std::cerr,"phasing homo");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    auto tmp = ((idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0)) ? this->mendel_passf : this->mendel_passm;

    for(int i = 0 ; i < tmp.size(); i++) {

        auto mendel_pas = tmp[i];
        auto result = results_for_variant[mendel_pas];
        if(result->ID == ".") {
//            auto tmp3 = 5;
            continue;
        }
        if (result->pos == 23343) {
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
                    o_side = abs(side - 1);
                } else {
                    o_side = abs(side);
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
            read->set_covered_call(s1_call->block_id, o_side, mendel_pas, false);
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
//    extract_lst(idx2, idx1, read, results_for_variant, this->prev_contig_variant_count);
    auto tttm = infoSet->blocks_reverse_info[92];
    auto tmpjj = read->side0_support[10];
    auto tmpjj2 = read->side1_support[10];
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
    free(read);
}


void ChromoPhaser::phase_with_homo2(int idx1, int idx2, int side, InfoSet* infoSet) {
    logging(std::cerr,"phasing homo");
//    auto idx2_phase_block = this->phased_blocks_info[idx2];
    auto read = new PInfo(SPECIFIC_HOMO_BLOCK);
    auto tmp = ((idx1 == 0 && idx2 == 1) or (idx1 == 1 && idx2 == 0)) ? this->mendel_passf : this->mendel_passm;

    for(int i = 0 ; i < tmp.size(); i++) {
        auto mendel_pas = tmp[i];
        auto result = results_for_variant[mendel_pas];
        if(result->ID == ".") {
//            auto tmp3 = 5;
            continue;
        }
        if (result->pos == 481589) {
            int tmp3=1;
        }
        Call *s1_call = result->calls[idx1];
        Call *s2_call = result->calls[idx2];
        if(s1_call->block_id == 10) {
            int mmm = 33;
        }
//        if (!) {
//
//        }
//        if (s1_call->isHomo() || (!s2_call->isHomo() )) continue;
        if ((!s2_call->isPhased() && !s2_call->isHomo()) || (!s2_call->isHomo() && s2_call->isPhased() && s2_call->block_id != 1) || s1_call->isHomo()) continue;
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
                } else if (s2_call->allele1 == s1_call->allele2) {
                    o_side = abs(side - 1);
                }
            } else {
                if (s2_call->allele2 == s1_call->allele1) {
                    o_side = abs(side - 1);
                } else if (s2_call->allele2 == s1_call->allele2) {
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
            if(s1_call->block_id == 9) {
                int tmp4 = 1;
            }
            read->set_covered_call(s1_call->block_id, o_side, mendel_pas, false);
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
    free(read);
}


void ChromoPhaser::check_mendel(int idx1, int idx2, int idx3) {
    if (!mendel_cs.empty()) return;
    for (int i = 0; i < results_for_variant.size(); ++i) {
        if (idx2 == -1 || idx3 == -1) {
            mendel_cs.push_back(i);
            mendel_passm.push_back(i);
            mendel_passf.push_back(i);
            continue;
        }
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
        if (c1 == -1 || c1 == -1 || f1 == -1 || f2 == -1 || m1 == -1 || m2 == -1) {
            mendel_cs.push_back(i);
            continue;
        } else {
            if (((c1 == f1 && c2 == m1) || (c1 == f1 && c2 == m2) || (c1 == f2 && c2 == m1) || (c1 == f2 && c2 == m2) ||
                 (c2 == f1 && c1 == m1) || (c2 == f1 && c1 == m2) || (c2 == f2 && c1 == m1) || (c2 == f2 && c1 == m2))){
                mendel_passm.push_back(i);
                mendel_passf.push_back(i);
        }
//        if(c1 == -1) {
//            mendel_cs.push_back(i);
//            continue;
//        } else if (f1 == -1 && m1 != -1) {
//            if (c1 == m1 || c1 == m2 || c2 == m1 || c2 == m2){
//                mendel_passm.push_back(i);
//                continue;
//            }
//        } else if (f1 != -1 && m1 == -1) {
//            if (c1 == f1 || c1 == f2 || c2 == f1 || c2 == f2){
//                mendel_passf.push_back(i);
//                continue;
//            }
////            mendel_passf.push_back(i);
////            continue;
//        } else if (f1 == -1 && m1 == -1){
//            mendel_cs.push_back(i);
//        } else {
//            if (((c1 == f1 && c2 == m1) || (c1 == f1 && c2 == m2) || (c1 == f2 && c2 == m1) || (c1 == f2 && c2 == m2) ||
//                 (c2 == f1 && c1 == m1) || (c2 == f1 && c1 == m2) || (c2 == f2 && c1 == m1) || (c2 == f2 && c1 == m2))){
//                mendel_passm.push_back(i);
//                mendel_passf.push_back(i);
//            }else if (c1 == m1 || c1 == m2 || c2 == m1 || c2 == m2) {
//                mendel_passm.push_back(i);
//            } else if (c1 == f1 || c1 == f2 || c2 == f1 || c2 == f2) {
//                mendel_passf.push_back(i);
//            }else
//                mendel_cs.push_back(i);
        }
    }
}

void extract_lst(int pos,int idx, PInfo* it,std::vector<std::shared_ptr<VcfRecord>>& result_for_variant, int prev_count) {

    int lst_count = 1;
    int prev_total_lst = 0;
    std::string lst1;
    std::string qual_str;
    if (it->side1_support.size() >20) {
        lst1.append("20");
    } else{
        lst1.append(std::to_string(it->side1_support.size()));
    }
    lst1.append(" ");
    lst1.append(std::to_string(pos));
    lst1.append("xxxx ");
    int final_sub = it->side1_support.size() % 20;
//    lst1.append()
    std::string prev_pos;
    std::string prev_hap;
    char prev_qual;
    for (auto item: it->side0_support) {

        if (lst_count > 20) {
            lst1.append(qual_str);
            lst1.append(" ");
            lst1.append(std::to_string(60));
            std::cout << lst1 << std::endl;
            lst1 = "";
            qual_str = "";
            prev_total_lst += 20;
            if ( it->side1_support.size() - prev_total_lst >= 20) {
                lst1.append("21");
            } else {
                lst1.append(std::to_string(it->side1_support.size() - prev_total_lst + 1));
            }
            lst1.append(" ");
            lst1.append(std::to_string(pos));
            lst1.append("xxxx ");
            if (prev_hap != "") {
                lst1.append(prev_pos);
                lst1.append(" ");
                lst1.append(prev_hap);
                lst1.append(" ");
                qual_str = prev_qual;
            }
            lst_count = 1;
        }
            auto item2_second = it->side1_support[item.first];
            auto flag = item.second->size() >= item2_second->size();
            auto current_item = flag ? item.second : item2_second;
//        auto prev_pos = 0;
            if(current_item == nullptr) continue;
            if (current_item->empty()) continue;
            if ((*current_item)[0] == 235) {
                int jjj=33;
            }
            lst1.append(std::to_string((*current_item)[0] + prev_count + 1));
            lst1.append(" ");
            auto call = result_for_variant[(*current_item)[0]]->calls[idx];
            if (flag) {
                lst1.append(std::to_string(call->allele1));
            } else {
                lst1.append(std::to_string(call->allele2));
            }
            lst1.append(" ");
//        lst1.append(std::string((*current_item->front()));
            char c = 50 + current_item->size();
            if (c >= 60)
                c = '<';
            qual_str = qual_str + c;
            if (lst_count == 20) {
                prev_pos = std::to_string((*current_item)[0] + prev_count + 1);
                prev_hap = flag ? std::to_string(call->allele1):std::to_string(call->allele2);
                prev_qual = c;
            }


        lst_count += 1;

//        qual_str.append(std::to_string(char(c)));
    }
    lst1.append(qual_str);
    lst1.append(" ");
    lst1.append(std::to_string(60));
    std::cout << lst1 << std::endl;
//        if (lst2_count != 0)
//            std::cout<<lst2_count<<" "<<lst2<<std::endl;


//    std::string lst1;
//    lst1.append(" "+std::to_string(pos));
//    if (pos == 21056) {
//        int lkk = 9;
//    }
//    int lst1_count = 0;
//    int l1_variant_count = 0;
//    for (auto item: it->side0_support) {
//        auto item2 = it->side1_support[item.first];
//        auto prev_pos = 0;
//        std::string prev_str;
//        if(item.second == nullptr) continue;
//        if (item.second->empty()) continue;
//        l1_variant_count = l1_variant_count + item.second->size();
//        for (auto i : *(item).second){
//            if (i - prev_pos != 1) {
//                lst1.append(prev_str + " ");
//                prev_str = "";
//                prev_str.append(std::to_string(i)+" ");
//                lst1_count++;
//            }
//            prev_pos = i;
//            auto call = result_for_variant[prev_pos]->calls[idx];
//            prev_str.append(std::to_string(call->allele1));
//        }
//        lst1.append(prev_str);
//    }
//    lst1.append(" "+ std::string(l1_variant_count,'G')+" "+ std::to_string(60));
//
//
//
//    std::string lst2;
//    int l2_variant_count = 0;
//    lst2.append(" "+std::to_string(pos));
//    int lst2_count = 0;
//    for (auto item: it->side1_support) {
//        auto prev_pos = 0;
//        std::string prev_str;
//        if(item.second == nullptr) continue;
//        if (item.second->empty()) continue;
//        l2_variant_count = l2_variant_count + item.second->size();
//        for (auto i : *(item).second){
//            if (i - prev_pos != 1) {
//                lst2.append(prev_str + " ");
//                prev_str = "";
//                prev_str.append(std::to_string( i)+" ");
//                lst2_count++;
//            }
//            prev_pos = i;
//            auto call = result_for_variant[prev_pos]->calls[idx];
//            prev_str.append(std::to_string(call->allele2));
//        }
//        lst2.append(prev_str);
////        lst2.append(" "+ std::string(item.second->size(),'G')+" "+ std::to_string(60));
//    }
//    lst2.append(" "+ std::string(l2_variant_count,'G')+" "+ std::to_string(60));
//
////    if (lst1_count != 0)
////        std::cout<<lst1_count<<" "<<lst1<<std::endl;
////    if (lst2_count != 0)
////        std::cout<<lst2_count<<" "<<lst2<<std::endl;
//    if (l1_variant_count >= l2_variant_count) {
//        if (lst1_count != 0)
//            std::cout<<lst1_count<<" "<<lst1<<std::endl;
//    } else {
//        if (lst2_count != 0)
//            std::cout<<lst2_count<<" "<<lst2<<std::endl;
//    }
}
