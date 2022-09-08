//
// Created by yyh on 3/4/2019.
//

#include "chromo_phaser.h"
#include "htslib/vcf.h"
#include "type.h"
#include "util.h"
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <string>

// TODO clarify between variant count and block count

Phaser::Phaser(const std::string &fnvcf, const std::string &fnout, const std::string &frecon)
{
    frvcf = new VCFReader(fnvcf.data());
    fwvcf = new VCFWriter(frvcf->header, fnout.data(), frecon.data());
    frvcf->sample_count = bcf_hdr_nsamples(this->frvcf->header);
    frvcf->sample_to_idx = new std::unordered_map<std::string, int>();
    for (int i = 0; i < frvcf->sample_count; ++i) {
        frvcf->sample_to_idx->emplace(frvcf->header->samples[i],i);
    }
}

Phaser::~Phaser()
{
    delete frvcf;
    delete fwvcf;
//    delete sample_to_idx;
}


int Phaser::load_contig_blocks(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    std::vector<int> prev_block_ids;
    prev_block_ids.reserve(chromo_phaser->sample_count);
    for(int i = 0; i < chromo_phaser->sample_count; i++) {
            prev_block_ids.push_back(0);
    }
    while (true)
    {
        auto result = std::make_shared<VcfRecord>();
        for(int i = 0; i < chromo_phaser->sample_count; i++) {
            auto item = result->calls[i];
            if(item->isHomo()) {
                item->block_id = prev_block_ids[i];
            } else {
                prev_block_ids[i] = item->block_id;
            }
        }
        status = this->frvcf->get_next_record_contig(result, true);
        if (status < 0) //eof, new chromosome 
            break;
        else if (status > 0) //homo 
            continue;
        if (result->pos == 13849720) {
            int tmp = 0;
        }
        chromo_phaser->add_result(result);
    }
    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
    return status;
}


void Phaser::phasing()
{
    auto nsmp = bcf_hdr_nsamples(this->frvcf->header);
    uint prev_count = 0;
    for (uint rid = 0; rid < frvcf->contigs_count; rid++)
    {
        if (!this->contigs.empty() && std::find(this->contigs.begin(), this->contigs.end(), frvcf->contigs[rid]) == this->contigs.end()) {
            if (frvcf->jump_to_contig(rid) != 0)
                break;
            auto *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], nsmp);
            load_contig_blocks(chromo_phaser);
            prev_count = prev_count + chromo_phaser->variant_count;
            continue;
        }
//            continue;
        if (frvcf->jump_to_contig(rid) != 0)
            break;
        auto *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], nsmp);
        this->chromoPhaser = chromo_phaser;
        this->chromoPhaser->set_prev_contig_variant_count(prev_count);
        std::string mess = "Reading contig " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess);
        load_contig_blocks(chromo_phaser);
        std::string mess2 = "phasing haplotype for " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess2);
        phasing_by_chrom();
        fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf, prev_count);
        prev_count = prev_count + chromo_phaser->variant_count;
        delete chromo_phaser;
    }
}

std::vector<std::vector<std::string>> Phaser::get_down_to_up() const {
    return this->frvcf->up_to_down;
}

std::vector<std::vector<std::string>> Phaser::get_up_to_down() const {
    return this->frvcf->down_to_up;
}

std::unordered_map<std::string, int>* Phaser::get_sample_to_idx() const {
    return frvcf->sample_to_idx;
}

//void Phaser::phasing_xy() const {
//    int s_idx=-1, f_idx=-1, m_idx=-1;
//    for(auto it: get_up_to_down()) {
//        s_idx = (*get_sample_to_idx())[it[0]];
//        if(it[1] != EMPTY_ID)
//            f_idx = (*get_sample_to_idx())[it[1]];
//        if(it[2] != EMPTY_ID)
//            m_idx = (*get_sample_to_idx())[it[2]];
//        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
////        chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//        chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
//        chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
//        chromoPhaser->phase_with_homo(s_idx, f_idx,0);
//        chromoPhaser->phase_with_homo(s_idx, m_idx,1);
//    }
//    if(ONLY_CHILD) return;
//    for(auto it: get_down_to_up()) {
//        s_idx = (*get_sample_to_idx())[it[0]];
////        chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//        if(it[1] != EMPTY_ID) f_idx = (*get_sample_to_idx())[it[1]];
//        if(it[2] != EMPTY_ID) m_idx = (*get_sample_to_idx())[it[2]];
//        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//        if(f_idx != -1) {
//            chromoPhaser->phase_with_hete(f_idx, s_idx, 0);
//            chromoPhaser->phase_with_homo(f_idx, s_idx,0);
//        }
//        if (m_idx != -1) {
//            chromoPhaser->phase_with_hete(m_idx, s_idx, 0);
//            chromoPhaser->phase_with_homo(m_idx, s_idx,0);
//        }
//    }
//
//}

void Phaser::correct() const {

}

void Phaser::phasing_by_chrom() const
{
    int s_idx=-1, f_idx=-1, m_idx=-1;
    bool is_child_male;
    int i = 0;
    while (i != 1) {
        InfoSet* hete_reads;
        InfoSet* home_reads;
        for(auto it: get_up_to_down()) {
            s_idx = (*get_sample_to_idx())[it[0]];
            if(it[1] != EMPTY_ID)
                f_idx = (*get_sample_to_idx())[it[1]];
            if(it[2] != EMPTY_ID)
                m_idx = (*get_sample_to_idx())[it[2]];
            is_child_male = it[3] == "male";
//        male and y
            hete_reads = new InfoSet();
            home_reads = new InfoSet();
            if(is_child_male && chromoPhaser->is_y()) {
                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
            } else if (is_child_male && chromoPhaser->is_x()) { //male and x
                chromoPhaser->phase_with_hete(s_idx, m_idx, 1,hete_reads);
                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
            } else if (!is_child_male && chromoPhaser->is_y()) continue; // female and y continue.
            else {
                if(s_idx != -1 && m_idx != -1 && f_idx != -1)
                    chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
                chromoPhaser->phase_with_hete(s_idx, m_idx, 1,hete_reads);
                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
                chromoPhaser->extend(s_idx,hete_reads,0, 0);
//                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
//                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
//                chromoPhaser->extend(s_idx,home_reads,0,1);
            }

//            chromoPhaser->extend(s_idx,hete_reads,0);
            free(hete_reads);
            free(home_reads);
        }
        if(ONLY_CHILD) return;
        for(auto it: get_down_to_up()) {
            s_idx = (*get_sample_to_idx())[it[0]];
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            if(it[1] != EMPTY_ID) f_idx = (*get_sample_to_idx())[it[1]];
            if(it[2] != EMPTY_ID) m_idx = (*get_sample_to_idx())[it[2]];
            is_child_male = it[3] == "male";
//        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            if(f_idx != -1) {
                if (!is_child_male && chromoPhaser->is_y()) continue;
                hete_reads = new InfoSet();
                home_reads = new InfoSet();
                chromoPhaser->phase_with_hete(f_idx, s_idx, 0, hete_reads);
                chromoPhaser->extend(f_idx,hete_reads,0, 0);
//                chromoPhaser->phase_with_homo(f_idx, s_idx,0, home_reads);
//                chromoPhaser->extend(f_idx,home_reads,0,1);
                free(hete_reads);
                free(home_reads);
            }
            if (m_idx != -1) {
                if (chromoPhaser->is_y()) continue;
                hete_reads = new InfoSet();
                home_reads = new InfoSet();
                chromoPhaser->phase_with_hete(m_idx, s_idx, 0, hete_reads);
                chromoPhaser->extend(m_idx,hete_reads,0, 0);
//                chromoPhaser->phase_with_homo(m_idx, s_idx,0, home_reads);
//                chromoPhaser->extend(m_idx,home_reads,0,1);
                free(hete_reads);
                free(home_reads);
            }
        }
        i++;
    }

    i = 0;

    while (i != 0) {
        InfoSet* hete_reads;
        InfoSet* home_reads;
        for(auto it: get_up_to_down()) {
            s_idx = (*get_sample_to_idx())[it[0]];
            if(it[1] != EMPTY_ID)
                f_idx = (*get_sample_to_idx())[it[1]];
            if(it[2] != EMPTY_ID)
                m_idx = (*get_sample_to_idx())[it[2]];
            is_child_male = it[3] == "male";
//        male and y
            hete_reads = new InfoSet();
            home_reads = new InfoSet();
            if(is_child_male && chromoPhaser->is_y()) {
                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
            } else if (is_child_male && chromoPhaser->is_x()) { //male and x
                chromoPhaser->phase_with_hete(s_idx, m_idx, 1,hete_reads);
                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
            } else if (!is_child_male && chromoPhaser->is_y()) continue; // female and y continue.
            else {
                if(s_idx != -1 && m_idx != -1 && f_idx != -1)
                    chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
//                chromoPhaser->phase_with_hete(s_idx, m_idx, 1,hete_reads);
                chromoPhaser->extend(s_idx,hete_reads,0,0);
//                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
//                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
//                chromoPhaser->extend(s_idx,home_reads,0);
            }

//            chromoPhaser->extend(s_idx,hete_reads,0);
            free(hete_reads);
            free(home_reads);
        }
        if(ONLY_CHILD) return;
        for(auto it: get_down_to_up()) {
            s_idx = (*get_sample_to_idx())[it[0]];
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            if(it[1] != EMPTY_ID) f_idx = (*get_sample_to_idx())[it[1]];
            if(it[2] != EMPTY_ID) m_idx = (*get_sample_to_idx())[it[2]];
            is_child_male = it[3] == "male";
//        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            if(f_idx != -1) {
                if (!is_child_male && chromoPhaser->is_y()) continue;
                hete_reads = new InfoSet();
                home_reads = new InfoSet();
//                chromoPhaser->phase_with_hete(f_idx, s_idx, 0, hete_reads);
//                chromoPhaser->extend(f_idx,hete_reads,0,0);
//                chromoPhaser->phase_with_homo(f_idx, s_idx,0, home_reads);
//                chromoPhaser->extend(f_idx,home_reads,0);
                free(hete_reads);
                free(home_reads);
            }
            if (m_idx != -1) {
                if (chromoPhaser->is_y()) continue;
                hete_reads = new InfoSet();
                home_reads = new InfoSet();
//                chromoPhaser->phase_with_hete(m_idx, s_idx, 0, hete_reads);
//                chromoPhaser->extend(m_idx,hete_reads,0,0);
//                chromoPhaser->phase_with_homo(m_idx, s_idx,0, home_reads);
//                chromoPhaser->extend(m_idx,home_reads,0);
                free(hete_reads);
                free(home_reads);
            }
        }

//        for(auto it: get_up_to_down()) {
//            s_idx = (*get_sample_to_idx())[it[0]];
//            if(it[1] != EMPTY_ID)
//                f_idx = (*get_sample_to_idx())[it[1]];
//            if(it[2] != EMPTY_ID)
//                m_idx = (*get_sample_to_idx())[it[2]];
//            is_child_male = it[3] == "male";
////        male and y
//            hete_reads = new InfoSet();
//            home_reads = new InfoSet();
//            if(is_child_male && chromoPhaser->is_y()) {
//                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
//                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
//            } else if (is_child_male && chromoPhaser->is_x()) { //male and x
//                chromoPhaser->phase_with_hete(s_idx, m_idx, 0,hete_reads);
//                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
//            } else if (!is_child_male && chromoPhaser->is_y()) continue; // female and y continue.
//            else {
////                if(s_idx != -1 && m_idx != -1 && f_idx != -1)
////                    chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
////                chromoPhaser->phase_with_hete(s_idx, f_idx, 0,hete_reads);
////                chromoPhaser->phase_with_hete(s_idx, m_idx, 1,hete_reads);
////                chromoPhaser->extend(s_idx,hete_reads,0);
//                chromoPhaser->phase_with_homo(s_idx, f_idx,0, home_reads);
//                chromoPhaser->phase_with_homo(s_idx, m_idx,1, home_reads);
//                chromoPhaser->extend(s_idx,home_reads,0,1);
//                chromoPhaser->correct_conflict(s_idx);
//            }
////            chromoPhaser->extend(s_idx,hete_reads,0);
//            free(hete_reads);
//            free(home_reads);
//        }
//
//        if(ONLY_CHILD) return;
//        for(auto it: get_down_to_up()) {
//            s_idx = (*get_sample_to_idx())[it[0]];
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//            if(it[1] != EMPTY_ID) f_idx = (*get_sample_to_idx())[it[1]];
//            if(it[2] != EMPTY_ID) m_idx = (*get_sample_to_idx())[it[2]];
//            is_child_male = it[3] == "male";
////            if(s_idx != -1 && m_idx != -1 && f_idx != -1)
////                chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//            if(f_idx != -1) {
//                if (!is_child_male && chromoPhaser->is_y()) continue;
//                hete_reads = new InfoSet();
//                home_reads = new InfoSet();
////                chromoPhaser->phase_with_hete(f_idx, s_idx, 0, hete_reads);
////                chromoPhaser->extend(f_idx,hete_reads,0);
//                chromoPhaser->phase_with_homo2(f_idx, s_idx,0, home_reads);
//                chromoPhaser->extend(f_idx,home_reads,0,1);
//                free(hete_reads);
//                free(home_reads);
//            }
//            if (m_idx != -1) {
//                if (chromoPhaser->is_y()) continue;
//                hete_reads = new InfoSet();
//                home_reads = new InfoSet();
////                chromoPhaser->phase_with_hete(m_idx, s_idx, 1, hete_reads);
////                chromoPhaser->extend(m_idx,hete_reads,0);
//                chromoPhaser->phase_with_homo2(m_idx, s_idx,1, home_reads);
//                chromoPhaser->extend(m_idx,home_reads,0,1);
//                free(hete_reads);
//                free(home_reads);
//            }
//        }

        i++;
    }

//
//    for(auto it: get_up_to_down()) {
//        s_idx = (*get_sample_to_idx())[it[0]];
//        if(it[1] != EMPTY_ID)
//            f_idx = (*get_sample_to_idx())[it[1]];
//        if(it[2] != EMPTY_ID)
//            m_idx = (*get_sample_to_idx())[it[2]];
//        is_child_male = it[3] == "male";
////        male and y
//        if(is_child_male && chromoPhaser->is_y()) {
//            chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
//            chromoPhaser->phase_with_homo(s_idx, f_idx,0);
//        } else if (is_child_male && chromoPhaser->is_x()) { //male and x
//            chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
//            chromoPhaser->phase_with_homo(s_idx, m_idx,1);
//        } else if (!is_child_male && chromoPhaser->is_y()) continue; // female and y continue.
//        else {
////            if(s_idx != -1 && m_idx != -1 && f_idx != -1)
////                chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//            chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
//            chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
//            chromoPhaser->phase_with_homo(s_idx, f_idx,0);
//            chromoPhaser->phase_with_homo(s_idx, m_idx,1);
//            chromoPhaser->correct_conflict(s_idx);
//        }
//    }
////
//    for(auto it: get_down_to_up()) {
//        s_idx = (*get_sample_to_idx())[it[0]];
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//        if(it[1] != EMPTY_ID) f_idx = (*get_sample_to_idx())[it[1]];
//        if(it[2] != EMPTY_ID) m_idx = (*get_sample_to_idx())[it[2]];
//        is_child_male = it[3] == "male";
////        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
////            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
//        if(f_idx != -1) {
//            if (!is_child_male && chromoPhaser->is_y()) continue;
//            chromoPhaser->phase_with_hete(f_idx, s_idx, 0);
//            chromoPhaser->phase_with_homo(f_idx, s_idx,0);
//        }
//        if (m_idx != -1) {
//            if (chromoPhaser->is_y()) continue;
//            chromoPhaser->phase_with_hete(m_idx, s_idx, 0);
//            chromoPhaser->phase_with_homo(m_idx, s_idx,0);
//        }
//    }
}