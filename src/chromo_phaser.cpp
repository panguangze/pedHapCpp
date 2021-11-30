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
#define EMPTY_ID "null"

// TODO clarify between variant count and block count

Phaser::Phaser(const std::string &fnvcf, const std::string &fnout)
{
    frvcf = new VCFReader(fnvcf.data());
    fwvcf = new VCFWriter(frvcf->header, fnout.data());
    sample_count = bcf_hdr_nsamples(this->frvcf->header);
    sample_to_idx = new std::unordered_map<std::string, int>();
    for (int i = 0; i < sample_count; ++i) {
        sample_to_idx->emplace(frvcf->header->samples[i],i);
    }
}

Phaser::~Phaser()
{
    delete frvcf;
    delete fwvcf;
    delete sample_to_idx;
}


int Phaser::load_contig_blocks(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    while (true)
    {
        auto result = std::make_shared<VcfRecord>();
        status = this->frvcf->get_next_record_contig(result, true);
        if (status < 0) //eof, new chromosome 
            break;
        else if (status > 0) //homo 
            continue;
        chromo_phaser->add_result(result);
    }
    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
    return status;
}


void Phaser::phasing()
{
    for (uint rid = 0; rid < frvcf->contigs_count; rid++)
    {
        if (frvcf->jump_to_contig(rid) != 0)
            break;
        auto nsmp = bcf_hdr_nsamples(this->frvcf->header);
        auto *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], nsmp);
        this->chromoPhaser = chromo_phaser;
        std::string mess = "Reading contig " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess);
        load_contig_blocks(chromo_phaser);
        std::string mess2 = "phasing haplotype for " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess2);
        phasing_by_chrom();
        fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf);
        delete chromo_phaser;
    }
}


void Phaser::phasing_xy() const {
    int s_idx=-1, f_idx=-1, m_idx=-1;
    for(auto it: up_to_down) {
        s_idx = (*sample_to_idx)[it[0]];
        if(it[1] != EMPTY_ID)
            f_idx = (*sample_to_idx)[it[1]];
        if(it[2] != EMPTY_ID)
            m_idx = (*sample_to_idx)[it[2]];
        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
        chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
        chromoPhaser->phase_with_homo(s_idx, f_idx,0);
        chromoPhaser->phase_with_homo(s_idx, m_idx,1);
    }
    if(ONLY_CHILD) return;
    for(auto it: down_to_up) {
        s_idx = (*sample_to_idx)[it[0]];
        chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        if(it[1] != EMPTY_ID) f_idx = (*sample_to_idx)[it[1]];
        if(it[2] != EMPTY_ID) m_idx = (*sample_to_idx)[it[2]];
        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        if(f_idx != -1) {
            chromoPhaser->phase_with_hete(f_idx, s_idx, 0);
            chromoPhaser->phase_with_homo(f_idx, s_idx,0);
        }
        if (m_idx != -1) {
            chromoPhaser->phase_with_hete(m_idx, s_idx, 0);
            chromoPhaser->phase_with_homo(m_idx, s_idx,0);
        }
    }

}

void Phaser::phasing_by_chrom() const
{
    int s_idx=-1, f_idx=-1, m_idx=-1;
    bool is_child_male;
    for(auto it: up_to_down) {
        s_idx = (*sample_to_idx)[it[0]];
        if(it[1] != EMPTY_ID)
            f_idx = (*sample_to_idx)[it[1]];
        if(it[2] != EMPTY_ID)
            m_idx = (*sample_to_idx)[it[2]];
        is_child_male = it[3] == "male";
//        male and y
        if(is_child_male && chromoPhaser->is_y()) {
            chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
            chromoPhaser->phase_with_homo(s_idx, f_idx,0);
        } else if (is_child_male && chromoPhaser->is_x()) { //male and x
            chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
            chromoPhaser->phase_with_homo(s_idx, m_idx,1);
        } else if (!is_child_male && chromoPhaser->is_y()) continue; // female and y continue.
        else {
            if(s_idx != -1 && m_idx != -1 && f_idx != -1)
                chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
            chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
            chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
            chromoPhaser->phase_with_homo(s_idx, f_idx,0);
            chromoPhaser->phase_with_homo(s_idx, m_idx,1);
        }
    }
    if(ONLY_CHILD) return;
    for(auto it: down_to_up) {
        s_idx = (*sample_to_idx)[it[0]];
//            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        if(it[1] != EMPTY_ID) f_idx = (*sample_to_idx)[it[1]];
        if(it[2] != EMPTY_ID) m_idx = (*sample_to_idx)[it[2]];
        is_child_male = it[3] == "male";
        if(s_idx != -1 && m_idx != -1 && f_idx != -1)
            chromoPhaser->check_mendel(s_idx, f_idx, m_idx);
        if(f_idx != -1) {
            if (!is_child_male && chromoPhaser->is_y()) continue;
            chromoPhaser->phase_with_hete(f_idx, s_idx, 0);
            chromoPhaser->phase_with_homo(f_idx, s_idx,0);
        }
        if (m_idx != -1) {
            if (chromoPhaser->is_y()) continue;
            chromoPhaser->phase_with_hete(m_idx, s_idx, 0);
            chromoPhaser->phase_with_homo(m_idx, s_idx,0);
        }
    }
}