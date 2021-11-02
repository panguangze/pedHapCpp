//
// Created by yyh on 3/4/2019.
//

#include "phaser.h"
#include "htslib/vcf.h"
#include "type.h"
#include "util.h"
#include <iostream>
#include <cstdlib>
#include <unordered_map>

// TODO clarify between variant count and block count

Phaser::Phaser(const std::string &fnvcf, const std::string &fnout)
{
    frvcf = new VCFReader(fnvcf.data());
    fwvcf = new VCFWriter(frvcf->header, fnout.data());
    sample_count = bcf_hdr_nsamples(this->frvcf->header);
    sample_to_idx = new std::unordered_map<std::string, int>();
    coverage = 30;  //deprecated
    bool use_secondary = false;
    threshold = 1e-5;
    for (int i = 0; i < sample_count; ++i) {
        sample_to_idx->emplace(frvcf->header->samples[i],i);
    }
}

Phaser::~Phaser()
{
    delete frvcf;
    delete fwvcf;
}


//int Phaser::load_contig_records(ChromoPhaser *chromo_phaser)
//{
//    int status = 0;
//    while (true)
//    {
//        VcfRecord result;
//        int status = this->frvcf->get_next_record_contig(result, false);
//        if (status < 0)
//            break;
//        else if (status > 0)
//            continue;
//        chromo_phaser->results_for_variant.push_back(result);
//    }
//
//
//    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
//    for (int i = 0; i < chromo_phaser->variant_count; i++)
//    {
//        chromo_phaser->variant_to_block_id[i] = i;
//    }
//    chromo_phaser->init_block_count = chromo_phaser->variant_count;
//    return status;
//}


int Phaser::load_contig_blocks(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    while (true)
    {
        auto result = std::make_shared<VcfRecord>();
        int status = this->frvcf->get_next_record_contig(result, true);
        if (status < 0) //eof, new chromosome 
            break;
        else if (status > 0) //homo 
            continue;
        chromo_phaser->add_result(result);
    }

    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();

//    uint block_count = 0;
//    std::vector<std::unordered_map<uint, uint>* > ps2block_ids;
//    for (int i = 0; i < chromo_phaser->variant_count; i++)
//    {
//        auto result = chromo_phaser->results_for_variant[i];
//        for(int j = 0; j < result->calls.size(); j++) {
//            if (ps2block_ids.size() <= j) {
//                auto tmp = new std::unordered_map<uint,uint>();
//                ps2block_ids.push_back(tmp);
//            }
//            auto tcall = result->calls[j];
//            uint ps = tcall.ps;
//            if (ps == 0) //not phased
//            {
//                (*(chromo_phaser->variant_to_block_id[i]))[j] = i;
//                block_count++;
//            }
//            else {      //phased
//                if (ps2block_ids[j]->count(ps) == 0)
//                {   // not met before
//                    (*(ps2block_ids[j]))[ps] = i;
//                    (*(chromo_phaser->variant_to_block_id[i]))[j] = (*(ps2block_ids[j]))[ps];
//                    block_count++;
//                }
//                else {
//                    (*(chromo_phaser->variant_to_block_id[i]))[j] = (*(ps2block_ids[j]))[ps];
//                }
//            }
//        }
//    }
//    chromo_phaser->init_block_count = block_count;
    return status;
}


void Phaser::phasing()
{
    for (uint rid = 0; rid < frvcf->contigs_count; rid++)
    {
        if (frvcf->jump_to_contig(rid) != 0)
            break;
        auto nsmp = bcf_hdr_nsamples(this->frvcf->header);
        ChromoPhaser *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], nsmp);
        this->chromoPhaser = chromo_phaser;
        std::string mess = "Reading contig " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess);
        load_contig_blocks(chromo_phaser);
        std::string mess2 = "phasing haplotype for " + std::string(frvcf->contigs[rid]);
        logging(std::clog, mess2);
        std::vector<std::string> names;
        names.emplace_back("child_12");
        names.emplace_back("father_12");
        names.emplace_back("mother_12");
        phasing_by_chrom(names);
        fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf);
        delete chromo_phaser;
    }
}

void Phaser::phasing_with_hete(std::string& s1, std::string& s2, int side, float t1, float t2){
    
}


void Phaser::phasing_by_chrom(std::vector<std::string>& trios)
{
//    while (chromo_phaser->phased->rest_blk_count > 0)
//    {
//        if (chromo_phaser->phased->rest_blk_count > chromo_phaser->init_block_count)
//            break;
//        if (OPERATION == MODE_10X)
//            chromo_phaser->phased->update_phasing_info(MAX_BARCODE_SPANNING);
//        else
//            {
//                if (KEEP_PS)
//                    chromo_phaser->phased->update_phasing_info_keep_phased();
//                else
//                    chromo_phaser->phased->update_phasing_info();
//            }
//        spectral->solver();
//    }
    int s_idx = (*sample_to_idx)[trios[0]];
    int f_idx = (*sample_to_idx)[trios[1]];
    int m_idx = (*sample_to_idx)[trios[2]];
//
//    auto s_blocks = chromo_phaser->phased_blocks_info[s_idx];
//    auto f_blocks = chromo_phaser->phased_blocks_info[f_idx];
//    auto m_blocks = chromo_phaser->phased_blocks_info[m_idx];
//
//    for(auto it : *s_blocks) {
//        auto s_block_start = it.second->start_variant_idx;
//        auto e_block_start = it.second->end_variant_idx;
//
//    }
    logging(std::cerr, "hete1");
    chromoPhaser->phase_with_hete(s_idx, f_idx, 0);
    logging(std::cerr, "hete2");
    chromoPhaser->phase_with_hete(s_idx, m_idx, 1);
    logging(std::cerr, "homo1");
    chromoPhaser->phase_with_homo(s_idx, f_idx,0);
    chromoPhaser->phase_with_homo(s_idx, m_idx,1);
}