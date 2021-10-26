//
// Created by yonghanyu2 on 10/10/2018.
//

#include "vcf_io.h"
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

VCFReader::VCFReader(const char *filename)
: iter(nullptr), tbx_index(nullptr)
{
    size_t l = strlen(filename);
    if (l < 1)
    {
        std::cerr << "Error: Check input vcf file" << std::endl;
        exit(-1);
    }
    this->filename = new char[l + 1];
    strcpy(this->filename, filename);
    vcf_file = bcf_open(filename, "r");
    if (vcf_file == nullptr)
    {
        std::cerr << "Error: Fail to load VCF file" << std::endl;
    }
    tbx_index = tbx_index_load(filename);
    if (tbx_index == nullptr)
    {
        std::cerr<< "Error: Fail to load tabix file, check whether it exists." << std::endl;
        exit(-1);
    }
    header = bcf_hdr_read(vcf_file);
    if (header == nullptr)
    {
        std::cerr << "Error: corrupted vcf header." << std::endl;
        exit(-1);
    }
    tmp = {0, 0, nullptr};
    get_contigs();
    buffer = bcf_init();
    //count_contigs();
    //reset();
}

VCFReader::~VCFReader()
{
    delete []filename;
    bcf_close(vcf_file);
    bcf_hdr_destroy(header);
    contigs.clear();
    var_count.clear();
    hts_itr_destroy(iter);
    tbx_destroy(tbx_index);
    free(tmp.s);
    bcf_destroy(buffer);
}

template<typename T> struct map_init_helper
{
    T& data;
    explicit map_init_helper(T& d) : data(d) {}
    map_init_helper& operator() (typename T::key_type const& key, typename T::mapped_type const& value)
    {
        data[key] = value;
        return *this;
    }
};

template<typename T> map_init_helper<T> map_init(T& item)
{
    return map_init_helper<T>(item);
}

VCFWriter::VCFWriter(const bcf_hdr_t *hdr, const char *outfile)
{

    fp = bcf_open(outfile, "w");
    map_init(this->FilterN) (filter_type::NOINFO, "INFO_NOT_ENOUGH") (filter_type::POOLRESULT, "POOL_SPECTRAL_RESULT") (filter_type::TENXINCONSISTENCY, "10X_PHASING_INCONSISTENCY") (filter_type::PASS, "PASS");
    map_init(this->FormatN) (GT, "GT") (PS, "PS");

    header = bcf_hdr_dup(hdr);
    sample_count = bcf_hdr_nsamples(header);
    header_init();
    ngt = 2 * sample_count;

    gt = new int[ngt];
}

void VCFWriter::header_init()
{
    char GT[] = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    char PS[] = "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing Block No.\">";

    std::string line_;
    kstring_t kstring = {0, 0, nullptr};
    bcf_hdr_format(header, 0, &kstring);
    char *htxt = new char[kstring.l + 5000];
    char *hdr_txt = ks_release(&kstring);
    std::stringstream ss(hdr_txt);
    int str_l;
    bool Format = false; bool Filter = false;

    std::getline(ss, line_, '\n');
    strcpy(htxt, line_.data()); strcat(htxt, "\n");

    while(std::getline(ss, line_, '\n'))
    {
        if (line_[0] != '#')
            break;
        bcf_hrec_t *hrec = bcf_hdr_parse_line(header, line_.c_str(), &str_l);
        if (hrec == nullptr)
        {
            strcat(htxt, line_.c_str());
            break;
        }
        else
        {
            if (strcmp(hrec->key, "FILTER") == 0)
            {
                if (!Filter)
                {
                    Filter = true;
                    strcat(htxt, "\n");
                }
                if ((hrec->nkeys != 0) && (strcmp(hrec->vals[0], "PASS") == 0))
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
                // format
            else if (strcmp(hrec->key, "FORMAT") == 0)
            {
                if (!Format)
                {
                    Format = true;
                    strcat(htxt, GT);
                    strcat(htxt, "\n");
                    strcat(htxt, PS);
                    strcat(htxt, "\n");
                }
                if (hrec->nkeys != 0 && strcmp(hrec->vals[0], "GT") == 0)
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                if (hrec->nkeys != 0 && strcmp(hrec->vals[0], "PS") == 0)
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
            else
            {
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
        }
        bcf_hrec_destroy(hrec);
    }
    bcf_hdr_destroy(header);
    header = bcf_hdr_init("w");
    bcf_hdr_parse(header, htxt);
    int rr = bcf_hdr_write(fp, header);

    delete []htxt;
    free(hdr_txt);
}

VCFWriter::~VCFWriter()
{
    bcf_hdr_destroy(header);
    bcf_close(fp);
    delete [] gt;
}

void VCFWriter::write_nxt_record(bcf1_t *record, std::shared_ptr<VcfRecord> result, std::vector<int>& ps_nos)
{
    int i, j,nsmpl = bcf_hdr_nsamples(header);
    int *gt_arr = nullptr, *ps_arr = nullptr, ngt_arr =0, nps_arr = 0;
    bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
    bcf_get_format_int32(header, record, "PS", &ps_arr, &nps_arr);

    for (i = 0; i < nsmpl; i++) {
        int ref;
        int alt;
        int *ptr = gt_arr + i * 2;
        int *p_ptr = ps_arr + i;
        auto tcall = result->calls[i];
        int allele0 = bcf_gt_allele(ptr[0]);
        int allele1 = bcf_gt_allele(ptr[1]);
        int temp = bcf_alleles2gt(allele0, allele1);;
        bcf_gt2alleles(temp, &ref, &alt);
        if(tcall->need_flip) {
            if (tcall->phased) {
                gt[2*i] = bcf_gt_phased(ref);
                gt[2*i + 1] = bcf_gt_phased(alt);
            } else {
                gt[2*i] = bcf_gt_unphased(ref);
                gt[2*i + 1] = bcf_gt_unphased(alt);
            }
        } else {
            if (tcall->phased) {
                gt[2*i] = bcf_gt_phased(alt);
                gt[2*i + 1] = bcf_gt_phased(ref);
            } else {
                gt[2*i] = bcf_gt_unphased(alt);
                gt[2*i + 1] = bcf_gt_unphased(ref);
            }
        }
        p_ptr[i] = ps_nos[i];
        bcf_update_genotypes(this->header, record, gt_arr, ngt);
        bcf_update_format_int32(this->header, record, "PS", ps_arr, sample_count);
    }
    free(gt_arr); free(ps_arr);

    bcf1_t *record_w = bcf_dup(record);
    bcf_destroy(record_w);
}



void VCFWriter::write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf)
{
    bcf1_t *record = bcf_init();
    uint blk_count = 0;
    int gap_count = 0;
    std::unordered_map<ptr_PhasedBlock, uint > encountered_phased_block;
    frvcf.jump_to_contig(frvcf.curr_contig);
    std::vector<std::unordered_map<uint, uint>> idx2pos;
    for (int i = 0; i<sample_count;i++) {
        std::unordered_map<uint, uint> tmp;
        idx2pos.push_back(tmp);
    }
    for (uint idx = 0; idx < chromo_phaser->variant_count; idx++)
    {
        frvcf.get_next_record(record);
        auto result = chromo_phaser->results_for_variant[idx];
        bcf_translate(this->header, frvcf.header, record);
        std::vector<int> ps_nos;
        for (int j = 0; j < result->calls.size(); j++) {
            idx2pos[j][idx] = record->pos + 1;
            auto tcall = result->calls[j];
            if (tcall->phased)
            {
                if (gap_count >= 30)
                {
                    blk_count++;
                }
                ptr_PhasedBlock block = tcall->block.lock();
                if (block.get() == nullptr)
                {
                    ps_nos.push_back(0);
                    continue;
                }
                //already met
                if (encountered_phased_block.count(block) != 0)
                {
                    ps_nos.push_back(idx2pos[j][block->start_variant_idx]);
//                    uint blk_no = encountered_phased_block[block];
//                    write_nxt_record(record, result, idx2pos[block->start_variant_idx]);
                }
                else
                {
                    encountered_phased_block.emplace(block, ++blk_count);
//                    write_nxt_record(record, result, idx2pos[block->start_variant_idx]);
                    ps_nos.push_back(idx2pos[j][block->start_variant_idx]);
                }

                gap_count = 0;
            }
            else
            {
                gap_count++;
//                write_nxt_record(record, result, 0);
                ps_nos.push_back(0);
            }
        }
        write_nxt_record(record, result, ps_nos);
    }
    encountered_phased_block.clear();
    bcf_destroy(record);
}

