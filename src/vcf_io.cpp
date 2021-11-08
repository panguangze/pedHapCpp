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
    map_init(this->FormatN) (GT, "GT") (PS, "PS");

    header = bcf_hdr_dup(hdr);
    sample_count = bcf_hdr_nsamples(header);
    header_init();
    ngt = 2 * sample_count;

    gt = new int[ngt];
    ps = new int[sample_count];
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
//    int *gt_arr = nullptr, *ps_arr = nullptr, ngt_arr =0, nps_arr = 0;
//    bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
//    bcf_get_format_int32(header, record, "PS", &ps_arr, &nps_arr);
//    if (ps_arr == nullptr) {
//        ps_arr = (int *) malloc(ngt/2);
//    }

    for (i = 0; i < nsmpl; i++) {
        auto tcall = result->calls[i];
        if (tcall->isPhased()) {
            gt[2*i] = bcf_gt_phased(tcall->allele1);
            gt[2*i + 1] = bcf_gt_phased(tcall->allele2);
        } else {
            gt[2*i] = bcf_gt_unphased(tcall->allele1);
            gt[2*i + 1] = bcf_gt_unphased(tcall->allele2);
        }
        ps[i] = ps_nos[i];
    }
    bcf_update_genotypes(this->header, record, gt, ngt);
    bcf_update_format_int32(this->header, record, "PS", ps, sample_count);
    bcf1_t *record_w = bcf_dup(record);
    int rr = bcf_write(this->fp, this->header, record_w);
    bcf_destroy(record_w);
}



void VCFWriter::write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf)
{
    bcf1_t *record = bcf_init();
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
        for (auto tcall : result->calls) {
            ps_nos.push_back(tcall->block_id);
        }
        write_nxt_record(record, result, ps_nos);
    }
    bcf_destroy(record);
}

