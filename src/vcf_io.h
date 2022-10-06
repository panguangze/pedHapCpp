//
// Created by yonghanyu2 on 10/10/2018.
//

#ifndef SPECHAP_VCF_IO_H
#define SPECHAP_VCF_IO_H

#include "htslib/vcf.h"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "htslib/tbx.h"
#include "results.h"
#include "type.h"
#include "util.h"
extern bool IS_DEBUG;
extern int MIN_SNP_RECOM;
enum format_n {
    GT, PS
};

class VCFReader {
public:
    htsFile *vcf_file;
    bcf_hdr_t *header;
    int curr_contig, curr_bcf_contig;
    char *filename;
    int contigs_count;
    std::vector<std::string> contigs;   //rid to contig name appeared in header
    std::map<std::string, int> bcf_contig_map;
    std::map<uint, uint> var_count;      //rid to contig_variant_count appeared in vcf file
    hts_itr_t *iter;
    tbx_t *tbx_index;
    kstring_t tmp;
    std::vector<std::vector<std::string>> up_to_down;
    std::vector<std::vector<std::string>> down_to_up;
    std::unordered_map<std::string, int>* sample_to_idx;
    uint sample_count;

    bcf1_t *buffer;

public:
    explicit VCFReader(const char *filename);

    ~VCFReader();

    inline void reset() {
        bcf_close(vcf_file);
        bcf_hdr_destroy(header);
        vcf_file = bcf_open(filename, "r");
        header = bcf_hdr_read(vcf_file);
        curr_contig = 0;
        hts_itr_destroy(iter);
        tbx_destroy(tbx_index);
    }

    inline const void get_contigs() {
        const char **temp = tbx_seqnames(tbx_index, &contigs_count);
        for (int i = 0; i < contigs_count; i++) {
            contigs.emplace_back(std::string(temp[i]));
            temp[i] = nullptr;
        }
        free(temp);

        int bcf_contig_count;
        temp = bcf_hdr_seqnames(header, &bcf_contig_count);
        for (int i = 0; i < bcf_contig_count; i++) {
            bcf_contig_map[std::string(temp[i])] = i;
            temp[i] = nullptr;
        }
        free(temp);
    }

    inline int get_next_record(bcf1_t *record) {
        if (iter == nullptr)
            return -1;
        int status = tbx_itr_next(vcf_file, tbx_index, iter, &tmp);
        status = vcf_parse1(&tmp, header, record);
        if (status != 0) //eof or corruption in file
            return status;
        if (record->rid != this->curr_bcf_contig)
            return -1;
        bcf_unpack(record, BCF_UN_ALL);
        return status;
    }

    inline int get_next_record_contig(std::shared_ptr<VcfRecord> &result, bool use_gt) {
        //bcf_empty(record);
        //no such contig or iter not initialized
        if (iter == nullptr)
            return -1;
        //bcf1_t *r = bcf_init();
        int status = tbx_itr_next(vcf_file, tbx_index, iter, &tmp);
        if (status < 0)
            return status;

        status = vcf_parse1(&tmp, header, buffer);
        if (status != 0) //eof or corruption in file
            return -1;
        //new contig
        int ninfo_arr = 0, ninfo = 0;
        char * info = NULL;

        ninfo = bcf_get_info_string(header, buffer, "SVTYPE", &info, &ninfo_arr);

        if (ninfo < 0)
            result->bnd = false;
        else
            result->bnd = true;

        if (buffer->rid != this->curr_bcf_contig)
            return -1;
        //bcf_subset_format(header, record);
        bcf_unpack(buffer, BCF_UN_ALL);


        int i, j, ngt,nps,npl, nsmpl = bcf_hdr_nsamples(header);
        int *gt_arr = nullptr, *ps_arr = nullptr, ngt_arr =0, nps_arr = 0, *pl_arr = nullptr, npl_arr = 0;
        ngt = bcf_get_genotypes(header, buffer, &gt_arr, &ngt_arr);
        nps = bcf_get_format_int32(header, buffer, "PS", &ps_arr, &nps_arr);
        npl = bcf_get_format_int32(header, buffer, "PL", &pl_arr, &npl_arr);
        result->pos = buffer->pos;
        result->ID = buffer->rid;
        int max_ploidy = ngt / nsmpl;
        if (buffer->pos == 196110441) {
            auto mmmm = 9;
        }
//        if(nps <0 ) return -1;
        for (i = 0; i < nsmpl; i++) {
            int *ptr = gt_arr + i * max_ploidy;
            int *p_ptr = (ps_arr == nullptr? nullptr : ps_arr + i);
            int *pl_ptr = (pl_arr == nullptr? nullptr : pl_arr + i);
            auto call = new Call();
            call->allele1 = bcf_gt_allele(ptr[0]);
            call->allele2 = bcf_gt_allele(ptr[1]);
            if (pl_ptr != nullptr && *pl_ptr <= 0  && call->allele1 <0 && call->allele2 < 2) {
                call->allele1 = 0;
                call->allele2 = 0;
                call->un_sure = true;
            }
//            if(call->allele1 == -1) call->allele1 = 0;
//            if(call->allele2 == -1) call->allele2 = 0;
            p_ptr == nullptr || *p_ptr <=0 ? call->block_id = 0 : call->block_id = *p_ptr;
            call->pos = buffer->pos;
            result->calls.push_back(call);
        }
        free(gt_arr); free(ps_arr);
        //bcf_destroy(r);
        return 0;
    }

    inline int jump_to_contig(int contig_id) {
        curr_contig = contig_id;
        tbx_itr_destroy(iter);
        iter = nullptr;
        iter = tbx_itr_querys(tbx_index, contigs[contig_id].c_str());
        if (iter == nullptr) // no such contig
            return -1;
        curr_bcf_contig = bcf_contig_map[contigs[curr_contig]];
        return 0; //success
    }
};

class VCFWriter {
private:
    bcf_hdr_t *header;
    htsFile *fp;
    std::ofstream debugFile;
    int sample_count;
    std::map<format_n, const char *> FormatN;
    int ngt;
    int *gt;
    int *ps;
//    std::vector<long> dup_region;

private:
    void header_init();

public:
    VCFWriter(const bcf_hdr_t *hdr, const char *outfile, const char *reconFile);

    ~VCFWriter();
    std::vector<std::string> extractAltRef(char *als);
    void write_nxt_record(bcf1_t *record, std::shared_ptr<VcfRecord>, std::vector<int>& blk_no);
    void write_recom_duo(bcf1_t *record, const std::shared_ptr<VcfRecord>& result, int cidx, int pidx, int allele_idx, int *conflictFlag, int *conflictCount, std::string& recom_content, std::vector<std::string> & dup_region,std::vector<std::string> &simple_repeat, uint prev_count);
    void write_nxt_record_debug(bcf1_t *record, std::shared_ptr<VcfRecord> result, std::vector<int>& ps_nos, int *conflictFlag, std::vector<std::string> & dup_region);
    void write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf, uint prev_count);
    std::vector<std::string> parse_segdup(std::string contig);
    std::vector<std::string> parse_simple_repeat(std::string contig);

    std::string segdup_file;
    std::string simple_repeat_file;

    const std::string &getSimpleRepeatFile() const;

    void setSimpleRepeatFile(const std::string &simpleRepeatFile);
};


#endif //SPECHAP_VCF_IO_H
