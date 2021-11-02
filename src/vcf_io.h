//
// Created by yonghanyu2 on 10/10/2018.
//

#ifndef SPECHAP_VCF_IO_H
#define SPECHAP_VCF_IO_H

#include "htslib/vcf.h"
#include <map>
#include <vector>
#include <string>
#include "htslib/tbx.h"
#include "results.h"
#include "type.h"

typedef std::map<filter_type, std::string> filter_map_t;

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

    /*
    inline int get_AF0(float *af0, bcf1_t *record)
    {
        float **AF = nullptr;
        int32_t l_AF;
        bcf_get_info_float(header, record, "AF", &AF, &l_AF);
    }*/

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
        if (buffer->rid != this->curr_bcf_contig)
            return -1;
        //bcf_subset_format(header, record);
        bcf_unpack(buffer, BCF_UN_ALL);


        int i, j, ngt,nps, nsmpl = bcf_hdr_nsamples(header);
        int *gt_arr = nullptr, *ps_arr = nullptr, ngt_arr =0, nps_arr = 0;
        ngt = bcf_get_genotypes(header, buffer, &gt_arr, &ngt_arr);
        nps = bcf_get_format_int32(header, buffer, "PS", &ps_arr, &nps_arr);

        int max_ploidy = ngt / nsmpl;
        for (i = 0; i < nsmpl; i++) {
            int *ptr = gt_arr + i * max_ploidy;
            int *p_ptr = ps_arr + i;
            auto call = new Call();
            call->allele1 = bcf_gt_allele(ptr[0]);
            call->allele2 = bcf_gt_allele(ptr[1]);
            p_ptr == nullptr ? call->ps = 0 : call->ps = *p_ptr;
            call->pos = buffer->pos;
            call->phased = bcf_gt_is_phased(ptr[0]);
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
    int sample_count;
    std::map<filter_type, const char *> FilterN;
    std::map<format_n, const char *> FormatN;
    filter_map_t filter_map =
            {
                    {filter_type::NOINFO,                              "INFO_NOT_ENOUGH"},
                    {filter_type::POOLRESULT,                          "POOL_SPECTRAL_RESULT"},
                    {filter_type::TENXINCONSISTENCY,                   "10X_PHASING_INCONSISTENCY"},
                    {filter_type::PASS,                                "PASS"},
                    {filter_type::CONFILCTINGRESULT,                   "WINDOW_RESULT_INCONSISTENCY"},
                    {filter_type::TENX_ALLELE_FRENCENCY_FILTER,        "10X_ALLELE_FREQUENCY_FILTER"},
                    {filter_type::LOW_COVERAGE,                        "LOW_COVERAGE"},
                    {filter_type::TENX_QUAL_FILTER,                    "TENX_QUAL_FILTER"},
                    {filter_type::TENX_RESCUED_MOLECUE_HIGH_DIVERSITY, "TENX_RESCUED_MOLECUE_HIGH_DIVERSITY"}
            };
    int ngt;
    int *gt;

private:
    void header_init();

public:
    VCFWriter(const bcf_hdr_t *hdr, const char *outfile);

    ~VCFWriter();

    void write_nxt_record(bcf1_t *record, std::shared_ptr<VcfRecord>, std::vector<int>& blk_no);

    void write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf);
};


#endif //SPECHAP_VCF_IO_H
