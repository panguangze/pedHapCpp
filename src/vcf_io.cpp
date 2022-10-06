//
// Created by yonghanyu2 on 10/10/2018.
//

#include "vcf_io.h"
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#define EMPTY_ID "null"


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
    delete sample_to_idx;
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

VCFWriter::VCFWriter(const bcf_hdr_t *hdr, const char *outfile, const char *reconFile) : debugFile(reconFile,std::ofstream::app)
{
//    std::ofstream df("debug.txt");
//    debugFile("debug.txt");
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
//        if (tcall->un_sure) {
//            gt[2*i] = bcf_gt_unphased(-1);
//            gt[2*i + 1] = bcf_gt_unphased(-1);
//        }
        if (tcall->isPhased() || tcall->isHomo()) {
            gt[2*i] = bcf_gt_phased(tcall->allele1);
            gt[2*i + 1] = bcf_gt_phased(tcall->allele2);
        } else {
            gt[2*i] = bcf_gt_unphased(tcall->allele1);
            gt[2*i + 1] = bcf_gt_unphased(tcall->allele2);
        }
        if(result->calls[i]->isHomo()) {
            ps[i] = 0;
//            ps[i] = ps_nos[i];
        } else{
            ps[i] = ps_nos[i];
        }
    }
    bcf_update_genotypes(this->header, record, gt, ngt);
    bcf_update_format_int32(this->header, record, "PS", ps, sample_count);
    bcf1_t *record_w = bcf_dup(record);
    int rr = bcf_write(this->fp, this->header, record_w);
    bcf_destroy(record_w);
}
std::vector<std::string> VCFWriter::extractAltRef(char *als) {
    std::vector<std::string> res;
    std::string str = "";
    int i = 0;
    int j = 0;
    while(j != 2) {
        if(als[i] == '\0') {
            res.push_back(str);
            str.clear();
            j++;
        } else {
            str += als[i];
        }
        i++;
    }
    return res;
}

void VCFWriter::write_recom_duo(bcf1_t *record, const std::shared_ptr<VcfRecord>& result, int cidx, int pidx, int allele_idx,
                                int *conflictFlag, int *conflictCount, std::string &recom_content, std::vector<std::string> &dup_region, std::vector<std::string> &simple_repeat, uint prev_count) {
//    int i, j,nsmpl = bcf_hdr_nsamples(header);
    auto ccall = result->calls[cidx];
    auto pcall = result->calls[pidx];
    if(record->pos == 196110562) {
        int tmp = 0;
    }
    if (!ccall->isHomo() && ccall->isPhased() && ccall->block_id != 1) {
        *conflictFlag += 1;
        return;
    }
    if (ccall->un_sure || pcall->un_sure) {
//        if (*conflictFlag == 1) {
//            *conflictFlag = 2;
//        }
        return;
    }
    if (!ccall->isPhased() && !ccall->isHomo()) {
//        *conflictFlag = 3;
        return;
    }
//    if ((!ccall->isPhased() && !ccall->isHomo()) || (!ccall->isHomo() && ccall->isPhased() && ccall->block_id != 1)) {
//        *conflictFlag = 0;
//        return;
//    }
    if (pcall->isHomo() || pcall->block_id != 1) {
        return;
    }
//    auto pcall = result->calls[pidx];
    auto c1 = ccall->allele1;
    auto c2 = ccall->allele2;
    auto p1 = pcall->allele1;
    auto p2 = pcall->allele2;
    if (allele_idx == 1) {
        c1 = ccall->allele2;
        c2 = ccall->allele1;
    }
//    if (ccall->un_sure || pcall->un_sure) {
//        *conflictFlag = 2;
//        return;
//    }
    if (c1 == -1 || c2 == -1 || p1 == -1 || p2 == -1) {
        return;
    }
    std::string line;
//    line.append("\t");
//    line.append(std::to_string(prev_count));
    line.append(std::to_string(result->pos + 1));
//    if (*conflictFlag == 2) {
//        line.append("_un_sure_GT");
//    }
//    if (*conflictFlag == 3) {
//        line.append("_un_sure_GT");
//    }
    if (c1 != p1 && c1 == p2) {
        if (*conflictFlag >=2 ) {
            recom_content.append("\n");
        }
        *conflictFlag = 0;
        line.append("_");
        if (allele_idx == 0) {
            line.append(std::to_string(c1));
            line.append("|");
            line.append(std::to_string(c2));
        } else {
            line.append(std::to_string(c2));
            line.append("|");
            line.append(std::to_string(c1));
        }
//        line.append(std::to_string(c1));
//        line.append("|");
//        line.append(std::to_string(c2));
        line.append("_");
        line.append(std::to_string(p1));
        line.append("|");
        line.append(std::to_string(p2));
//        std::string line = std::to_string(record->pos);
//        line.append("\t");
        line.append("_");
        line.append("conflict");
        int *ends = nullptr;
        int end_n = 0;
        if (result->bnd) {
            bcf_get_info_int32(header, record, "END", &ends, &end_n);
            line.append("_BND");
            if (*ends - record->pos >= 100000) {
                line.append("LONG");
            }
        }
//        if (result->bnd && check_contains(dup_region, record->pos, *ends) && (c1 != 0 or c2 != 0)) {
//            if (*ends - record->pos >= 100000) {
//                line.append("\ttoo long");
//            }
//        }
        recom_content.append(line + "\t");
//        conflictCount++;
        free(ends);
    } else {
//        if (*conflictFlag < 2) {
//            *conflictFlag = 2;
//        } else {
//            *conflictFlag = 0;
//        }
        *conflictFlag += 1;

//        if (*conflictFlag == 1) {
//            recom_content.append("\n");
//        }
//        if ((c1 != c2)) {
//            line.append("\t");
//            *conflictFlag = 0;
//        }
//        std::vector<std::string> altRef = extractAltRef(record->d.als);
//        line.append("\t");
//        line.append(altRef[1]);
    }


//    if(ccall->isPhased() || ccall->isHomo()) {
//        line.append("|");
//    } else {
//        line.append("/");
//    }
//    line.append(std::to_string(ccall->allele2));
//    line.append(":");
//    for (i = 0; i < nsmpl; i++) {
//        line.append("\t");
//        auto tcall = result->calls[i];
//        line.append(std::to_string(tcall->allele1));
//        if(tcall->isPhased() || tcall->isHomo()) {
//            line.append("|");
//        } else {
//            line.append("/");
//        }
//        line.append(std::to_string(tcall->allele2));
//        line.append(":");
//        if (tcall->isHomo()) {
//            line.append("1");
//        } else {
//            line.append(std::to_string(ps_nos[i]));
//        }
//    }
//    auto child_ps = std::to_string(ps_nos[0]);


//    if (result->bnd) {

//        auto m1 = result->calls[2]->allele1;
//        auto m2 = result->calls[2]->allele2;
//        if (((c1 == f1 && c2 == m1) || (c1 == f1 && c2 == m2) || (c1 == f2 && c2 == m1) || (c1 == f2 && c2 == m2) ||
//             (c2 == f1 && c1 == m1) || (c2 == f1 && c1 == m2) || (c2 == f2 && c1 == m1) || (c2 == f2 && c1 == m2))){
//
//        } else {
//            if (!result->bnd) {
//                return;
//            }
//        }
//        if(!result->bnd && !(((c1 == f1 && c2 == m1) || (c1 == f1 || c2 == m2)) || ((c2 == f1 || c2 == f2) && (c1 == m1 || c1 == m2)))) {
//            return;
//        }
//                    if (c1 != f1 || c2 !=m1)
}

void VCFWriter::write_nxt_record_debug(bcf1_t *record, std::shared_ptr<VcfRecord> result, std::vector<int>& ps_nos, int *conflictFlag, std::vector<std::string> & dup_region)
{
    int i, j,nsmpl = bcf_hdr_nsamples(header);
    std::string line = std::to_string(record->pos);
    for (i = 0; i < nsmpl; i++) {
        line.append("\t");
        auto tcall = result->calls[i];
        line.append(std::to_string(tcall->allele1));
        if(tcall->isPhased() || tcall->isHomo()) {
            line.append("|");
        } else {
            line.append("/");
        }
        line.append(std::to_string(tcall->allele2));
        line.append(":");
        if (tcall->isHomo()) {
            line.append("1");
        } else {
            line.append(std::to_string(ps_nos[i]));
        }
    }
    auto child_ps = std::to_string(ps_nos[0]);

    if (true) {

//    if (result->bnd) {
        auto c1 = result->calls[0]->allele1;
        auto c2 = result->calls[0]->allele2;
        auto f1 = result->calls[1]->allele1;
        auto f2 = result->calls[1]->allele2;
        auto m1 = result->calls[2]->allele1;
        auto m2 = result->calls[2]->allele2;
//        if (((c1 == f1 && c2 == m1) || (c1 == f1 && c2 == m2) || (c1 == f2 && c2 == m1) || (c1 == f2 && c2 == m2) ||
//             (c2 == f1 && c1 == m1) || (c2 == f1 && c1 == m2) || (c2 == f2 && c1 == m1) || (c2 == f2 && c1 == m2))){
//
//        } else {
//            if (!result->bnd) {
//                return;
//            }
//        }
//        if(!result->bnd && !(((c1 == f1 && c2 == m1) || (c1 == f1 || c2 == m2)) || ((c2 == f1 || c2 == f2) && (c1 == m1 || c1 == m2)))) {
//            return;
//        }
//                    if (c1 != f1 || c2 !=m1)
        if (c1 != f1) {
            line.append("\t");
            line.append("conflict");
            if (*conflictFlag == 0) {
                this->debugFile << "xxxxxxxxxxxxxxxxx\n";
            }
            *conflictFlag = 1;
        } else {
            if ((c1 != c2)) {
                line.append("\t");
                *conflictFlag = 0;
            }
            std::vector<std::string> altRef = extractAltRef(record->d.als);
            line.append("\t");
            line.append(altRef[1]);
            if ((c1 == c2) || !result->calls[0]->isPhased() || !result->calls[1]->isPhased()) {

            } else if (c1 != f1)
                this->debugFile << line << "\n";
            int *ends = nullptr;
            int end_n;
            if (result->bnd) {
                bcf_get_info_int32(header, record, "END", &ends, &end_n);
            }
            if (result->bnd  && f1 == 0 && f2 == 0 && m1 == 0 &&
                m2 == 0 && (c1 != 0 or c2 != 0)) {
                if (check_contains(dup_region, record->pos, *ends))
                    line.append("_NAHR");
                int zero = 0;
                int *start = &zero ,*end = &zero;
                if (check_contains_repeat(dup_region,record->pos,start,end)) {
                    line.append("SIMPLER_");
                    line.append(std::to_string(*start));
                    line.append("_");
                    line.append(std::to_string(*end));line.append("_");
                }
                if (*ends - record->pos >= 100000) {
                    line.append("\ttoo long");
                }
                this->debugFile << line << "\n";
            }
            if (result->bnd) free(ends);
        }
    }
//    this->debugFile<<line<<"\n";
}



void VCFWriter::write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf, uint prev_count)
{
    bcf1_t *record = bcf_init();
    frvcf.jump_to_contig(frvcf.curr_contig);
    std::vector<std::unordered_map<uint, uint>> idx2pos;
    for (int i = 0; i<sample_count;i++) {
        std::unordered_map<uint, uint> tmp;
        idx2pos.push_back(tmp);
    }
    int *conflictFlag, v = 0;
    conflictFlag = &v;
    auto dup_region = parse_segdup(frvcf.contigs[frvcf.curr_contig]);
    auto simple_repeat_region = parse_simple_repeat(frvcf.contigs[frvcf.curr_contig]);
//    this->debugFile<<"chrom\t"<<contig;
    std::unordered_map<std::string, std::string> sample_name2recom;
    int* fConflictFlags = new int [frvcf.up_to_down.size()];
    int* mConflictFlags = new int [frvcf.up_to_down.size()];
    int* fConflictCounts = new int [frvcf.up_to_down.size()];
    int* mConflictCounts = new int [frvcf.up_to_down.size()];
//    int* prev2Snp = new int [frvcf.up_to_down.size()];
    int trio_idx = 0;
    int s_idx=-1, f_idx=-1, m_idx=-1;

    for(auto it: (frvcf.up_to_down)) {
        s_idx = (*frvcf.sample_to_idx)[it[0]];
        sample_name2recom[it[0]] = "";
        if (it[1] != EMPTY_ID) {
            f_idx = (*frvcf.sample_to_idx)[it[1]];
            sample_name2recom[it[1]] = "";
        }
        if (it[2] != EMPTY_ID) {
            m_idx = (*frvcf.sample_to_idx)[it[2]];
//            sample_name2recom[it[2]];
            sample_name2recom[it[2]] = "";
        }
    }

    for (uint idx = 0; idx < chromo_phaser->variant_count; idx++)
    {
        prev_count = prev_count + 1;
        trio_idx = 0;
        frvcf.get_next_record(record);
        auto result = chromo_phaser->results_for_variant[idx];
        if (result->pos == 90234140) {
            int tmp = 99;
        }
        bcf_translate(this->header, frvcf.header, record);
        std::vector<int> ps_nos;
        for (auto tcall : result->calls) {
            ps_nos.push_back(tcall->block_id);
        }

        if(IS_DEBUG) {
            for(auto it: (frvcf.up_to_down)) {
                s_idx = (*frvcf.sample_to_idx)[it[0]];
                if(it[1] != EMPTY_ID && it[2] != EMPTY_ID) {
                    std::string line;
                    if (result->bnd && result->calls[1]->allele1 <= 0 && result->calls[1]->allele2 <= 0
                        && result->calls[2]->allele1 <= 0 && result->calls[2]->allele2 <= 0) {
                        if (result->calls[0]->isPhased() && result->calls[0]->allele1 == 1) {
                            int *ends = nullptr;
                            int end_n = 0;
                            bcf_get_info_int32(header, record, "END", &ends, &end_n);
                            line.append("\n");
                            line.append(std::to_string(result->pos+1));
                            line.append("_");
                            line.append(std::to_string(result->calls[0]->allele1));
                            line.append("|");
                            line.append(std::to_string(result->calls[0]->allele2));
                            line.append("_");
                            line.append(std::to_string(result->calls[1]->allele1));
                            line.append("|");
                            line.append(std::to_string(result->calls[1]->allele2));
                            line.append("_conflict_");
                            if (check_contains(dup_region, record->pos, *ends)) {
                                line.append("NAHR");
                            } else {
                                int zero = 0;
                                int *start = &zero ,*end = &zero;
                                if (check_contains_repeat(simple_repeat_region,record->pos,start,end)) {
                                    line.append("SIMPLER_");
                                    line.append(std::to_string(*start));
                                    line.append("_");
                                    line.append(std::to_string(*end));line.append("_");
                                }
                                line.append("UNNAHR");
                            }
                            sample_name2recom[it[1]].append(line);
                            free(ends);
                        } else if (result->calls[0]->isPhased() && result->calls[0]->allele2 == 1) {
                            int *ends = nullptr;
                            int end_n = 0;
                            bcf_get_info_int32(header, record, "END", &ends, &end_n);
                            line.append("\n");
                            line.append(std::to_string(result->pos+1));
                            line.append("_");
                            line.append(std::to_string(result->calls[0]->allele1));
                            line.append("|");
                            line.append(std::to_string(result->calls[0]->allele2));
                            line.append("_");
                            line.append(std::to_string(result->calls[2]->allele1));
                            line.append("|");
                            line.append(std::to_string(result->calls[2]->allele2));
                            line.append("_conflict_");
                            if (check_contains(dup_region, record->pos, *ends)) {
                                line.append("NAHR");
                            } else {
                                int zero = 0;
                                int *start = &zero ,*end = &zero;
                                if (check_contains_repeat(simple_repeat_region,record->pos,start,end)) {
                                    line.append("SIMPLER_");
                                    line.append(std::to_string(*start));
                                    line.append("_");
                                    line.append(std::to_string(*end));line.append("_");
                                }
                                line.append("UNNAHR");
                            }
                            sample_name2recom[it[2]].append(line);
                            free(ends);
                        } else if(!result->calls[0]->isPhased() && !result->calls[0]->isHomo()) {
                            int *ends = nullptr;
                            int end_n = 0;
                            bcf_get_info_int32(header, record, "END", &ends, &end_n);
                            line.append("\n");
                            line.append(std::to_string(result->pos+1));
                            line.append("_");
                            line.append(std::to_string(result->calls[0]->allele1));
                            line.append("\\");
                            line.append(std::to_string(result->calls[0]->allele2));
                            line.append("_");
                            line.append(std::to_string(result->calls[2]->allele1));
                            line.append("|");
                            line.append(std::to_string(result->calls[2]->allele2));
                            line.append("_conflict_");
                            if (check_contains(dup_region, record->pos, *ends)) {
                                line.append("NAHR\n");
                            } else {
                                int zero = 0;
                                int *start = &zero ,*end = &zero;
                                if (check_contains_repeat(simple_repeat_region,record->pos,start,end)) {
                                    line.append("SIMPLER_");
                                    line.append(std::to_string(*start));
                                    line.append("_");
                                    line.append(std::to_string(*end));line.append("_");
                                }
                                line.append("UNNAHR\n");
                            }
                            sample_name2recom[it[0]].append(line);
                            free(ends);
                        }
                        continue;
                    }

//                    if ((result->calls[0]->allele1 == 1 || result->calls[0]->allele2 == 1)
//                        && result->calls[1]->allele1 <= 0 && result->calls[1]->allele2 <= 0
//                        && result->calls[2]->allele1 <= 0 && result->calls[2]->allele2 <= 0) {
//                        sample_name2recom[it[0]].append()
//                    }
                }
                if (it[1] != EMPTY_ID) {
                    f_idx = (*frvcf.sample_to_idx)[it[1]];
                    write_recom_duo(record,result,s_idx,f_idx,0,fConflictFlags+trio_idx, fConflictCounts+trio_idx, sample_name2recom[it[1]],dup_region, simple_repeat_region, prev_count);
                }
                if (it[2] != EMPTY_ID) {
                    m_idx = (*frvcf.sample_to_idx)[it[2]];
                    write_recom_duo(record,result,s_idx,m_idx,1,mConflictFlags+trio_idx, mConflictCounts+trio_idx, sample_name2recom[it[2]],dup_region,simple_repeat_region, prev_count);
                }
//                write_nxt_record_debug(record, result, ps_nos, conflictFlag, dup_region);
                trio_idx++;
            }
            write_nxt_record(record, result, ps_nos);
        } else {
            write_nxt_record(record, result, ps_nos);
        }
    }
    int i = 0;
    for (const auto& it : sample_name2recom) {
        if (i == 0)
            this->debugFile<<"chrom\t"<<contig;
        else
            this->debugFile<<"\nchrom\t"<<contig;
        i++;
        this->debugFile<<"\t"<<it.first<<it.second<<"\n";
    }
    free(fConflictFlags);
    free(mConflictFlags);
    free(mConflictCounts);
    free(fConflictCounts);
    bcf_destroy(record);
}

std::vector<std::string> VCFWriter::parse_segdup(std::string contig) {
    std::vector<std::string> result;
    std::ifstream infile(segdup_file);
    int start,end;
    std::string line;
    while (std::getline(infile, line)) {
        auto vs = split(line, '\t');
        if (vs[0] != contig || vs[3] != contig) continue;
        result.push_back(line);
//        result.push_back(std::stoi(vs[2]));
    }
    infile.close();
    return result;
}

std::vector<std::string> VCFWriter::parse_simple_repeat(std::string contig) {
    std::vector<std::string> result;
    std::ifstream infile(segdup_file);
    int start,end;
    std::string line;
    while (std::getline(infile, line)) {
        auto vs = split(line, '\t');
        if (vs[0] != contig || vs[3] != contig) continue;
        result.push_back(line);
//        result.push_back(std::stoi(vs[2]));
    }
    infile.close();
    return result;
}

const std::string &VCFWriter::getSimpleRepeatFile() const {
    return simple_repeat_file;
}

void VCFWriter::setSimpleRepeatFile(const std::string &simpleRepeatFile) {
    simple_repeat_file = simpleRepeatFile;
}

