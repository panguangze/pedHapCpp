//
// Created by yyh on 3/4/2019.
//

#ifndef SPECHAP_PHASER_H
#define SPECHAP_PHASER_H

#include "type.h"
#include "vcf_io.h"
#include "util.h"

extern bool HYBRID;
extern bool KEEP_PS;
extern int MAX_HIC_INSERTION;
extern int RECURSIVE_LIMIT;

class Phaser
{
public:
    std::unordered_map<std::string, int>* sample_to_idx;
    uint sample_count;
    explicit Phaser(const std::string & fnvcf, const std::string & fnout);
    ~Phaser();
    void phasing();

private:
    void sort_frag_file(std::string file_name);
    double threshold;
    VCFReader *frvcf;
    VCFWriter *fwvcf;
    int coverage;
    void phasing_by_chrom(uint var_count, ChromoPhaser *chromo_phaser, std::vector<std::string>&);
    int load_contig_records(ChromoPhaser *chromo_phaser);
    int load_contig_blocks(ChromoPhaser *chromo_phaser);
};

#endif //SPECHAP_PHASER_H
