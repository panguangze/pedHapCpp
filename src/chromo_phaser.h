//
// Created by yyh on 3/4/2019.
//

#ifndef SPECHAP_PHASER_H
#define SPECHAP_PHASER_H

#include "type.h"
#include "vcf_io.h"
#include "util.h"

extern bool ONLY_CHILD;
class Phaser
{
public:
    std::unordered_map<std::string, int>* sample_to_idx;
    uint sample_count;
    ChromoPhaser* chromoPhaser{};
    std::vector<std::vector<std::string>> up_to_down;
    std::vector<std::vector<std::string>> down_to_up;
    explicit Phaser(const std::string & fnvcf, const std::string & fnout);
    ~Phaser();
    void phasing();

private:
    VCFReader *frvcf;
    VCFWriter *fwvcf;
    void phasing_by_chrom() const;
    int load_contig_blocks(ChromoPhaser *chromo_phaser);
};

#endif //SPECHAP_PHASER_H
