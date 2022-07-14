//
// Created by yyh on 3/4/2019.
//

#ifndef SPECHAP_PHASER_H
#define SPECHAP_PHASER_H

#include <utility>

#include "type.h"
#include "vcf_io.h"
#include "util.h"
#define EMPTY_ID "null"


extern bool ONLY_CHILD;
extern bool XY;
extern bool IS_MALE;
class Phaser
{
public:
//    std::unordered_map<std::string, int>* sample_to_idx;
    uint sample_count;
    ChromoPhaser* chromoPhaser{};
//    std::vector<std::vector<std::string>> up_to_down;
//    std::vector<std::vector<std::string>> down_to_up;
    std::vector<std::string> contigs;
    explicit Phaser(const std::string & fnvcf, const std::string & fnout, const std::string & frecon);
    ~Phaser();
    void phasing();

    VCFReader *frvcf;
    VCFWriter *fwvcf;
    void phasing_by_chrom() const;
    void correct() const;
    int load_contig_blocks(ChromoPhaser *chromo_phaser);

    void phasing_xy() const;
    void set_contigs(std::vector<std::string> contig_vec) {
        this->contigs = std::move(contig_vec);
    }
    void set_segdup(const std::string & dupfile) {
        this->fwvcf->segdup_file = dupfile;
    }
    void set_trio(std::vector<std::string> trio, bool is_up_to_down) {
        if (is_up_to_down) this->frvcf->up_to_down.push_back(trio);
        else this->frvcf->down_to_up.push_back(trio);
    }
    std::vector<std::vector<std::string>> get_up_to_down() const;
    std::vector<std::vector<std::string>> get_down_to_up() const;
    std::unordered_map<std::string, int>* get_sample_to_idx() const;
};

#endif //SPECHAP_PHASER_H
