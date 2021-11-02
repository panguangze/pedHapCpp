
// Note that the index maintains different meaning in different scope
// Variant idx: idx of variants on chromosome according to its order in VCF file, 0 based
// Variant idx in fragment: order according to VCF file, 1 based
// Block idx: always start from 0
// Block Matrix idx: same as above

#ifndef SPECHAP_TYPE_H
#define SPECHAP_TYPE_H

#include "results.h"
#include "phaseInfo.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <cmath>
#include <map>
typedef std::pair<int, double> allele_info;
typedef std::pair<uint, allele_info> snp_info;
typedef std::vector<snp_info> snp_container;

#define FRAG_HIC    0
#define FRAG_10X    1
#define FRAG_NORMAL 2

#define MODE_10X        0
#define MODE_HIC        1
#define MODE_PE         2
#define MODE_PACBIO     3
#define MODE_NANOPORE   4
#define MODE_HYBRID     5
class Fragment
{
public:
    Fragment() = default;
    Fragment(const Fragment &rhs);
    ~Fragment() = default;
    uint start, end;
    int type = FRAG_NORMAL;               
    snp_container snps;         //<variant_idx, allele>
    std::string barcode;
    double read_qual = 0;
    int insertion_size = 0;
    bool rescued = false;
    double dm = 0;
    inline void reset() { this->snps.clear(); }
    inline void insert(snp_info &snp) { snps.push_back(snp); } //potential bug here
    inline void update_start_end()
    {
        start = this->snps[0].first;
        end = this->snps.back().first;
    }
};


class PhasedBlock : public std::enable_shared_from_this<PhasedBlock>
{
public:
    uint start_variant_idx;                                             //variant idx
    uint end_variant_idx;
    std::unordered_map<uint, std::shared_ptr<Call>> calls;           //key variant idx
    std::set<uint> variant_idxes;
    uint block_id = 0;                                                        //for output

public:///=//public:

    PhasedBlock();
    PhasedBlock(uint start_idx);
    PhasedBlock(uint start_idx, std::shared_ptr<VcfRecord> result);
    PhasedBlock(const PhasedBlock &rhs);
    PhasedBlock(const PhasedBlock &rhs, std::unordered_map<uint, std::shared_ptr<VcfRecord>> &results_dup);
    ~PhasedBlock();

    //getter
    inline std::size_t size() { return  variant_idxes.size(); }
    inline uint get_start_position() { return start_variant_idx; }
    inline void add_call(Call *call) {
        this->calls.emplace(call->ps, call);
        this->end_variant_idx = call->pos;
    }
    //setter
    void flip(int idx);                                                          //flip the haplotype of the block


    //utility
//    inline bool contain_variant(uint idx) { return calls.count(idx) != 0; }
//    inline void insert_variant(uint variant_idx, std::shared_ptr<VcfRecord> result)
//    {
//        calls[variant_idx] = result;
//        variant_idxes.insert(variant_idx);
//    }
};


typedef std::shared_ptr<PhasedBlock> ptr_PhasedBlock;
template class std::unordered_map<uint, ptr_PhasedBlock>;
typedef std::unordered_map<uint ,ptr_PhasedBlock> sample_PhasedBlocks;
class ChromoPhaser
{
public:
    //initialized parameter
    uint chr_id;
    std::string chr_name;
    uint variant_count;         //variant no to phase
    uint init_block_count;
    uint sample_count;

    //to be calculated afterwards

    //ChromoPhaser is obliged to manage the following two container
    std::vector<std::shared_ptr<VcfRecord>> results_for_variant;
    std::vector<std::unordered_map<uint, uint>*> variant_to_block_id;
    std::vector<sample_PhasedBlocks* > phased_blocks_info;

public:
    ChromoPhaser(const uint &chr_id, const std::string &chr_name, int nsmp);
    ~ChromoPhaser();
    void setSampleCount(int count);
    void add_result(std::shared_ptr<VcfRecord> result);
    void phase_with_hete(int idx1, int idx2, int side);
    void phase_with_homo(int idx1, int idx2, int side);
    void extend(int idx, InfoSet& infoSet);
//    inline uint get_var_pos(uint idx) { return results_for_variant[idx]->get_pos(); }


};

#endif

