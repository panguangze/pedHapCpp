//
// Created by yyh on 1/7/2019.
//

#ifndef SPECHAP_RESULTS_H
#define SPECHAP_RESULTS_H
#include <bitset>
#include <memory>
#include <unordered_map>
#include "vector"

enum class filter_type
        {
            PASS,
            NOINFO,
            POOLRESULT,
            LOW_COVERAGE,
            TENXINCONSISTENCY,
            CONFILCTINGRESULT,
            TENX_ALLELE_FRENCENCY_FILTER,
            TENX_QUAL_FILTER,
            TENX_RESCUED_MOLECUE_HIGH_DIVERSITY
        };



typedef std::bitset<1> shap;
typedef unsigned int uint;

template <typename T>
bool is_uninitialized(std::weak_ptr<T> const& weak)
{
    using wt = std::weak_ptr<T>;
    return !weak.owner_before(wt{}) && !wt{}.owner_before(weak);
}


class PhasedBlock;

class Call : public std::enable_shared_from_this<Call>{
public:
    bool phased;
    int allele1;
    int allele2;
    int ps;
    int pos;
    std::weak_ptr<PhasedBlock> block;
    bool need_flip;
    Call();
    ~Call() = default;
    inline void flip() {
        need_flip = !this->need_flip;
    }
    inline bool isHomo() {
        return allele1 == allele2;
    }
    inline bool isPhased() {
        return ps!=0;
    }
};


class VcfRecord{
public:
    std::unordered_map<std::string,uint> sample_to_index;
    std::vector<Call* > calls;
    uint pos;
    std::string ID;
    std::string ref;
    std::string alts;

    VcfRecord();
    ~VcfRecord() = default;
};

//class CoverInfo{
//public:
//    uint ps;
//    std::unordered_map<uint, >
//};

#endif //SPECHAP_RESULTS_H