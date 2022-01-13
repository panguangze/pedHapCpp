//
// Created by yyh on 1/7/2019.
//

#ifndef SPECHAP_RESULTS_H
#define SPECHAP_RESULTS_H
#include <bitset>
#include <memory>
#include <unordered_map>
#include "vector"

class Call : public std::enable_shared_from_this<Call>{
public:
    int allele1;
    int allele2;
    int block_id;
    int pos;
    bool need_flip;
    Call();
    ~Call() = default;
    inline bool isHomo() {
        return allele1 == allele2;
    }
    inline bool isPhased() {
        return block_id > 0;
    }
    inline void flip() {
        auto t = allele1;
        allele1 = allele2;
        allele2 = t;
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
    bool bnd;

    VcfRecord();
    ~VcfRecord() = default;
};


#endif //SPECHAP_RESULTS_H