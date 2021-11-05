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

class ChromoPhaser
{
public:
    //initialized parameter
    uint chr_id;
    std::string chr_name;
    uint variant_count;         //variant no to phase
    uint sample_count;
    std::vector<int> mendel_cs;
    std::vector<int> mendel_pass;
    std::vector<std::shared_ptr<VcfRecord>> results_for_variant;
public:
    ChromoPhaser(const uint &chr_id, const std::string &chr_name, int nsmp);
    ~ChromoPhaser() = default;
    void add_result(const std::shared_ptr<VcfRecord>& result);
    void phase_with_hete(int idx1, int idx2, int side);
    void phase_with_homo(int idx1, int idx2, int side);
    bool check_mendel(int idx1, int idx2, int idx3);
    void extend(int idx, InfoSet& infoSet, int side);
//    inline uint get_var_pos(uint idx) { return results_for_variant[idx]->get_pos(); }
};

#endif

