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
    std::vector<int> mendel_passc;
    std::vector<int> mendel_passf;
    std::vector<int> mendel_passm;
    std::vector<std::shared_ptr<VcfRecord>> results_for_variant;
    std::set<int> conflicts1;
    std::set<int> conflicts2;
public:
    ChromoPhaser(const uint &chr_id, const std::string &chr_name, int nsmp);
    ~ChromoPhaser() = default;
    void correct_conflict(int idx);
    void add_result(const std::shared_ptr<VcfRecord>& result);
    void phase_with_hete(int idx1, int idx2, int side, InfoSet* infoSet);
    void phase_with_homo(int idx1, int idx2, int side, InfoSet* infoSet);
    void phase_with_homo2(int idx1, int idx2, int side, InfoSet* infoSet);
    void check_mendel(int idx1, int idx2, int idx3);
    void extend(int idx, InfoSet* infoSet, int side);
    inline bool is_x() const {
        return chr_name.find('x') != std::string::npos || chr_name.find('X') != std::string::npos;
    }

    inline bool is_y() const {
        return chr_name.find('y') != std::string::npos || chr_name.find('Y') != std::string::npos;
    }
//    inline uint get_var_pos(uint idx) { return results_for_variant[idx]->get_pos(); }
};

#endif

