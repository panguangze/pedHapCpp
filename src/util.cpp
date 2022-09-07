#include "util.h"
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>


void logging(std::ostream &stream, const std::string & message)
{
    std::time_t t = std::time(nullptr);

    stream << "[pedHap " << std::put_time(std::localtime(&t), "%Y:%m:%d %H:%M:%S]") << message << std::endl;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
        // elems.push_back(std::move(item)); // if C++11 (based on comment from @mchiasson)
    }
    return elems;
}

bool check_contains(std::vector<std::string> regions, long pos, long end_pos) {
    for (int i = 0; i < regions.size(); i++) {
        auto vs = split(regions[i], '\t');
        auto s1 = std::stoi(vs[1]);
        auto e1 = std::stoi(vs[2]);
        auto s2 = std::stoi(vs[4]);
        auto e2 = std::stoi(vs[5]);
        if (pos >= s1 && pos <= e1 && end_pos >= s2 && end_pos <= e2){
            auto t = regions[i];
            auto t1 = regions[i+1];
            return true;
        }
    }
    return false;
}

bool check_contains_repeat(std::vector<std::string> regions, long pos, int *start, int *end) {
    for (int i = 0; i < regions.size(); i++) {
        auto vs = split(regions[i], '\t');
        auto s1 = std::stoi(vs[1]);
        auto e1 = std::stoi(vs[2]);
        if (pos >= s1 && pos <= e1){
            auto t = regions[i];
            auto t1 = regions[i+1];
            *start = s1;
            *end = e1;
            return true;
        }
    }
    return false;
}