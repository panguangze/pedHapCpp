#ifndef SPECHAP_UTIL_H
#define SPECHAP_UTIL_H

#include <ios>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>


void logging(std::ostream &stream, const std::string & message);
std::vector<std::string> split(const std::string &s, char delim);
bool check_contains(std::vector<std::string> regions, long pos, long end_pos);
//bool check_mendel(au)
#endif