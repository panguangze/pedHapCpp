#include <iostream>
#include <cstring>
#include <fstream>
#include "src/cxxopts.hpp"
#include "vcflib/Variant.h"
using namespace std;

//TODO this is a tmp
int main(int argc, char *argv[]) {
    cxxopts::Options options("localhap", "Local Haplotype constructer");
    options.add_options()
            ("vcf", "operate: check or solve", cxxopts::value<std::string>())
            ("ped", "If true, across chromo integration will take into consideration", cxxopts::value<bool>()->default_value("false"))
            ("out", "Junction database", cxxopts::value<std::string>())
            ("t1", "Input lh file", cxxopts::value<int>())
            ("t2", "Checked local hap input file, lh format", cxxopts::value<int>())
            ("bnd", "ILP out file prefix, only for check", cxxopts::value<bool>())
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    result["op"].as<std::string>().c_str();
    std::cout << result["op"].as<std::string>().c_str() << std::endl;
    auto pedFileName = result["ped"].as<std::string>();
    auto vcfFileName = result["vcf"].as<std::string>();
    auto t1 = result["t1"].as<int>();
    auto t2 = result["t2"].as<int>();
    auto bnd = result["bnd"].as<bool>();

    std::ifstream ped(pedFileName);
    std::string fa, mo, ch;
    while (ped >> fa >> mo >> ch) {

    }

    vcflib::VariantCallFile variantFile;

    variantFile.open(vcfFileName);
    if (!variantFile.is_open()) {
        exit(1);
    }
    vcflib::Variant var( variantFile );
    vector<string> samples = variantFile.sampleNames;

    
}
