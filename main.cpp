#include <iostream>
#include <cstring>
#include <fstream>
#include "src/cxxopts.hpp"
#include "vcflib/Variant.h"
#include "src/phaser.h"
using namespace std;

//TODO this is a tmp
int main(int argc, char *argv[]) {
    cxxopts::Options options("pedhap", "pedigree hap");
    options.add_options()
            ("vcf", "operate: check or solve", cxxopts::value<std::string>())
            ("ped", "If true, across chromo integration will take into consideration", cxxopts::value<std::string>())
            ("out", "Junction database", cxxopts::value<std::string>())
            ("t1", "Input lh file", cxxopts::value<float>()->default_value("0.6"))
            ("t2", "Checked local hap input file, lh format", cxxopts::value<float>()->default_value("0"))
            ("bnd", "ILP out file prefix, only for check", cxxopts::value<bool>()->default_value("1"))
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    auto pedFileName = result["ped"].as<std::string>();
    auto vcfFileName = result["vcf"].as<std::string>();
    auto vcfOut = result["out"].as<std::string>();
    auto t1 = result["t1"].as<float>();
    auto t2 = result["t2"].as<float>();
    auto bnd = result["bnd"].as<bool>();

    std::ifstream ped(pedFileName);
    std::string fa, mo, ch;
    while (ped >> fa >> mo >> ch) {

    }

//    vcflib::VariantCallFile variantFile;
//
//    variantFile.open(vcfFileName);
//    if (!variantFile.is_open()) {
//        exit(1);
//    }
//    vcflib::Variant var( variantFile );
//    vector<string> samples = variantFile.sampleNames;

    Phaser *phaser = new Phaser(vcfFileName, vcfOut);
    phaser->phasing();
    delete phaser;
}
