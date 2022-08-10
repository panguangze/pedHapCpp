#include <iostream>
#include <cstring>
#include <fstream>
#include "src/cxxopts.hpp"
#include "src/chromo_phaser.h"
#include "src/util.h"
using namespace std;
float T1;
float T2;
int P_ENSURE_BLOCK;
int P_ENSURE_SIDE;
int CORRECT_SCORE = 2;
int ERROR_SCORE = 1;
bool ONLY_CHILD;
bool XY;
bool IS_MALE;
bool IS_DEBUG;
int BND_RANGE;
int MIN_SNP_RECOM;
bool HETE;
int CHILD = 0;

//TODO this is a tmp
int main(int argc, char *argv[]) {
    cxxopts::Options options("pedhap", "pedigree hap");
    options.add_options()
            ("vcf", "sorted vcf, gz", cxxopts::value<std::string>())
            ("ped", "Pedigree file generated by python file", cxxopts::value<std::string>())
            ("out", "Output phased vcf file", cxxopts::value<std::string>())
            ("t1", "Block Merge Support site count for hete", cxxopts::value<float>()->default_value("2"))
            ("t2", "Block Merge Support site count for homo", cxxopts::value<float>()->default_value("2"))
            ("oc", "Only child with be phased", cxxopts::value<bool>()->default_value("false"))
            ("bnd_range", "Phasing bnd position range", cxxopts::value<int>()->default_value("100"))
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("xy", "Whether phasing sexual chromosome", cxxopts::value<bool>()->default_value("false"))
            ("sexual", "Sexual of child", cxxopts::value<std::string>())
            ("debug", "debug", cxxopts::value<bool>()->default_value("true"))
            ("hete", "hete", cxxopts::value<bool>()->default_value("true"))
            ("child", "child", cxxopts::value<int>()->default_value("0"))
            ("seg_dup", "Reference segment duplication regions", cxxopts::value<std::string>())
            ("ref", "Reference fasta", cxxopts::value<std::string>())
            ("contigs", "phase contigs", cxxopts::value<std::string>())
//            ("snp_min_support", "Minimum continuous snp site support for home recombination event", cxxopts::value<int>()->default_value("2"))
            ("homo_recom", "Output homology recombination events to a file", cxxopts::value<std::string>())
            ("help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    std::vector<std::string> contig_vec;
    if (result.count("contigs")) {
        contig_vec = split(result["contigs"].as<std::string>(), ',');
    }
    std::string seg_dup;
    if (result.count("seg_dup")) {
        seg_dup = result["seg_dup"].as<std::string>();
    }
    if (result.count("debug")) {
        IS_DEBUG = result["debug"].as<bool>();
    }
    if (result.count("hete")) {
        HETE = result["hete"].as<bool>();
    }
    if (result.count("child")) {
        CHILD = result["child"].as<int>();
    }
    std::string homoOut;
    if (result.count("homo_recom")) {
        homoOut = result["homo_recom"].as<std::string>();
    }
    auto pedFileName = result["ped"].as<std::string>();
    auto vcfFileName = result["vcf"].as<std::string>();
    auto vcfOut = result["out"].as<std::string>();
    XY = result["xy"].as<bool>();
    T1 = result["t1"].as<float>();
    T2 = result["t2"].as<float>();
    BND_RANGE = result["bnd_range"].as<int>();
//    MIN_SNP_RECOM = result["snp_min_support"].as<int>();
    ONLY_CHILD = result["oc"].as<bool>();
    auto *phaser = new Phaser(vcfFileName, vcfOut, homoOut);
    phaser->set_contigs(contig_vec);
    phaser->set_segdup(seg_dup);

    std::ifstream ped(pedFileName);
    std::string fa, mo, ch, gender;
    std::vector<std::string> trio;
    bool flag = true;
    while (ped >> ch >> fa >> mo) {
        if(ch == "#" && flag) {
            flag = false;
            continue;
        }
        trio.push_back(ch);
        trio.push_back(fa);
        trio.push_back(mo);
        trio.push_back(gender);
        phaser->set_trio(trio,flag);
//        if(flag)
//            phaser->set_trio(trio,flag);
//        else
//            phaser->down_to_up.push_back(trio);
        trio.clear();
    }
    phaser->phasing();
    delete phaser;
}