#include <iostream>
#include <cstring>
#include <fstream>
#include "src/cxxopts.hpp"
#include "src/chromo_phaser.h"
using namespace std;
float T1;
float T2;
int P_ENSURE_BLOCK;
int P_ENSURE_SIDE;
bool ONLY_CHILD;

//TODO this is a tmp
int main(int argc, char *argv[]) {
    cxxopts::Options options("pedhap", "pedigree hap");
    options.add_options()
            ("vcf", "operate: check or solve", cxxopts::value<std::string>())
            ("ped", "If true, across chromo integration will take into consideration", cxxopts::value<std::string>())
            ("out", "Junction database", cxxopts::value<std::string>())
            ("t1", "Input lh file", cxxopts::value<float>()->default_value("5"))
            ("t2", "Checked local hap input file, lh format", cxxopts::value<float>()->default_value("8"))
            ("oc", "Checked local hap input file, lh format", cxxopts::value<bool>()->default_value("1"))
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
    T1 = result["t1"].as<float>();
    T2 = result["t2"].as<float>();
    auto bnd = result["bnd"].as<bool>();
    ONLY_CHILD = result["oc"].as<bool>();
    auto *phaser = new Phaser(vcfFileName, vcfOut);

    std::ifstream ped(pedFileName);
    std::string fa, mo, ch;
    std::vector<std::string> trio;
    bool flag = true;
    while (ped >> ch >> fa >> mo) {
        if(ch == "#" && flag) {
            flag = false;
        }
        trio.push_back(ch);
        trio.push_back(fa);
        trio.push_back(mo);
        if(flag)
            phaser->up_to_down.push_back(trio);
        else
            phaser->down_to_up.push_back(trio);
        trio.clear();
    }
    phaser->phasing();
    delete phaser;
}