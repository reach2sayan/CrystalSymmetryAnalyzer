#ifndef __SYMMETRY_ANALYZER_BONDS_HPP__
#define __SYMMETRY_ANALYZER_BONDS_HPP__

#include <string>
#include <unordered_map>

const std::unordered_map<std::string, double> bond_lengths = {
    {"C-C", 1.89},   {"C-H", 1.25},   {"C-O", 1.60},   {"C-N", 1.80},
    {"C-S", 2.00},   {"C-P", 1.95},   {"C-Si", 2.00},  {"C-F", 1.45},
    {"C-Cl", 1.90},  {"C-Br", 2.10},  {"C-I", 2.30},   {"H-C", 1.25},
    {"H-H", 0.80},   {"H-O", 1.15},   {"H-N", 1.22},   {"H-S", 1.70},
    {"H-P", 1.50},   {"H-Si", 1.50},  {"H-F", 0.90},   {"H-Cl", 1.40},
    {"H-Br", 1.60},  {"H-I", 1.80},   {"O-C", 1.60},   {"O-H", 1.15},
    {"O-O", 1.50},   {"O-N", 1.55},   {"O-S", 2.00},   {"O-P", 1.90},
    {"O-Si", 2.00},  {"O-F", 1.40},   {"O-Cl", 1.80},  {"O-Br", 2.00},
    {"O-I", 2.30},   {"N-C", 1.80},   {"N-H", 1.22},   {"N-O", 1.55},
    {"N-N", 1.56},   {"N-S", 2.05},   {"N-P", 1.95},   {"N-Si", 1.95},
    {"N-F", 1.35},   {"N-Cl", 1.80},  {"N-Br", 1.98},  {"N-I", 2.17},
    {"S-C", 2.00},   {"S-H", 1.79},   {"S-O", 2.00},   {"S-N", 2.05},
    {"S-S", 2.50},   {"S-P", 2.32},   {"S-Si", 2.36},  {"S-F", 1.82},
    {"S-Cl", 2.27},  {"S-Br", 2.45},  {"S-I", 2.64},   {"P-C", 1.95},
    {"P-H", 1.50},   {"P-O", 1.90},   {"P-N", 1.95},   {"P-S", 2.32},
    {"P-P", 2.64},   {"P-Si", 2.18},  {"P-F", 1.64},   {"P-Cl", 2.09},
    {"P-Br", 2.27},  {"P-I", 2.46},   {"Si-C", 2.00},  {"Si-H", 1.50},
    {"Si-O", 2.00},  {"Si-N", 1.95},  {"Si-S", 2.36},  {"Si-P", 2.18},
    {"Si-Si", 2.22}, {"Si-F", 1.68},  {"Si-Cl", 2.13}, {"Si-Br", 2.31},
    {"Si-I", 2.50},  {"F-C", 1.45},   {"F-H", 0.90},   {"F-O", 1.40},
    {"F-N", 1.35},   {"F-S", 1.82},   {"F-P", 1.64},   {"F-Si", 1.68},
    {"F-F", 1.14},   {"F-Cl", 1.59},  {"F-Br", 1.77},  {"F-I", 1.96},
    {"Cl-C", 1.90},  {"Cl-H", 1.40},  {"Cl-O", 1.80},  {"Cl-N", 1.80},
    {"Cl-S", 2.27},  {"Cl-P", 2.09},  {"Cl-Si", 2.13}, {"Cl-F", 1.59},
    {"Cl-Cl", 2.04}, {"Cl-Br", 2.22}, {"Cl-I", 2.41},  {"Br-C", 2.10},
    {"Br-H", 1.60},  {"Br-O", 2.00},  {"Br-N", 1.98},  {"Br-S", 2.45},
    {"Br-P", 2.27},  {"Br-Si", 2.31}, {"Br-F", 1.77},  {"Br-Cl", 2.22},
    {"Br-Br", 2.40}, {"Br-I", 2.59},  {"I-C", 2.30},   {"I-H", 1.80},
    {"I-O", 2.30},   {"I-N", 2.17},   {"I-S", 2.64},   {"I-P", 2.46},
    {"I-Si", 2.50},  {"I-F", 1.96},   {"I-Cl", 2.41},  {"I-Br", 2.59},
    {"I-I", 2.78}};

#endif
