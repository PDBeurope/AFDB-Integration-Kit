/**
 * High-throughput C++ implementation of ipSAE calculation.
 * Optimized for ColabFold (AlphaFold2) inputs.
 * 
 * Compile with:
 *   g++ -O3 -march=native -fopenmp -I deps/eigen-3.4.0 -I deps ipsae_cpp.cpp -o ipsae_cpp
 * 
 * Usage:
 *   ./ipsae_cpp <pae_json> <pdb_file> <pae_cutoff> <dist_cutoff>
 *   ./ipsae_cpp --batch <data_dir> <pae_cutoff> <dist_cutoff> [--summary <csv_file>] [--workers N]
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <cstring>

#include <Eigen/Dense>
#include "json.hpp"

#include <omp.h>

using json = nlohmann::json;
using namespace Eigen;
namespace fs = std::filesystem;

// ============================================================================
// Data Structures
// ============================================================================

struct ChainPairResult {
    std::string chain1;
    std::string chain2;
    float ipsae_d0res_asym;
    float ipsae_d0chn_asym;
    float ipsae_d0dom_asym;
    float iptm_d0chn_asym;
    float pDockQ;
    float pDockQ2;
    float LIS;
    int n0res;
    int n0chn;
    int n0dom;
    float d0res;
    float d0chn;
    float d0dom;
    int nres1;
    int nres2;
    int dist_nres1;
    int dist_nres2;
};

struct StructureData {
    std::string pdb_path;
    std::string pae_path;
    MatrixXf cb_coords;  // (n_res, 3)
    MatrixXf pae_matrix; // (n_res, n_res)
    VectorXf plddt;      // (n_res)
    std::vector<char> chain_ids;
    std::vector<int> res_nums;
    std::vector<std::string> unique_chains;
    float iptm_af;
    int n_res;
    bool valid;
};

struct StructureResult {
    std::string pdb_path;
    float pae_cutoff;
    float dist_cutoff;
    float iptm_af;
    std::vector<ChainPairResult> chain_pair_results;
    double processing_time_ms;
};

// ============================================================================
// Helper Functions
// ============================================================================

inline float calc_d0(float L) {
    if (L > 27.0f) {
        return std::max(1.0f, 1.24f * std::pow(L - 15.0f, 1.0f/3.0f) - 1.8f);
    }
    return 1.0f;
}

inline float ptm_func(float pae, float d0) {
    float ratio = pae / d0;
    return 1.0f / (1.0f + ratio * ratio);
}

// ============================================================================
// Parsing Functions
// ============================================================================

StructureData parse_pdb(const std::string& pdb_path) {
    StructureData data;
    data.pdb_path = pdb_path;
    data.valid = false;
    
    std::ifstream file(pdb_path);
    if (!file.is_open()) {
        std::cerr << "Cannot open PDB file: " << pdb_path << std::endl;
        return data;
    }
    
    std::vector<std::array<float, 3>> ca_coords;
    std::vector<std::array<float, 3>> cb_coords;
    std::vector<char> ca_chains;
    std::vector<int> ca_resnums;
    std::unordered_map<std::string, std::array<float, 3>> cb_map;
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) != "ATOM") continue;
        if (line.size() < 54) continue;
        
        std::string atom_name = line.substr(12, 4);
        // Trim whitespace
        atom_name.erase(0, atom_name.find_first_not_of(" "));
        atom_name.erase(atom_name.find_last_not_of(" ") + 1);
        
        std::string res_name = line.substr(17, 3);
        char chain_id = line[21];
        int res_num = std::stoi(line.substr(22, 4));
        float x = std::stof(line.substr(30, 8));
        float y = std::stof(line.substr(38, 8));
        float z = std::stof(line.substr(46, 8));
        
        if (atom_name == "CA") {
            ca_coords.push_back({x, y, z});
            ca_chains.push_back(chain_id);
            ca_resnums.push_back(res_num);
        } else if (atom_name == "CB") {
            std::string key = std::string(1, chain_id) + "_" + std::to_string(res_num);
            cb_map[key] = {x, y, z};
        }
    }
    
    int n_res = ca_coords.size();
    if (n_res == 0) {
        std::cerr << "No CA atoms found in: " << pdb_path << std::endl;
        return data;
    }
    
    // Build CB coordinates (use CA for GLY where CB is missing)
    cb_coords.resize(n_res);
    for (int i = 0; i < n_res; i++) {
        std::string key = std::string(1, ca_chains[i]) + "_" + std::to_string(ca_resnums[i]);
        auto it = cb_map.find(key);
        if (it != cb_map.end()) {
            cb_coords[i] = it->second;
        } else {
            cb_coords[i] = ca_coords[i];  // Use CA for GLY
        }
    }
    
    // Convert to Eigen matrices
    data.cb_coords.resize(n_res, 3);
    for (int i = 0; i < n_res; i++) {
        data.cb_coords(i, 0) = cb_coords[i][0];
        data.cb_coords(i, 1) = cb_coords[i][1];
        data.cb_coords(i, 2) = cb_coords[i][2];
    }
    
    data.chain_ids = ca_chains;
    data.res_nums = ca_resnums;
    data.n_res = n_res;
    
    // Get unique chains in order
    std::set<char> seen;
    for (char c : ca_chains) {
        if (seen.find(c) == seen.end()) {
            data.unique_chains.push_back(std::string(1, c));
            seen.insert(c);
        }
    }
    
    data.valid = true;
    return data;
}

bool parse_json(StructureData& data, const std::string& json_path) {
    data.pae_path = json_path;
    
    std::ifstream file(json_path);
    if (!file.is_open()) {
        std::cerr << "Cannot open JSON file: " << json_path << std::endl;
        return false;
    }
    
    json j;
    try {
        file >> j;
    } catch (const std::exception& e) {
        std::cerr << "JSON parse error in " << json_path << ": " << e.what() << std::endl;
        return false;
    }
    
    // Parse PAE matrix
    std::vector<std::vector<float>> pae_vec;
    if (j.contains("pae")) {
        pae_vec = j["pae"].get<std::vector<std::vector<float>>>();
    } else if (j.contains("predicted_aligned_error")) {
        pae_vec = j["predicted_aligned_error"].get<std::vector<std::vector<float>>>();
    } else {
        std::cerr << "No PAE data in: " << json_path << std::endl;
        return false;
    }
    
    int pae_size = pae_vec.size();
    
    // Validate size match
    if (pae_size != data.n_res) {
        if (pae_size > data.n_res) {
            // Truncate PAE matrix
            pae_size = data.n_res;
        } else {
            std::cerr << "Size mismatch in " << json_path << ": PDB=" << data.n_res 
                      << " PAE=" << pae_vec.size() << std::endl;
            return false;
        }
    }
    
    data.pae_matrix.resize(pae_size, pae_size);
    for (int i = 0; i < pae_size; i++) {
        for (int j = 0; j < pae_size; j++) {
            data.pae_matrix(i, j) = pae_vec[i][j];
        }
    }
    
    // Parse pLDDT
    if (j.contains("plddt")) {
        std::vector<float> plddt_vec = j["plddt"].get<std::vector<float>>();
        data.plddt.resize(std::min((int)plddt_vec.size(), data.n_res));
        for (int i = 0; i < data.plddt.size(); i++) {
            data.plddt(i) = plddt_vec[i];
        }
    } else {
        data.plddt = VectorXf::Zero(data.n_res);
    }
    
    // Parse iptm
    data.iptm_af = j.value("iptm", -1.0f);
    
    return true;
}

// ============================================================================
// Core Computation
// ============================================================================

StructureResult compute_ipsae(const StructureData& data, float pae_cutoff, float dist_cutoff) {
    auto start = std::chrono::high_resolution_clock::now();
    
    StructureResult result;
    result.pdb_path = data.pdb_path;
    result.pae_cutoff = pae_cutoff;
    result.dist_cutoff = dist_cutoff;
    result.iptm_af = data.iptm_af;
    
    int n_res = data.n_res;
    const float pDockQ_cutoff = 8.0f;
    
    // Compute distance matrix using Eigen (SIMD optimized)
    MatrixXf distances(n_res, n_res);
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_res; i++) {
        for (int j = 0; j < n_res; j++) {
            float dx = data.cb_coords(i, 0) - data.cb_coords(j, 0);
            float dy = data.cb_coords(i, 1) - data.cb_coords(j, 1);
            float dz = data.cb_coords(i, 2) - data.cb_coords(j, 2);
            distances(i, j) = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    
    // Pre-compute chain masks
    std::unordered_map<std::string, std::vector<int>> chain_indices;
    std::unordered_map<std::string, int> chain_lengths;
    
    for (const auto& chain : data.unique_chains) {
        chain_indices[chain] = std::vector<int>();
        for (int i = 0; i < n_res; i++) {
            if (data.chain_ids[i] == chain[0]) {
                chain_indices[chain].push_back(i);
            }
        }
        chain_lengths[chain] = chain_indices[chain].size();
    }
    
    // Process each chain pair
    for (const auto& chain1 : data.unique_chains) {
        for (const auto& chain2 : data.unique_chains) {
            if (chain1 == chain2) continue;
            
            const auto& idx1 = chain_indices[chain1];
            const auto& idx2 = chain_indices[chain2];
            
            int n0chn = chain_lengths[chain1] + chain_lengths[chain2];
            float d0chn = calc_d0(n0chn);
            
            // Pre-compute PTM matrix for d0chn
            MatrixXf ptm_matrix_d0chn(n_res, n_res);
            #pragma omp parallel for schedule(static) collapse(2)
            for (int i = 0; i < n_res; i++) {
                for (int j = 0; j < n_res; j++) {
                    ptm_matrix_d0chn(i, j) = ptm_func(data.pae_matrix(i, j), d0chn);
                }
            }
            
            // ================== pDockQ ==================
            std::vector<bool> contact_mask(n_res, false);
            int n_pairs = 0;
            
            for (int i : idx1) {
                for (int j : idx2) {
                    if (distances(i, j) <= pDockQ_cutoff) {
                        contact_mask[i] = true;
                        contact_mask[j] = true;
                        n_pairs++;
                    }
                }
            }
            
            float pDockQ_val = 0.0f;
            if (n_pairs > 0) {
                float sum_plddt = 0.0f;
                int count = 0;
                for (int i = 0; i < n_res; i++) {
                    if (contact_mask[i]) {
                        sum_plddt += data.plddt(i);
                        count++;
                    }
                }
                float mean_plddt = sum_plddt / count;
                float x = mean_plddt * std::log10(n_pairs);
                pDockQ_val = 0.724f / (1.0f + std::exp(-0.052f * (x - 152.611f))) + 0.018f;
            }
            
            // ================== pDockQ2 ==================
            float pDockQ2_val = 0.0f;
            if (n_pairs > 0) {
                float ptm_sum = 0.0f;
                int ptm_count = 0;
                
                for (int i : idx1) {
                    for (int j : idx2) {
                        if (distances(i, j) <= pDockQ_cutoff) {
                            ptm_sum += ptm_func(data.pae_matrix(i, j), 10.0f);
                            ptm_count++;
                        }
                    }
                }
                
                float mean_ptm = ptm_sum / ptm_count;
                float sum_plddt = 0.0f;
                int count = 0;
                for (int i = 0; i < n_res; i++) {
                    if (contact_mask[i]) {
                        sum_plddt += data.plddt(i);
                        count++;
                    }
                }
                float mean_plddt = sum_plddt / count;
                float x = mean_plddt * mean_ptm;
                pDockQ2_val = 1.31f / (1.0f + std::exp(-0.075f * (x - 84.733f))) + 0.005f;
            }
            
            // ================== LIS ==================
            float LIS_val = 0.0f;
            float lis_sum = 0.0f;
            int lis_count = 0;
            
            for (int i : idx1) {
                for (int j : idx2) {
                    float pae = data.pae_matrix(i, j);
                    if (pae < 12.0f) {
                        lis_sum += (12.0f - pae) / 12.0f;
                        lis_count++;
                    }
                }
            }
            
            if (lis_count > 0) {
                LIS_val = lis_sum / lis_count;
            }
            
            // ================== ipTM/ipSAE ==================
            VectorXf iptm_d0chn_byres = VectorXf::Zero(n_res);
            VectorXf ipsae_d0chn_byres = VectorXf::Zero(n_res);
            
            std::set<int> unique_res_chain1, unique_res_chain2;
            std::set<int> dist_unique_res_chain1, dist_unique_res_chain2;
            
            for (int i : idx1) {
                // ipTM: mean over all chain2 residues
                float iptm_sum = 0.0f;
                for (int j : idx2) {
                    iptm_sum += ptm_matrix_d0chn(i, j);
                }
                iptm_d0chn_byres(i) = iptm_sum / idx2.size();
                
                // ipSAE: mean over chain2 residues with PAE < cutoff
                float ipsae_sum = 0.0f;
                int ipsae_count = 0;
                
                for (int j : idx2) {
                    if (data.pae_matrix(i, j) < pae_cutoff) {
                        ipsae_sum += ptm_matrix_d0chn(i, j);
                        ipsae_count++;
                        unique_res_chain1.insert(data.res_nums[i]);
                        unique_res_chain2.insert(data.res_nums[j]);
                        
                        if (distances(i, j) < dist_cutoff) {
                            dist_unique_res_chain1.insert(data.res_nums[i]);
                            dist_unique_res_chain2.insert(data.res_nums[j]);
                        }
                    }
                }
                
                if (ipsae_count > 0) {
                    ipsae_d0chn_byres(i) = ipsae_sum / ipsae_count;
                }
            }
            
            int n0dom = unique_res_chain1.size() + unique_res_chain2.size();
            float d0dom = calc_d0(n0dom);
            
            // Compute d0dom-based scores
            VectorXf ipsae_d0dom_byres = VectorXf::Zero(n_res);
            VectorXf n0res_byres = VectorXf::Zero(n_res);
            VectorXf d0res_byres = VectorXf::Ones(n_res);
            VectorXf ipsae_d0res_byres = VectorXf::Zero(n_res);
            
            MatrixXf ptm_matrix_d0dom(n_res, n_res);
            if (n0dom > 0) {
                #pragma omp parallel for schedule(static) collapse(2)
                for (int i = 0; i < n_res; i++) {
                    for (int j = 0; j < n_res; j++) {
                        ptm_matrix_d0dom(i, j) = ptm_func(data.pae_matrix(i, j), d0dom);
                    }
                }
            }
            
            for (int i : idx1) {
                int valid_count = 0;
                float dom_sum = 0.0f;
                
                for (int j : idx2) {
                    if (data.pae_matrix(i, j) < pae_cutoff) {
                        valid_count++;
                        if (n0dom > 0) {
                            dom_sum += ptm_matrix_d0dom(i, j);
                        }
                    }
                }
                
                n0res_byres(i) = valid_count;
                
                if (valid_count > 0) {
                    if (n0dom > 0) {
                        ipsae_d0dom_byres(i) = dom_sum / valid_count;
                    }
                    
                    float d0_res = calc_d0(valid_count);
                    d0res_byres(i) = d0_res;
                    
                    float res_sum = 0.0f;
                    for (int j : idx2) {
                        if (data.pae_matrix(i, j) < pae_cutoff) {
                            res_sum += ptm_func(data.pae_matrix(i, j), d0_res);
                        }
                    }
                    ipsae_d0res_byres(i) = res_sum / valid_count;
                }
            }
            
            // Get max indices
            int max_idx_iptm = 0, max_idx_ipsae_chn = 0, max_idx_ipsae_dom = 0, max_idx_ipsae_res = 0;
            float max_val;
            
            max_val = iptm_d0chn_byres.maxCoeff(&max_idx_iptm);
            max_val = ipsae_d0chn_byres.maxCoeff(&max_idx_ipsae_chn);
            max_val = ipsae_d0dom_byres.maxCoeff(&max_idx_ipsae_dom);
            max_val = ipsae_d0res_byres.maxCoeff(&max_idx_ipsae_res);
            
            ChainPairResult cpr;
            cpr.chain1 = chain1;
            cpr.chain2 = chain2;
            cpr.ipsae_d0res_asym = ipsae_d0res_byres(max_idx_ipsae_res);
            cpr.ipsae_d0chn_asym = ipsae_d0chn_byres(max_idx_ipsae_chn);
            cpr.ipsae_d0dom_asym = ipsae_d0dom_byres(max_idx_ipsae_dom);
            cpr.iptm_d0chn_asym = iptm_d0chn_byres(max_idx_iptm);
            cpr.pDockQ = pDockQ_val;
            cpr.pDockQ2 = pDockQ2_val;
            cpr.LIS = LIS_val;
            cpr.n0res = (int)n0res_byres(max_idx_ipsae_res);
            cpr.n0chn = n0chn;
            cpr.n0dom = n0dom;
            cpr.d0res = d0res_byres(max_idx_ipsae_res);
            cpr.d0chn = d0chn;
            cpr.d0dom = d0dom;
            cpr.nres1 = unique_res_chain1.size();
            cpr.nres2 = unique_res_chain2.size();
            cpr.dist_nres1 = dist_unique_res_chain1.size();
            cpr.dist_nres2 = dist_unique_res_chain2.size();
            
            result.chain_pair_results.push_back(cpr);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    result.processing_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    return result;
}

// ============================================================================
// Output Functions
// ============================================================================

void write_txt_output(const StructureResult& result) {
    std::string pdb_stem = result.pdb_path;
    size_t pos = pdb_stem.rfind(".pdb");
    if (pos != std::string::npos) {
        pdb_stem = pdb_stem.substr(0, pos);
    }
    
    char pae_str[3], dist_str[3];
    snprintf(pae_str, 3, "%02d", (int)result.pae_cutoff);
    snprintf(dist_str, 3, "%02d", (int)result.dist_cutoff);
    
    std::string out_path = pdb_stem + "_" + pae_str + "_" + dist_str + ".txt";
    
    std::ofstream out(out_path);
    out << "\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n";
    
    // Group by chain pairs
    std::map<std::pair<std::string, std::string>, std::vector<ChainPairResult>> pair_results;
    for (const auto& r : result.chain_pair_results) {
        auto key = std::make_pair(std::min(r.chain1, r.chain2), std::max(r.chain1, r.chain2));
        pair_results[key].push_back(r);
    }
    
    for (const auto& [key, results] : pair_results) {
        for (const auto& r : results) {
            char line[512];
            snprintf(line, sizeof(line),
                "%s    %s     %3s  %3s  asym  %8.6f    %8.6f    %8.6f    %5.3f    %8.6f    %8.4f   %8.4f   %8.4f   %5d  %5d  %5d  %6.2f  %6.2f  %6.2f  %5d   %5d   %5d   %5d   %s\n",
                r.chain1.c_str(), r.chain2.c_str(), pae_str, dist_str,
                r.ipsae_d0res_asym, r.ipsae_d0chn_asym, r.ipsae_d0dom_asym,
                result.iptm_af, r.iptm_d0chn_asym, r.pDockQ, r.pDockQ2, r.LIS,
                r.n0res, r.n0chn, r.n0dom, r.d0res, r.d0chn, r.d0dom,
                r.nres1, r.nres2, r.dist_nres1, r.dist_nres2, pdb_stem.c_str());
            out << line;
        }
        
        // Write max line
        if (results.size() == 2) {
            const auto& r1 = results[0];
            const auto& r2 = results[1];
            
            float max_ipsae_d0res = std::max(r1.ipsae_d0res_asym, r2.ipsae_d0res_asym);
            int max_n0res = (r1.ipsae_d0res_asym >= r2.ipsae_d0res_asym) ? r1.n0res : r2.n0res;
            float max_d0res = (r1.ipsae_d0res_asym >= r2.ipsae_d0res_asym) ? r1.d0res : r2.d0res;
            float max_ipsae_d0chn = std::max(r1.ipsae_d0chn_asym, r2.ipsae_d0chn_asym);
            float max_ipsae_d0dom = std::max(r1.ipsae_d0dom_asym, r2.ipsae_d0dom_asym);
            int max_n0dom = (r1.ipsae_d0dom_asym >= r2.ipsae_d0dom_asym) ? r1.n0dom : r2.n0dom;
            float max_d0dom = (r1.ipsae_d0dom_asym >= r2.ipsae_d0dom_asym) ? r1.d0dom : r2.d0dom;
            float max_iptm = std::max(r1.iptm_d0chn_asym, r2.iptm_d0chn_asym);
            float avg_LIS = (r1.LIS + r2.LIS) / 2.0f;
            float max_pDockQ2 = std::max(r1.pDockQ2, r2.pDockQ2);
            
            char line[512];
            snprintf(line, sizeof(line),
                "%s    %s     %3s  %3s  max   %8.6f    %8.6f    %8.6f    %5.3f    %8.6f    %8.4f   %8.4f   %8.4f   %5d  %5d  %5d  %6.2f  %6.2f  %6.2f  %5d   %5d   %5d   %5d   %s\n",
                key.first.c_str(), key.second.c_str(), pae_str, dist_str,
                max_ipsae_d0res, max_ipsae_d0chn, max_ipsae_d0dom,
                result.iptm_af, max_iptm, r1.pDockQ, max_pDockQ2, avg_LIS,
                max_n0res, r1.n0chn, max_n0dom, max_d0res, r1.d0chn, max_d0dom,
                std::max(r1.nres1, r2.nres2), std::max(r1.nres2, r2.nres1),
                std::max(r1.dist_nres1, r2.dist_nres2), std::max(r1.dist_nres2, r2.dist_nres1),
                pdb_stem.c_str());
            out << line;
        }
    }
    
    out << "\n";
}

void write_csv_summary(const std::vector<StructureResult>& results, const std::string& csv_path) {
    // Collect all unique ordered pairs (first < second) across all results
    std::set<std::string> all_unique_pairs;
    for (const auto& result : results) {
        for (const auto& r : result.chain_pair_results) {
            std::string pair = (r.chain1 < r.chain2)
                ? r.chain1 + r.chain2
                : r.chain2 + r.chain1;
            all_unique_pairs.insert(pair);
        }
    }
    std::vector<std::string> pairs(all_unique_pairs.begin(), all_unique_pairs.end());
    std::sort(pairs.begin(), pairs.end());

    // Write header
    std::ofstream out(csv_path);
    out << "pdb_path,pae_cutoff,dist_cutoff,iptm_af";
    for (const auto& p : pairs) {
        std::string ab = p;
        std::string ba = std::string(1, p[1]) + std::string(1, p[0]);
        // Direction-dependent metrics
        out << ",ipsae_" << ab << ",ipsae_" << ba
            << ",ipsae_d0chn_" << ab << ",ipsae_d0chn_" << ba
            << ",ipsae_d0dom_" << ab << ",ipsae_d0dom_" << ba
            << ",iptm_d0chn_" << ab << ",iptm_d0chn_" << ba
            << ",pDockQ2_" << ab << ",pDockQ2_" << ba
            << ",LIS_" << ab << ",LIS_" << ba
            << ",n0res_" << ab << ",n0res_" << ba
            << ",n0dom_" << ab << ",n0dom_" << ba
            << ",d0res_" << ab << ",d0res_" << ba
            << ",d0dom_" << ab << ",d0dom_" << ba
            << ",nres1_" << ab << ",nres1_" << ba
            << ",nres2_" << ab << ",nres2_" << ba
            << ",dist_nres1_" << ab << ",dist_nres1_" << ba
            << ",dist_nres2_" << ab << ",dist_nres2_" << ba;
        // Symmetric metrics (same for both directions)
        out << ",pDockQ"
            << ",n0chn"
            << ",d0chn";
    }
    out << ",processing_time_ms\n";

    // Write data
    for (const auto& result : results) {
        // Build lookup: "AB" -> &ChainPairResult for A->B direction
        std::map<std::string, const ChainPairResult*> lookup;
        for (const auto& r : result.chain_pair_results) {
            lookup[r.chain1 + r.chain2] = &r;
        }

        out << result.pdb_path << "," << result.pae_cutoff << "," << result.dist_cutoff << "," << result.iptm_af;

        for (const auto& p : pairs) {
            std::string ab = p;
            std::string ba = std::string(1, p[1]) + std::string(1, p[0]);
            auto it_ab = lookup.find(ab);
            auto it_ba = lookup.find(ba);
            bool has_ab = it_ab != lookup.end();
            bool has_ba = it_ba != lookup.end();

            // Direction-dependent metrics: _AB then _BA
            // ipsae
            out << ","; if (has_ab) out << it_ab->second->ipsae_d0res_asym;
            out << ","; if (has_ba) out << it_ba->second->ipsae_d0res_asym;
            // ipsae_d0chn
            out << ","; if (has_ab) out << it_ab->second->ipsae_d0chn_asym;
            out << ","; if (has_ba) out << it_ba->second->ipsae_d0chn_asym;
            // ipsae_d0dom
            out << ","; if (has_ab) out << it_ab->second->ipsae_d0dom_asym;
            out << ","; if (has_ba) out << it_ba->second->ipsae_d0dom_asym;
            // iptm_d0chn
            out << ","; if (has_ab) out << it_ab->second->iptm_d0chn_asym;
            out << ","; if (has_ba) out << it_ba->second->iptm_d0chn_asym;
            // pDockQ2
            out << ","; if (has_ab) out << it_ab->second->pDockQ2;
            out << ","; if (has_ba) out << it_ba->second->pDockQ2;
            // LIS
            out << ","; if (has_ab) out << it_ab->second->LIS;
            out << ","; if (has_ba) out << it_ba->second->LIS;
            // n0res
            out << ","; if (has_ab) out << it_ab->second->n0res;
            out << ","; if (has_ba) out << it_ba->second->n0res;
            // n0dom
            out << ","; if (has_ab) out << it_ab->second->n0dom;
            out << ","; if (has_ba) out << it_ba->second->n0dom;
            // d0res
            out << ","; if (has_ab) out << it_ab->second->d0res;
            out << ","; if (has_ba) out << it_ba->second->d0res;
            // d0dom
            out << ","; if (has_ab) out << it_ab->second->d0dom;
            out << ","; if (has_ba) out << it_ba->second->d0dom;
            // nres1
            out << ","; if (has_ab) out << it_ab->second->nres1;
            out << ","; if (has_ba) out << it_ba->second->nres1;
            // nres2
            out << ","; if (has_ab) out << it_ab->second->nres2;
            out << ","; if (has_ba) out << it_ba->second->nres2;
            // dist_nres1
            out << ","; if (has_ab) out << it_ab->second->dist_nres1;
            out << ","; if (has_ba) out << it_ba->second->dist_nres1;
            // dist_nres2
            out << ","; if (has_ab) out << it_ab->second->dist_nres2;
            out << ","; if (has_ba) out << it_ba->second->dist_nres2;

            // Symmetric metrics (use either direction, they're identical)
            // pDockQ
            out << ","; if (has_ab) out << it_ab->second->pDockQ; else if (has_ba) out << it_ba->second->pDockQ;
            // n0chn
            out << ","; if (has_ab) out << it_ab->second->n0chn; else if (has_ba) out << it_ba->second->n0chn;
            // d0chn
            out << ","; if (has_ab) out << it_ab->second->d0chn; else if (has_ba) out << it_ba->second->d0chn;
        }

        out << "," << result.processing_time_ms << "\n";
    }
}

// ============================================================================
// Batch Processing
// ============================================================================

std::vector<std::pair<std::string, std::string>> find_paired_files(const std::string& data_dir) {
    std::vector<std::pair<std::string, std::string>> pairs;
    
    for (const auto& entry : fs::directory_iterator(data_dir)) {
        std::string path = entry.path().string();
        if (path.find("-meta_v1.json") != std::string::npos) {
            std::string pdb_path = path;
            size_t pos = pdb_path.find("-meta_v1.json");
            pdb_path.replace(pos, 13, "-model_v1.pdb");
            
            if (fs::exists(pdb_path)) {
                pairs.push_back({pdb_path, path});
            }
        }
    }
    
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::vector<StructureResult> process_batch(
    const std::vector<std::pair<std::string, std::string>>& file_pairs,
    float pae_cutoff, float dist_cutoff, int n_workers, bool quiet) {
    
    std::vector<StructureResult> results(file_pairs.size());
    std::atomic<int> completed{0};
    std::atomic<int> errors{0};
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    #pragma omp parallel for schedule(dynamic) num_threads(n_workers)
    for (size_t i = 0; i < file_pairs.size(); i++) {
        const auto& [pdb_path, pae_path] = file_pairs[i];
        
        StructureData data = parse_pdb(pdb_path);
        if (!data.valid) {
            errors++;
            completed++;
            continue;
        }
        
        if (!parse_json(data, pae_path)) {
            errors++;
            completed++;
            continue;
        }
        
        results[i] = compute_ipsae(data, pae_cutoff, dist_cutoff);
        
        int c = ++completed;
        if (!quiet && c % 100 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - total_start).count();
            double rate = c / elapsed;
            double eta = (file_pairs.size() - c) / rate;
            #pragma omp critical
            {
                std::cout << "  Progress: " << c << "/" << file_pairs.size() 
                          << " (" << std::fixed << std::setprecision(1) << rate << "/s, ETA: " 
                          << (int)eta << "s)" << std::endl;
            }
        }
    }
    
    // Filter out empty results
    std::vector<StructureResult> valid_results;
    for (const auto& r : results) {
        if (!r.pdb_path.empty()) {
            valid_results.push_back(r);
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(total_end - total_start).count();
    
    if (!quiet) {
        std::cout << "\nCompleted " << valid_results.size() << "/" << file_pairs.size() 
                  << " structures in " << std::fixed << std::setprecision(1) << total_time << "s" << std::endl;
        std::cout << "Throughput: " << std::fixed << std::setprecision(1) 
                  << valid_results.size() / total_time << " structures/sec" << std::endl;
    }
    
    return valid_results;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <pae_json> <pdb_file> <pae_cutoff> <dist_cutoff>" << std::endl;
        std::cerr << "   or: " << argv[0] << " --batch <data_dir> <pae_cutoff> <dist_cutoff> [--summary <csv>] [--workers N] [--quiet]" << std::endl;
        return 1;
    }
    
    if (std::string(argv[1]) == "--batch") {
        if (argc < 5) {
            std::cerr << "Usage: " << argv[0] << " --batch <data_dir> <pae_cutoff> <dist_cutoff> [options]" << std::endl;
            return 1;
        }
        
        std::string data_dir = argv[2];
        float pae_cutoff = std::stof(argv[3]);
        float dist_cutoff = std::stof(argv[4]);
        
        std::string summary_csv = "";
        int n_workers = omp_get_max_threads();
        bool quiet = false;
        bool no_individual = false;
        
        for (int i = 5; i < argc; i++) {
            if (std::string(argv[i]) == "--summary" && i + 1 < argc) {
                summary_csv = argv[++i];
            } else if (std::string(argv[i]) == "--workers" && i + 1 < argc) {
                n_workers = std::stoi(argv[++i]);
            } else if (std::string(argv[i]) == "--quiet") {
                quiet = true;
            } else if (std::string(argv[i]) == "--no-individual") {
                no_individual = true;
            }
        }
        
        auto file_pairs = find_paired_files(data_dir);
        
        if (!quiet) {
            std::cout << "Processing " << file_pairs.size() << " structures with " 
                      << n_workers << " workers..." << std::endl;
        }
        
        auto results = process_batch(file_pairs, pae_cutoff, dist_cutoff, n_workers, quiet);
        
        if (!no_individual) {
            for (const auto& result : results) {
                write_txt_output(result);
            }
        }
        
        if (!summary_csv.empty()) {
            write_csv_summary(results, summary_csv);
            if (!quiet) {
                std::cout << "Summary written to: " << summary_csv << std::endl;
            }
        }
        
    } else {
        // Single file mode
        std::string pae_path = argv[1];
        std::string pdb_path = argv[2];
        float pae_cutoff = std::stof(argv[3]);
        float dist_cutoff = std::stof(argv[4]);
        
        StructureData data = parse_pdb(pdb_path);
        if (!data.valid) {
            return 1;
        }
        
        if (!parse_json(data, pae_path)) {
            return 1;
        }
        
        StructureResult result = compute_ipsae(data, pae_cutoff, dist_cutoff);
        write_txt_output(result);
    }
    
    return 0;
}
