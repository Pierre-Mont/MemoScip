#include "memoscip_reader.h"

#include "scip/scip.h"
#include "scip/scip_cons.h"
#include "scip/cons_linear.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <lzma.h>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

static inline string trim(const string& s)
{
    size_t a = 0;
    while (a < s.size() && isspace((unsigned char)s[a])) ++a;
    size_t b = s.size();
    while (b > a && isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

// Custom streambuf for LZMA decompression
class LZMAStreambuf : public std::streambuf {
private:
    FILE* file_;
    lzma_stream stream_;
    std::vector<uint8_t> in_buf_;
    std::vector<uint8_t> out_buf_;
    static const size_t BUF_SIZE = 65536;

    int fill_buffer() {
        if (in_buf_.empty()) in_buf_.resize(BUF_SIZE);
        if (out_buf_.empty()) out_buf_.resize(BUF_SIZE);
        
        size_t nread = fread(in_buf_.data(), 1, BUF_SIZE, file_);
        if (nread == 0 && ferror(file_)) return -1;
        
        stream_.next_in = in_buf_.data();
        stream_.avail_in = nread;
        stream_.next_out = out_buf_.data();
        stream_.avail_out = out_buf_.size();
        
        lzma_ret ret = lzma_code(&stream_, nread == 0 ? LZMA_FINISH : LZMA_RUN);
        if (ret != LZMA_OK && ret != LZMA_STREAM_END) return -1;
        
        size_t produced = out_buf_.size() - stream_.avail_out;
        setg((char*)out_buf_.data(), (char*)out_buf_.data(), (char*)out_buf_.data() + produced);
        
        return produced > 0 ? *eback() : EOF;
    }

public:
    LZMAStreambuf(FILE* f) : file_(f), stream_(LZMA_STREAM_INIT) {
        lzma_ret ret = lzma_stream_decoder(&stream_, UINT64_MAX, LZMA_CONCATENATED);
        if (ret != LZMA_OK) throw std::runtime_error("lzma_stream_decoder failed");
        fill_buffer();
    }
    
    ~LZMAStreambuf() {
        lzma_end(&stream_);
    }

    int underflow() override {
        if (gptr() < egptr()) return *gptr();
        return fill_buffer();
    }
};

class LZMAIStream : public std::istream {
private:
    std::unique_ptr<LZMAStreambuf> buf_;
public:
    LZMAIStream(FILE* f) : std::istream(nullptr) {
        buf_ = std::make_unique<LZMAStreambuf>(f);
        rdbuf(buf_.get());
    }
};

std::unique_ptr<std::istream> open_instance_file(const std::string& path, std::string& err) {
    if (path.length() >= 3 && path.substr(path.length() - 3) == ".xz") {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) {
            err = "Failed to open .xz file: " + path;
            return nullptr;
        }
        try {
            return std::make_unique<LZMAIStream>(f);
        } catch (const std::exception& e) {
            fclose(f);
            err = std::string("LZMA error: ") + e.what();
            return nullptr;
        }
    } else {
        auto ifs = std::make_unique<std::ifstream>(path);
        if (!ifs->is_open()) {
            err = "Failed to open file: " + path;
            return nullptr;
        }
        return ifs;
    }
}


bool readKPCGInstance(SCIP* scip, const std::string& path, std::string& err)
{
    if (scip == nullptr) {
        err = "SCIP pointer is null";
        return false;
    }

    auto stream = open_instance_file(path, err);
    if (!stream) {
        return false;
    }
    std::istream& ifs = *stream;

    int n_items = -1;
    int capacity = -1;
    vector<tuple<int,int,int>> items;
    vector<pair<int,int>> edges;
    string mode;

    regex param_re(R"(^\s*param\s+(\w+)\s*:=\s*(\d+)\s*;\s*$)");

    string raw;
    while (std::getline(ifs, raw)) {
        string line = trim(raw);
        if (line.empty() || line[0] == '#')
            continue;

        smatch m;
        if (regex_match(line, m, param_re)) {
            string name = m[1];
            int val = stoi(m[2]);
            if (name == "n") n_items = val;
            else if (name == "c") capacity = val;
            continue;
        }

        if (line.rfind("param :", 0) == 0) {
            mode = "items";
            continue;
        }
        if (line.rfind("set E", 0) == 0) {
            mode = "edges";
            continue;
        }

        if (line == ";") {
            if (mode == "items") mode.clear();
            else if (mode == "edges") break;
            continue;
        }

        if (mode == "items") {
            // parse: idx profit weight
            istringstream iss(line);
            int idx, profit, weight;
            if (!(iss >> idx >> profit >> weight))
                continue;
            items.emplace_back(idx, profit, weight);
        }
        else if (mode == "edges") {
            istringstream iss(line);
            int i, j;
            if (!(iss >> i >> j))
                continue;
            if (i > j) swap(i, j);
            edges.emplace_back(i, j);
        }
    }

    if (n_items < 0 || capacity < 0) {
        err = "Missing param n and/or param c in instance file";
        return false;
    }

    if ((int)items.size() != n_items) {
        err = "Number of items does not match param n";
        return false;
    }

    // validate indices
    vector<int> indices;
    for (auto &t : items) indices.push_back(get<0>(t));
    sort(indices.begin(), indices.end());
    for (int i = 0; i < n_items; ++i) {
        if (indices[i] != i) {
            err = "Item indices must be 0..n-1";
            return false;
        }
    }

    // Create SCIP variables
    vector<SCIP_VAR*> var0(n_items, nullptr);
    vector<SCIP_VAR*> var1(n_items, nullptr);

    int total_weight = 0;
    for (auto &t : items) total_weight += get<2>(t);
    int knapsack_rhs = total_weight - capacity;

    for (int idx = 0; idx < n_items; ++idx) {
        int profit = get<1>(items[idx]);

        char name0[64];
        char name1[64];
        snprintf(name0, sizeof(name0), "x%d_0", idx);
        snprintf(name1, sizeof(name1), "x%d_1", idx);

        SCIP_RETCODE ret;
        SCIP_VAR* v0 = nullptr;
        SCIP_VAR* v1 = nullptr;

        ret = SCIPcreateVarBasic(scip, &v0, name0, 0.0, 1.0, (SCIP_Real)profit, SCIP_VARTYPE_BINARY);
        if (ret != SCIP_OKAY) { err = "SCIP error creating var"; return false; }
        ret = SCIPaddVar(scip, v0);
        if (ret != SCIP_OKAY) { err = "SCIP error adding var"; return false; }

        ret = SCIPcreateVarBasic(scip, &v1, name1, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY);
        if (ret != SCIP_OKAY) { err = "SCIP error creating var"; return false; }
        ret = SCIPaddVar(scip, v1);
        if (ret != SCIP_OKAY) { err = "SCIP error adding var"; return false; }

        var0[idx] = v0;
        var1[idx] = v1;

        SCIPreleaseVar(scip, &v0);
        SCIPreleaseVar(scip, &v1);
    }

    // knapsack constraint: sum weight * x0 >= knapsack_rhs
    {
        vector<SCIP_VAR*> consvars;
        vector<SCIP_Real> coefs;
        consvars.reserve(n_items);
        coefs.reserve(n_items);
        for (int idx = 0; idx < n_items; ++idx) {
            int weight = get<2>(items[idx]);
            consvars.push_back(var0[idx]);
            coefs.push_back((SCIP_Real)weight);
        }
        SCIP_CONS* cons = nullptr;
        SCIP_RETCODE ret = SCIPcreateConsBasicLinear(scip, &cons, "knapsack",
                                                     (int)consvars.size(), consvars.data(), coefs.data(),
                                                     (SCIP_Real)knapsack_rhs, SCIPinfinity(scip));
        if (ret != SCIP_OKAY) { err = "SCIP error creating knapsack constraint"; return false; }
        ret = SCIPaddCons(scip, cons);
        if (ret != SCIP_OKAY) { err = "SCIP error adding knapsack constraint"; return false; }
        {
            SCIP_RETCODE r2 = SCIPreleaseCons(scip, &cons);
            (void)r2;
        }
    }

    // conflict edges: x1_i + x1_j <= 1
    for (size_t k = 0; k < edges.size(); ++k) {
        int i = edges[k].first;
        int j = edges[k].second;
        SCIP_VAR* consvars[2] = { var1[i], var1[j] };
        SCIP_Real coefs[2] = { 1.0, 1.0 };
        char cname[128];
        snprintf(cname, sizeof(cname), "conflict_%zu", k);
        SCIP_CONS* cons = nullptr;
        SCIP_RETCODE ret = SCIPcreateConsBasicLinear(scip, &cons, cname, 2, consvars, coefs, -SCIPinfinity(scip), 1.0);
        if (ret != SCIP_OKAY) { err = "SCIP error creating conflict constraint"; return false; }
        ret = SCIPaddCons(scip, cons);
        if (ret != SCIP_OKAY) { err = "SCIP error adding conflict constraint"; return false; }
        {
            SCIP_RETCODE r2 = SCIPreleaseCons(scip, &cons);
            (void)r2;
        }
    }

    // domain constraints: x0 + x1 = 1
    for (int idx = 0; idx < n_items; ++idx) {
        SCIP_VAR* consvars[2] = { var0[idx], var1[idx] };
        SCIP_Real coefs[2] = { 1.0, 1.0 };
        char cname[128];
        snprintf(cname, sizeof(cname), "domain_constraint_%d", idx);
        SCIP_CONS* cons = nullptr;
        SCIP_RETCODE ret = SCIPcreateConsBasicLinear(scip, &cons, cname, 2, consvars, coefs, 1.0, 1.0);
        if (ret != SCIP_OKAY) { err = "SCIP error creating domain constraint"; return false; }
        ret = SCIPaddCons(scip, cons);
        if (ret != SCIP_OKAY) { err = "SCIP error adding domain constraint"; return false; }
        {
            SCIP_RETCODE r2 = SCIPreleaseCons(scip, &cons);
            (void)r2;
        }
    }

    return true;
}
