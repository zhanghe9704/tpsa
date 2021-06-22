//#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../include/da.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>

using std::complex;
using std::string;
using namespace std::complex_literals;

const std::string k_whitespace = " \n\r\t\v\f"; //\n newline, \r carriage return, \t horizontal tab, \v vertical tab, \f form feed.

// Trim the spaces at the head and the tail of the string
string ltrim_whitespace(string input_line) {
    string::size_type st = input_line.find_first_not_of(k_whitespace);
    return (st == string::npos)? "" : input_line.substr(st);
}

string rtrim_whitespace(string input_line) {
    string::size_type fi = input_line.find_last_not_of(k_whitespace);
    return (fi == string::npos)? "" : input_line.substr(0, fi+1);
}

std::string trim_whitespace(std::string input_line) {
    return rtrim_whitespace(ltrim_whitespace(input_line));
}

int read_da_from_file(string filename, DAVector& d) {
    d.reset();
    std::fstream input;
    input.open(filename, std::ios::in | std::ios::out | std::ios::app);
    string line;

    bool reading = false;
    bool read_success = false;
    int n_terms = 0;
    while(std::getline(input, line)) {
        if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
        if (!line.empty()) {
            line = trim_whitespace(line);
            if(reading) {
                int i = std::stoi(line.substr(line.length()-6));
                line.erase(line.end()-6, line.end());
                double elem = 0;
                std::vector<int> idx;

                std::istringstream iss(line);
                string word;
                iss >> word;
                iss >> word;
                elem = std::stod(word);
                while(!iss.eof()) {
                    int index = 0;
                    iss >> word;
                    index = std::stoi(word);
                    idx.push_back(index);
                }
                d.set_element(idx, elem);
                if(n_terms-i == 1) read_success = true;
            }
            else {
                auto pi = line.find_last_of("/");
                auto pf = line.find_last_of("[");
                if(pi!=string::npos && pf!=string::npos) {
                    n_terms = std::stoi(line.substr(pi+1,pf));
                    continue;
                }
                line.erase(std::remove(line.begin(), line.end(), '-'), line.end());
                if(line.empty()) {
                    reading = true;
                    continue;
                }
            }
        }
    }
    return read_success;
}

int read_cd_from_file(string filename, complex<DAVector>& cd) {
    DAVector& rl = get_real(cd);
    DAVector& img = get_imag(cd);
    rl.reset();
    img.reset();
    std::fstream input;
    input.open(filename, std::ios::in | std::ios::out | std::ios::app);
    string line;

    bool reading = false;
    bool read_success = false;
    int n_terms = 0;
    while(std::getline(input, line)) {
        if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
        if (!line.empty()) {
            line = trim_whitespace(line);
            if(reading) {
                int i = std::stoi(line.substr(line.length()-6));
                line.erase(line.end()-6, line.end());
                double elem_rl = 0;
                double elem_img = 0;
                std::vector<int> idx;

                std::istringstream iss(line);
                string word;
                iss >> word;
                iss >> word;
                elem_rl = std::stod(word);
                iss >> word;
                elem_img = std::stod(word);
                while(!iss.eof()) {
                    int index = 0;
                    iss >> word;
                    index = std::stoi(word);
                    idx.push_back(index);
                }
                if(std::abs(elem_rl)>=std::numeric_limits<double>::min()) rl.set_element(idx, elem_rl);
                if(std::abs(elem_img)>=std::numeric_limits<double>::min()) img.set_element(idx, elem_img);
                if(n_terms-i == 1) read_success = true;
            }
            else {
                auto pi = line.find_last_of("/");
                auto pf = line.find_last_of("[");
                if(pi!=string::npos && pf!=string::npos) {
                    n_terms = std::stoi(line.substr(pi+1,pf));
                    continue;
                }
                line.erase(std::remove(line.begin(), line.end(), '-'), line.end());
                if(line.empty()) {
                    reading = true;
                    continue;
                }
            }
        }
    }
    return read_success;
}

TEST_CASE("DA_INIT") {
    int da_dim = 3;
    int da_order = 4;
    int n_vec = 400;
    int init = da_init(da_order, da_dim, n_vec);
    REQUIRE( init  == 0);
}

TEST_CASE("TEMP") {
    DAVector x = sin(3+da[0]);
    DAVector y = 2+da[1];
    x.print();
    std::cout<<y+x*1i<<std::endl;
    string filename = "da_output.txt";
    read_da_from_file(filename, x);
    x.print();

    complex<DAVector> zz;
    read_cd_from_file("cd_output.txt", zz);
    std::cout<<zz<<std::endl;
    REQUIRE( 0  == 0);
}
