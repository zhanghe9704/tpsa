//#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "da.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::complex;
using std::string;
using namespace std::complex_literals;
using std::vector;

bool read_da_from_file_c(string filename, DAVector& d) {
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

DAVector devide_by_element_c(DAVector& t, DAVector& b) {
    DAVector r;
    int l = t.length();
    if(l<b.length()) l = b.length();
    // for(int i=0; i<DAVector::full_length(); ++i) {
    for(int i=0; i<l; ++i) {
        double te = t.element(i);
        double be = b.element(i);
        if(std::abs(be)>std::numeric_limits<double>::min()) {
            r.set_element(r.element_orders(i), te/be);
        }
        else {
            r.set_element(r.element_orders(i), te);
        }
    }
    return r;
}

bool compare_da_vectors_c(DAVector& a, DAVector& b, double eps) {
    DAVector r;
    r = a - b;
    r = devide_by_element_c(r, b);
    return r.iszero(eps);
}


TEST_CASE("INITIALIZE DA ENVIRONMENT") {
    int da_dim = 3;
    int da_order = 4;
    int n_vec = 400;
    int init = da_init(da_order, da_dim, n_vec);
    REQUIRE( init  == 0);
}

TEST_CASE("DA FUNCTIONS") {
    DAVector x = 1 + da[0] + 2*da[1] + 5*da[2];
    double eps = 1e-14;

    DAVector y = sqrt(x);
    SECTION("SQRT") {
        REQUIRE(compare_da_with_file("sqrt_da.txt", y, eps));
    }

    y = log(x);
    SECTION("LOG") {
        REQUIRE(compare_da_with_file("log_da.txt", y, eps));
    }

    y = da_int(y,1);
    DAVector m;
    read_da_from_file_c("da_int.txt", m);
    SECTION("DA_INT") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }
    y = log(x);
    y = da_der(y,1);
    read_da_from_file_c("da_der.txt", m);
    SECTION("DA_DER") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }

    y = pow(x,3);
    read_da_from_file_c("pow3_da.txt", m);
    SECTION("POW(x,3)") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }
    
    y = pow(x,0.3);
    SECTION("POW(x,0.3)") {
        REQUIRE(compare_da_with_file("pow0p3_da.txt", y, eps));
    }

    y = exp(x);
    SECTION("EXP") {
        REQUIRE(compare_da_with_file("exp_da.txt", y, eps));
    }

	DAVector z;
    da_substitute(y, 0, 1, z);

    SECTION("SUBSTITUTE A NUMBER") {
        REQUIRE(compare_da_with_file("substitute_number.txt", z, eps));
    }

    da_substitute(y, 0, x, z);

    SECTION("SUBSTITUTE A DA VECTOR") {
        REQUIRE(compare_da_with_file("substitute_da_vector.txt", z, eps));
    }

	std::vector<DAVector> lv(2);
    lv.at(0) = sin(x);
    lv.at(1) = cos(x);
    std::vector<unsigned int> idx{0,1};
    da_substitute(y, idx, lv, z);
    SECTION("SUBSTITUTE MULTIPLE DA VECTORS") {
        REQUIRE(compare_da_with_file("substitute_multiple_da_vectors.txt", z, eps));
    }

	std::vector<DAVector> lx(3);
    std::vector<DAVector> ly(3);
    lx.at(0) = x;
    lx.at(1) = y;
    lx.at(2) = sinh(x);

    da_substitute(lx, idx, lv, ly);

    SECTION("BUNCH SUBSTITUTION") {
        REQUIRE(compare_da_with_file("bunch_substitution_0.txt", ly.at(0), eps));
        REQUIRE(compare_da_with_file("bunch_substitution_1.txt", ly.at(1), eps));
        REQUIRE(compare_da_with_file("bunch_substitution_2.txt", ly.at(2), eps));
    }

	std::vector<DAVector> lu(3);
    lu.at(0) = sin(x);
    lu.at(1) = cos(x);
    lu.at(2) = tan(x);
    da_composition(lx, lu, ly);

    SECTION("DA COMPOSITION") {
        REQUIRE(compare_da_with_file("da_composition_0.txt", ly.at(0), eps));
        REQUIRE(compare_da_with_file("da_composition_1.txt", ly.at(1), eps));
        REQUIRE(compare_da_with_file("da_composition_2.txt", ly.at(2), eps));
    }

}

TEST_CASE("CD FUNCTIONS") {
    double eps = 1e-14;
    DAVector x1, x2, x3, x4;
    x1 = da[0] + 2*da[1] + 3*da[2];
    x2 = sin(x1);
    x1 = cos(x1);

    x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2];
    x4 = sin(x3);
    x3 = cos(x3);

    auto y1 = x1 + x2*1i;
    auto y2 = x3 + x4*1i;

    SECTION("CD FUNDAMENTAL CALCULATIONS") {
        auto r = y1+y2;
        REQUIRE(compare_cd_with_file("cd_calculation_0.txt", r, eps));
        r = y1-y2;
        REQUIRE(compare_cd_with_file("cd_calculation_1.txt", r, eps));
        r = y1*y2;
        REQUIRE(compare_cd_with_file("cd_calculation_2.txt", r, eps));
        r = y1/y2;
        REQUIRE(compare_cd_with_file("cd_calculation_3.txt", r, eps));
    }

    std::vector<DAVector> mmap;
    mmap.push_back(x1);
    mmap.push_back(x2);
    std::vector<complex<DAVector>> cnmap;
    cnmap.push_back(y1);
    cnmap.push_back(y2);
    cnmap.push_back(y1*y2);
    std::vector<complex<DAVector>> comap(2);
    cd_composition(mmap, cnmap, comap);

    SECTION("DA COMPOSITION CD") {
        REQUIRE(compare_cd_with_file("da_composition_cd_0.txt", comap.at(0), eps));
        REQUIRE(compare_cd_with_file("da_composition_cd_1.txt", comap.at(1), eps));
    }

	std::vector<complex<DAVector>> cmmap;
    cmmap.push_back(x1+1i*exp(x1));
    cmmap.push_back(x2+1i*exp(x2));

    cd_composition(cmmap, cnmap, comap);

    SECTION("CD COMPOSITION CD") {
        REQUIRE(compare_cd_with_file("cd_composition_cd_0.txt", comap.at(0), eps));
        REQUIRE(compare_cd_with_file("cd_composition_cd_1.txt", comap.at(1), eps));
    }

	mmap.push_back(x1+0.33*x2);
    cd_composition(cmmap, mmap, comap);

    SECTION("CD COMPOSITION DA") {
        REQUIRE(compare_cd_with_file("cd_composition_da_0.txt", comap.at(0), eps));
        REQUIRE(compare_cd_with_file("cd_composition_da_1.txt", comap.at(1), eps));
    }
}
