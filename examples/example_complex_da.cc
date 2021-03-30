#include "../include/da.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <complex>
using namespace std::complex_literals;
using std::complex;

int main() {

    unsigned int da_dim = 3;
    unsigned int da_order = 4;
    unsigned int n_vec = 400;

    //Initialize the DA domain.
    da_init(da_order, da_dim, n_vec);

    DAVector x1, x2, x3, x4;
    x1 = da[0] + 2*da[1] + 3*da[2];
    x2 = sin(x1);
    x1 = cos(x1);

    x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2];
    x4 = sin(x3);
    x3 = cos(x3);

    auto y1 = x1 + x2*1i;
    auto y2 = x3 + x4*1i;

    std::cout<<"y1: "<<std::endl<<y1<<std::endl;
    std::cout<<"y2: "<<std::endl<<y2<<std::endl;

    std::cout<<"y1+y2: "<<std::endl<<y1+y2<<std::endl;
    std::cout<<"y1-y2: "<<std::endl<<y1-y2<<std::endl;
    std::cout<<"y1*y2: "<<std::endl<<y1*y2<<std::endl;
    std::cout<<"y1/y2: "<<std::endl<<y1/y2<<std::endl;

    std::vector<DAVector> mmap;
    mmap.push_back(x1);
    mmap.push_back(x2);

    std::vector<complex<double>> nmap;
    nmap.push_back(4.2+0.3*1i);
    nmap.push_back(1/3.0 + sqrt(2)*1i);
    nmap.push_back(sin(0.7) + cos(0.4)*1i);

    std::vector<complex<double>> omap(2);
    da_composition(mmap, nmap, omap);
    std::cout<<"Composition of DA vectors with complex numbers."<<std::endl<<std::endl;
    for(auto& o: omap) std::cout<<o<<std::endl;

    std::vector<complex<DAVector>> cnmap;
    cnmap.push_back(y1);
    cnmap.push_back(y2);
    cnmap.push_back(y1*y2);

    std::vector<complex<DAVector>> comap(2);
    cd_composition(mmap, cnmap, comap);
    std::cout<<"Composition of DA vectors with complex DA vectors."<<std::endl<<std::endl;
    for(auto& o: comap) std::cout<<o<<std::endl;

    std::vector<complex<DAVector>> cmmap;
    cmmap.push_back(x1+1i*exp(x1));
    cmmap.push_back(x2+1i*exp(x2));

    cd_composition(cmmap, cnmap, comap);
    std::cout<<"Composition of complex DA vectors with complex DA vectors."<<std::endl<<std::endl;
    for(auto& o: comap) std::cout<<o<<std::endl;

    mmap.push_back(x1+0.33*x2);
    cd_composition(cmmap, mmap, comap);
    std::cout<<"Composition of complex DA vectors with DA vectors."<<std::endl<<std::endl;
    for(auto& o: comap) std::cout<<o<<std::endl;

    da_clear();

    return 0;
}

