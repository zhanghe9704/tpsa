#include "../include/da.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

int main() {

    unsigned int da_dim = 3;
    unsigned int da_order = 4;
    unsigned int n_vec = 400;

    //Initialize the DA domain.
    da_init(da_order, da_dim, n_vec);

    std::cout<<"Print out the base vector."<<std::endl<<std::endl;
    da[0].print();
    da[1].print();
    da[2].print();

    std::cout<<"Fundamental calculations of DA vectors."<<std::endl<<std::endl;
    DAVector x = 1 + da[0] + 2*da[1] + 5*da[2];
    x.print();

    DAVector y = exp(x);
    y.print();

    //Substitute a number for a base.
    std::cout<<"Substitute a number for a base."<<std::endl<<std::endl;
    DAVector z;
    da_substitute(y, 0, 1, z);
    z.print();

    //Substitute a DA vector for a base.
    std::cout<<"Substitute a DA vector for a base."<<std::endl<<std::endl;
    da_substitute(y, 0, x, z);
    z.print();

    //Substitute multiple DA vectors for bases at once.
    std::cout<<"Substitute multiple DA vectors for bases at once."<<std::endl<<std::endl;
    std::vector<DAVector> lv(2);
    lv.at(0) = sin(x);
    lv.at(1) = cos(x);
    std::vector<unsigned int> idx{0,1};
    da_substitute(y, idx, lv, z);
    z.print();

    std::cout<<"The norm of z is "<<z.norm()<<std::endl;
    std::cout<<"The weighted norm of z is "<<z.weighted_norm(0.1)<<std::endl;

    //Bunch processing for substitutions.
    std::cout<<"Bunch processing for substitutions."<<std::endl<<std::endl;
    std::vector<DAVector> lx(3);
    std::vector<DAVector> ly(3);
    lx.at(0) = x;
    lx.at(1) = y;
    lx.at(2) = sinh(x);

    da_substitute(lx, idx, lv, ly);
    ly.at(0).print();
    ly.at(1).print();
    ly.at(2).print();

    //Composition of DA vectors with numbers.
    std::cout<<"Composition of DA vectors with numbers."<<std::endl<<std::endl;
    std::vector<double> lm{0.1, 2, 1};
    std::vector<double> ln(3);
    da_composition(lx, lm, ln);
    for(auto& x: ln) std::cout<<x<<' ';
    std::cout<<std::endl<<std::endl;

    //Composition of DA vectors with DA vectors.
    std::cout<<"Composition of DA vectors with DA vectors."<<std::endl<<std::endl;
    std::vector<DAVector> lu(3);
    lu.at(0) = sin(x);
    lu.at(1) = cos(x);
    lu.at(2) = tan(x);
    da_composition(lx, lu, ly);
    ly.at(0).print();
    ly.at(1).print();
    ly.at(2).print();

	//Output a DA vector to file.
	std::ofstream fout;
	fout.open("da_output.txt");
	fout << z;
	fout.close();

	z /= 1e5;
	std::cout<<"z:"<<std::endl;
	z.print();
	std::cout<<"norm of z: "<<z.norm()<<std::endl;

	std::cout<<"DA eps: "<<DAVector::eps()<<std::endl;
	DAVector::set_eps(1e-20);
	std::cout<<"Reset DA eps: "<<DAVector::eps()<<std::endl;

    return 0;
}
