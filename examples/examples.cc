#include "../include/da.h"
#include <cmath>
#include <vector>

int main() {

    unsigned int da_dim = 3;
    unsigned int da_order = 4;
    unsigned int n_vec = 100;

    da_init(da_order, da_dim, n_vec);

    DAVector x = 1 + da[0] + 2*da[1] + 5*da[2];
    x.print();

    x = exp(x);
    x.print();

    DAVector y = cos(x);
    y.print();

    y += x;
    y.print();

    DAVector x2 = 0.5+da[1];
    y = sin(x)*exp(x2);
    y.print();

    DAVector z = log(y+0.5*da[2]);
    z.print();

    x = 3.14159265;
    x.print();

    DAVector x3 = atan(0.3+da[0]);
    x3.print();

    x = asin(x3);
    x.print();

    x = acos(x3);
    x.print();

    x = sinh(x3);
    x.print();

    x = cosh(x3);
    x.print();

    x = tanh(x3);
    x.print();

    x3.print();
    x = -x3;
    x.print();
    x = +x3;
    x.print();

    x3 = 0.3+da[0]+2*da[1]+3*da[2];
    x = asin(x3);
    x.print();

    x = erf(x3);
    x.print();

   std::vector<DAVector> lx, ly, lz;
   lx.push_back(da[0]);
   lx.push_back(da[1]);
   ly.push_back(1+da[1]);
   ly.push_back(2+da[0]);
   lz.push_back(da[0]);
   lz.push_back(da[1]);
   std::vector<unsigned int> idx{0,1};

   da_subscribe(lx, idx, ly, lz);
   lz.at(0).print();
   lz.at(1).print();


    return 0;
}
