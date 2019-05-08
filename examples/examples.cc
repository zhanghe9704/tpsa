#include "../include/da.h"

int main() {

    unsigned int da_dim = 3;
    unsigned int da_order = 4;
    unsigned int n_vec = 100;

    da_init(da_order, da_dim, n_vec);

    DAVector x = 1 + da[0];
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

    return 0;
}
