#include "../include/da.h"
#include "../include/tpsa_extend.h"

typedef DAVector DAV;
int main() {

    unsigned int da_dim = 3;
    unsigned int da_order = 10;
    unsigned int n_vec = 100;

    da_init(da_order, da_dim, n_vec);

    DAVector x = 1 - da[0] +  2*da[1] + 0.5*da[2];
    x.print();

    da_clear();

    return 0;
}
