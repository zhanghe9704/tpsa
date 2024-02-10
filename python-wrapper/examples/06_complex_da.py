import pybind11
import tpsa
from math import *

tpsa.da_init(4,3,1000)
da = tpsa.base()

x1 = da[0] + 2*da[1] + 3*da[2]
x2 = tpsa.sin(x1)
x1 = tpsa.cos(x1)

x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2]
x4 = tpsa.sin(x3)
x3 = tpsa.cos(x3)

y1 = tpsa.complex(x1, x2)
y2 = tpsa.complex(x3, x4)

mmap = tpsa.assign(2)
mmap[0] = x1
mmap[1] = x2

c1 = complex(4.2, 0.3)
c2 = complex(1/3.0, sqrt(2))
c3 = complex(sin(0.7), cos(0.4))

nmap = [c1, c2, c3]
omap = tpsa.da_composition(mmap, nmap)
print('Composition of DA vectors with complex numbers:')
for v in omap:
    print(v)

cnmap = tpsa.assign_cd(3)
cnmap[0] = y1
cnmap[1] = y2
cnmap[2] = y1*y2

comap = tpsa.assign_cd(2)
tpsa.cd_composition(mmap, cnmap, comap)
print('Composition of DA vectors with complex DA vectors:')
for v in comap:
    v.print()

cmmap = tpsa.assign_cd(2)
cmmap[0] = tpsa.complex(x1, tpsa.exp(x1))
cmmap[1] = tpsa.complex(x2, tpsa.exp(x2))

tpsa.cd_composition(cmmap, cnmap, comap)
print('Composition of complex DA vectors with complex DA vectors:')
for v in comap:
    v.print()

mmap.append(x1+0.33*x2)
tpsa.cd_composition(cmmap, mmap, comap)
print('Composition of complex DA vectors with DA vectors:')
for v in comap:
    v.print()

