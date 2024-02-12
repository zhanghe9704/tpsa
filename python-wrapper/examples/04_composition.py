# -*- coding: utf-8 -*-
"""
Composition of two function g(x) and f(x) as g(f(x)). 
g(x) should be DA vector(s). f(x) could be float(s) or DA vector(s).
"""

import tpsa

tpsa.da_init(4,2,400)
da = tpsa.base()

print("g(f(x)) with f(x) floats: ")
x=tpsa.assign(2)
x[0] = 1 + da[0] + 2*da[1]
x[1] = 0.5 + 3*da[0] + da[1]

y = [1.0, 2.0]

z = tpsa.da_composition(x, y)
print(z)

print("g(f(x)) with f(x) DA vectors: ")
y=tpsa.assign(2)
y[0] = 1 + 2*da[0] + da[1]
y[1] = 2 + da[0] + 0.5*da[1]

z=tpsa.assign(2)

tpsa.da_composition(x, y, z)
z[0].print()
z[1].print()