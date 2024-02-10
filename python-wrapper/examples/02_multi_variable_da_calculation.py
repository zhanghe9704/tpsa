# -*- coding: utf-8 -*-
"""
This example creates a DA vector that contains three variables 
and calculate sin(x). 
Results are printed to screen. 
"""

import tpsa
tpsa.da_init(4, 3, 100)
da = tpsa.base()
x = 1+da[0]+0.5*da[1]+2*da[2]
y = tpsa.sin(x)
x.print()
y.print()