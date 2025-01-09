# -*- coding: utf-8 -*-
"""
This example shows how to create a DA vector x than only contains one variable 
and calculate sin(x). 
Results are printed to screen.  
"""

import tpsa

tpsa.da_init(4, 1, 100)
da = tpsa.base()
x = 1+da[0]
y = tpsa.sin(x)
x.print()
y.print()
