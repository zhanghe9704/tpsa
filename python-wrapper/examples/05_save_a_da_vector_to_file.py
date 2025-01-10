# -*- coding: utf-8 -*-
"""
This example shows how to save a da vector to a file. 
"""

import tpsa
import sys
tpsa.da_init(4, 1, 100)
da = tpsa.base()
x = tpsa.sin(1+da[0])

f = open("da_output.txt", "w")
ori_stdout = sys.stdout
sys.stdout = f
tpsa.print(x)
sys.stdout = ori_stdout
f.close()
