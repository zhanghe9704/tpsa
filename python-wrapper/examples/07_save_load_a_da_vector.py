# -*- coding: utf-8 -*-
"""
This example shows how to save a da vector to a file and load a da vector from a file. 
"""

import tpsa
import sys
tpsa.da_init(4, 2, 100) # 4 variables, 2nd order, 100 da vectors.
da = tpsa.base()    # da is the base vector. 
x = tpsa.sin(1+da[0]+0.3*da[1]*da[1])
tpsa.print(x)

# Save the da vector x to a file.
f = open("da_output.txt", "w")
ori_stdout = sys.stdout
sys.stdout = f
tpsa.print(x)
sys.stdout = ori_stdout
f.close()

# Load the da vector from the file.
y = tpsa.assign() # Create an empty da vector y to save the loaded da vector.
tpsa.read_da_from_file("da_output.txt", y)
tpsa.print(y)
