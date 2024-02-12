# -*- coding: utf-8 -*-
"""
Substitute a float for a variable in a DA vector. 
Substitute a DA vector for a variable in a DA vector.
Substitute DA vectors for variables in a DA vector.
Substitute DA vectors for variables in DA vectors. 
"""

import tpsa

tpsa.da_init(4, 3, 400)
da = tpsa.base()

print("1. substitute a float")
x = 1 + da[0] + da[1] + da[2]
x.print()

#Create an empty DA vector y to save the result of subsstitution. 
y = tpsa.assign()
y.print()
#Substitute 1.5 for the 1st variable in x and output to y. 
tpsa.da_substitute_const(x, 0, 1.5, y) 
y.print()

print("2. substitute a DA vector")
z = 1.3 + 0.5*da[1]
#Substitute z for the 3rd variable in x and output to y. 
tpsa.da_substitute(x, 2, z, y)
y.print()


print("3. substitute DA vectors")
idx = [0,1]
#Create a DAVectorList and assigan values into it. 
ly = tpsa.DAVectorList()
ly.append(1.3 + 0.5*da[1])
ly.append(0.4+1.2*da[2])
#Substitute two DA vectors in ly for the 1st and 2nd variables in x respectively 
# and output to y. 
tpsa.da_substitute(x,idx,ly,y)
y.print()


print("4. substitute DA vectors (bunch processing)")
#Another way to assign values into a DAVectorList. 
lx = tpsa.assign(2)
lx[0] = 1 + da[0]
lx[1] = 2 + da[1]
#Create lz for the result
lz = tpsa.assign(2)
tpsa.da_substitute(lx, idx, ly, lz)
lz[0].print()
lz[1].print()




