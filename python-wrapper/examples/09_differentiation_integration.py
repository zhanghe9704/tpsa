import tpsa
import sys
tpsa.da_init(4, 2, 100) # 4 variables, 2nd order, 100 da vectors.
da = tpsa.base()    # da is the base vector. 
x = tpsa.sin(1+da[0]+0.3*da[1]*da[1])
tpsa.print(x)

# Take the derivative of x with respect to the 1st variable.
y = tpsa.da_der(x, 0)
tpsa.print(y)

#
z = tpsa.da_int(x, 0)
tpsa.print(z)