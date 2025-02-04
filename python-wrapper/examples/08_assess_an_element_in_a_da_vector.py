import tpsa
import sys
tpsa.da_init(4, 2, 100) # 4 variables, 2nd order, 100 da vectors.
da = tpsa.base()    # da is the base vector. 
x = tpsa.sin(1+da[0]+0.3*da[1]*da[1])
tpsa.print(x)

l = x.length()
print('The number of elements in x: ', l) # length of the da vector x.

# Obtain the orders of the elements in x and the values of the coefficients.
for i in range(l):
    orders, elem = x.index_element(i)
    print("The orders of the element: ", orders)
    print("The value of the coefficient", elem)

# Alternative way to obtain the values of the coefficients using the orders.
orders = [0,0]
elem = x.element(orders)
print("Coefficient of the constant element: ", elem)

# Alternative way to obtain the values of the coefficients using the index.
elem = x.element(3)
print("Coefficient of the 3rd element: ", elem)
