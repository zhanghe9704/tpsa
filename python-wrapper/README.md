# pyTPSA - PYTHON TPSA Lib

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/zhanghe9704/tpsa-python/blob/master/LICENSE.md)

## About this code

This code allows users to do computations using Truncated Power Series Algebra (TPSA) and/or Differential Algebra (DA) in python.

For TPSA and DA, please refer to chapter 8 in [*Lecture Notes on Special Topics in Accelerator Physics*](http://inspirehep.net/record/595287/files/slac-pub-9574.pdf)  by  Prof. Alex Chao  and chapter 2 in [*Modern Map Methods in Particle Beam Physics*](http://bt.pa.msu.edu/cgi-bin/display.pl?name=AIEP108book) by Prof. Martin Berz. 

## How to compile and install

To compile the codes, you will need a C++ compiler that supports C++ 11 standard, make, python3, pybind11 and scikit-build-core.  If you do not have the above programs in your system, you can install them using the following commands. 

```
sudo apt install build-essential 
sudo apt install python python-dev python-pip  
```

(Note: pybind11 works with python2.7, too. So in principle you can compile the code for python 2.7 and use the lib in python 2.7. However, I have not tried it. )

Assuming you are in the "tpsa" folder, go to the "python-wrapper" folder to compile and install the lib using PIP. 

```shell
cd python-wrapper
pip install .
```

**WARNING**: do NOT run "pip install python-wrapper" under the "tpsa" folder because cmake will find the wrong CMakeList.txt.  

## Use the lib in python3

### DA environment initialization and simple DA computation

In side python, load the lib by

`>>> import tpsa`

Initialize the DA environment by `

```
>>> tpsa.da_init(3, 2, 1000)
0
```

The above command initialize a DA environment that holds at most 1000 DA vectors of two variables up to order three. Then let us define the bases as da:

`>>> da = tpsa.base()`

Each individual base can be accessed by da[i], where i is the index of the base. In the above case, i ranges in [0,1]. For example,`

```
>>> da[3]
<tpsa.DAVector object at 0x7f1ca31921f0>
```

After we defined the bases, we can define other DA vectors and perform computations using the bases. 

```
>>> x = 0.5 + da[0] + 2.1*da[1]
>>> y = tpsa.sin(x)
>>> y.print()
         V [27]              Base  [ 10 / 10 ]
----------------------------------------------
  4.794255386042030e-01     0 0     0
  8.775825618903728e-01     1 0     1
  1.842923379969783e+00     0 1     2
 -2.397127693021015e-01     2 0     3
 -1.006793631068826e+00     1 1     4
 -1.057133312622268e+00     0 2     5
 -1.462637603150621e-01     3 0     6
 -9.214616899848915e-01     2 1     7
 -1.935069548968272e+00     1 2     8
 -1.354548684277791e+00     0 3     9
```

The list of mathematical operators and functions supported can be found in the following. 

### Substitute for variables in DA vectors

Let us  see a 1D example first. Assume we want to calculate the Taylor expansion of sin(x) at x=1 up to the forth order. We can do it as follows. 

```
>>> import tpsa
>>> tpsa.da_init(4,1,100)
0
>>> da = tpsa.base()
>>> x = tpsa.sin(1+da[0])
>>> x.print()
 I           V [16]              Base  [ 5 / 5 ]
------------------------------------------------
 1   8.414709848078965e-01     0     0
 2   5.403023058681398e-01     1     1
 3  -4.207354924039483e-01     2     2
 4  -9.005038431135663e-02     3     3
 5   3.506129103366235e-02     4     4
```

Now if we want to evaluate sin(1.2) using the above extension, we need to substitute 0.2 for the variable (base) in the DA vector x. 

```
>>> y = tpsa.da_substitute_const(x, 0, 0.2)
>>> y.print()
 I           V [20]              Base  [ 1 / 5 ]
------------------------------------------------
 1   9.320377212765296e-01     0     0
```

In the first line, we substitute 0.2 for the 1st variable in x and store the result in y. Generally, y is a DA vector, considering x can have more than one bases. But for this specific case, y is just a real number. 

Besides substitute a real number, we can substitute a DA vector, too. 

```
>>> v = 0.2 + da[0]
>>> y = tpsa.da_substitute(x, 0, v)
>>> y.print()
 I           V [26]              Base  [ 5 / 5 ]
------------------------------------------------
 1   9.320377212765296e-01     0     0
 2   3.623240241022748e-01     1     1
 3  -4.663510131426833e-01     2     2
 4  -6.200135148442674e-02     3     3
 5   3.506129103366235e-02     4     4
```

If we have a multiple-dimensional DA vector, we can substitute DA vectors for more than one bases at once. 

```
>>> import tpsa
>>> tpsa.da_init(4, 3, 100)
0
>>> da = tpsa.base()
>>> x = 1.0 + da[0] + 2*da[1] + 0.2*da[2]
>>> v = tpsa.assign(2)
>>> v[0] = 1 + 0.5*da[0]+ 3*da[1] + 2*da[2]
>>> v[1] = 2 + da[1] + da[2]
>>> idx = [0,1]
>>> y = tpsa.da_substitute(x, idx, v)
>>> y.print()
 I           V [32]              Base  [ 4 / 35 ]
------------------------------------------------
 1   6.000000000000000e+00     0 0 0     0
 2   5.000000000000000e-01     1 0 0     1
 3   5.000000000000000e+00     0 1 0     2
 4   4.200000000000000e+00     0 0 1     3
```

In the above code, we substitute DA vectors stored in "v" for the bases determined by "idx". The size of "idx" should be equal to the size of "v".  

Please note that here we use the following  command to create a DAVectorList:

```
>>> v = tpsa.assign(2)
```

The function assign() without argument create a DA vector and assign(n) with a integer n create a DAVectorList that contains n DA vectors with all coefficients zero. 

 The substitution can be carried out in bunches, which means "x" and "y" can be a DAVectorLists that store multiple DA vectors.The substitution will be performed on all the DA vectors in "x" and the results are saved into "y".  An example is shown as follows. 

```
>>> import tpsa
>>> tpsa.da_init(4, 3, 100)
0
>>> da = tpsa.base()
>>> x = tpsa.assign(2)
>>> x[0] = 1.0 + da[0] + 2*da[1] + 0.2*da[2]
>>> x[1] = 0.5 + da[2]
>>> v = tpsa.assign(2)
>>> v[0] = 1 + 0.5*da[0]+ 3*da[1] + 2*da[2]
>>> v[1] = 2 + da[1] + da[2]
>>> idx = [0,1]
>>> y = tpsa.assign(2)
>>> tpsa.da_substitute(x, idx, v, y)
>>> y[0].print()
 I           V [60]              Base  [ 4 / 35 ]
------------------------------------------------
 1   6.000000000000000e+00     0 0 0     0
 2   5.000000000000000e-01     1 0 0     1
 3   5.000000000000000e+00     0 1 0     2
 4   4.200000000000000e+00     0 0 1     3

>>> y[1].print()
 I           V [61]              Base  [ 4 / 35 ]
------------------------------------------------
 1   5.000000000000000e-01     0 0 0     0
 2   1.000000000000000e+00     0 0 1     3
```

### Composition of DA vectors

In the n dimensional DA domain, if we want to substitute n DA vectors for all the n  bases in m DA vectors, we can use the function da_composition. 

```
>>> import tpsa
>>> tpsa.da_init(4, 3, 100)
0
>>> da = tpsa.base()
>>> x = tpsa.assign(2)
>>> x[0] = 1.0 + da[0] + 2*da[1] + 0.2*da[2]
>>> x[1] = 0.5 + da[2]
>>> v = tpsa.assign(3)
>>> v[0] = 1 + 0.5*da[0]+ 3*da[1] + 2*da[2]
>>> v[1] = 0.5 + da[2]
>>> v[2] = da[0] + 3.3*da[1]
>>> y = tpsa.assign(2)
>>> tpsa.da_composition(x, v, y)
>>> y[0].print()
 I           V [60]              Base  [ 4 / 35 ]
------------------------------------------------
 1   3.000000000000000e+00     0 0 0     0
 2   7.000000000000000e-01     1 0 0     1
 3   3.660000000000000e+00     0 1 0     2
 4   4.000000000000000e+00     0 0 1     3

>>> y[1].print()
 I           V [61]              Base  [ 3 / 35 ]
------------------------------------------------
 1   5.000000000000000e-01     0 0 0     0
 2   1.000000000000000e+00     1 0 0     1
 3   3.300000000000000e+00     0 1 0     2
```

Besides DA vectors, we can also substitute numbers for the bases and the results are numbers too. 

```
>>> v = [1.0, 2.2, 3]
>>> y = tpsa.da_composition(x, v)
>>> print(y)
[7.0, 3.5]
```

In both the examples above, the size of the second arguments has to be equal to the dimension of the DA domain. The size of the first arguments can vary and the size of the result is always equal to  the size of the first arguments. 

For more examples of using this lib, please check out the files in the **examples** folder.

### Operators and functions overload

Currently, the tpsa lib supports the following operators and math functions. 

- Math operator overloaded: (DA - DA vector, CD - complex DA vector)

  | Left hand | Operator | Right hand |
  | :-------: | :------: | :--------: |
  |   DA/CD   |    +     |   DA/CD    |
  |  double   |    +     |   DA/CD    |
  |   DA/CD   |    +     |   double   |
  |           |    +     |   DA/CD    |
  |   DA/CD   |    -     |   DA/CD    |
  |   DA/CD   |    -     |   double   |
  |  double   |    -     |   DA/CD    |
  |           |    -     |   DA/CD    |
  |   DA/CD   |    *     |   DA/CD    |
  |   DA/CD   |    *     |   double   |
  |  double   |    *     |   DA/CD    |
  |   DA/CD   |    /     |   DA/CD    |
  |   DA/CD   |    /     |   double   |
  |  double   |    /     |   DA/CD    |
  |   DA/CD   |    =     |   DA/CD    |
  |   DA/CD   |    =     |   double   |
  |   DA/CD   |    +=    |   DA/CD    |
  |   DA/CD   |    +=    |   double   |
  |   DA/CD   |    -=    |   DA/CD    |
  |   DA/CD   |    -=    |   double   |
  |   DA/CD   |    *=    |   DA/CD    |
  |   DA/CD   |    *=    |   double   |
  |   DA/CD   |    /=    |   DA/CD    |
  |   DA/CD   |    /=    |   double   |

  Math functions overloaded:

  - sqrt
  - exp
  - log
  - sin
  - cos
  - tan
  - asin
  - acos
  - atan
  - sinh
  - cosh
  - tanh
  - pow
  - abs
  - erf (DA only)

## Known issues

1. When some temporary variables in the C++ lib go out of scope, the memory of them are not released immediately while they do in pure C++ environment, although eventually they will be released a few steps after the function call finishes. This means we may need a larger DA vector pool in Python than in C++. Usually a pool size of a few hundred to a few thousand should be large enough, which is not a problem for a modern personal computer.   

2. We observed "segmentation fault" error in the following two scenarios: (1) When running a TPSA script in the command line using 

   ```shell
   python my_tpsa_scrpt.py
   ```

   a segmentation fault may happen after the script is carried out. (2) When using the TPSA lib in a Python interactive environment, a segmentation fault may happen after the following command is called:

   ```shell
   exit()
   ```

    to quit the Python environment. Unfortunately, I do not know what exactly happen in the above scenarios. we suspect the segmentation error is caused by destructors in C++ to reset the memory slots of the exiting DA vectors. However, we observed the segmentation fault can be avoided by enlarge the memory pool when initializing the TPSA environment using the following command:

   ```shell
   tpsa.da_init(DA_Order, Variable_Number, DA_Vector_Number)
   ```

   Putting a large number for "DA_Vector_Number" tends to prevent the segmentation error from happening. 

## Acknowledgement

Thanks to Dr. Lingyun Yang for providing his tpsa code. 

Thanks to pybind11 developers.

## Contact

Contact the author by hezhang.AT.jlab.org. 
