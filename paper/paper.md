---
title: 'cppTPSA/pyTPSA: a C++/Python package for truncated power series algebra'
tags:
  - differential algebra
  - truncated power series algebra
  - Python
  - C++
  - accelerator physics
  - astronomy
  - beam dynamics
authors:
  - name: He Zhang^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0001-7701-4118
    affiliation: 1 
affiliations:
 - name: Thomas Jefferson National Accelerator Facility, Newport News, VA 23606, USA
   index: 1
date: 15 July 2021
bibliography: paper.bib
---

# Summary

The truncated power series algebra (TPSA), also referred to as differential algebra (DA), is a well-established and widely used method in particle accelerator physics and astronomy. The most straightforward usage of TPSA/DA is to calculate the Taylor expansion  of a given function at a specific point up to order $n$. In recent years, as the application of TPSA/TA has been extended to other fields, a reusable implementation of TPSA/DA as a modern C++ library or other high level programming language like Python has become desirable. The cppTPSA package implements TPSA/DA in C++11 and provides developers a convenient library with which to build advanced TPSA/DA-based methods. A Python 3 library, pyTPSA, has also been developed based on the C++ lib.  

# Background

In the following, we give a very brief introduction on TPSA/DA from a practical computational perspective. Please refer to @AIEP108book and @chao2002lecture for the complete theory with more details. 

The fundamental concept in DA is the DA vector. To make this concept easier to understand, we can consider a DA vector as the Taylor expansion of a function at a specific point.  

Considering a function $f(\mathbf{x})$ and its Taylor expansion $f_{\mathrm{T}}(\mathbf{x}_0)$  at the point $\mathbf{x}_0$ up to the order $n$, we can define  an equivalence relation between the Taylor expansion and the DA vector as follows

$$ [f]_n = f_{\mathrm{T}}(\mathbf{x}_0) = \sum {C_{n_1,n_2, ..., n_v}} \cdot d_1^{n_1} \cdot \dots \cdot d_v^{n_v}, $$ where $\mathbf{x} = (x_1, x_2, \dots, x_v)$, and $n \ge n_1 + n_2 + \dots + n_v$. Here $d_i$ is a special number: it represents a small variance in $x_i$. Generally one can define a DA vector by directly setting values to respective terms, without defining the function $f$. The addition and multiplication of two DA vectors can be defined straightforwardly. To add two DA vectors, we simply add  the coefficients of the like terms. To multiply two DA vectors, we multiply each term in the first one with all the terms in the second one and combine like terms while ignoring all terms above order $n$. So given two DA vectors $[a]_n$ and $[b]_n$ and a scalar c, we have the following formulae:

\begin{eqnarray}
[a]_{n}+[b]_{n} & := & [a+b]_{n},\nonumber \\
c\cdot[a]_{n} & := & [c\cdot a]_{n},\label{eq:addmlt}\\
{}[a]_{n}\cdot[b]_{n} & := & [a\cdot b]_{n},\nonumber 
\end{eqnarray}

According to the fixed point theorem [@AIEP108book], the inverse of a DA vector that is not infinitely small can be calculated iteratively in a limited number of iterations. 

The derivation operator $\partial_v$ with respect to the $v^{\mathrm{th}}$ variable can be defined as 

$$ \partial_v[a]_n = \left[ \frac{\partial}{\partial x_v} a \right]_{n-1}, $$

which can be carried out term by term on $[a]_n$. The operator $\partial_v$ satisfies the chain rule:

$$ \partial_v([a]\cdot [b]) = [a]\cdot (\partial_v [b]) + (\partial_v [a])\cdot [b]. $$

The inverse operator $\partial^{-1}_v$ can also be defined and carried out easily in a term-by-term manner. Once the fundamental operators are defined, the DA vector can be used in calculations just as a number. More sophisticated methods using DA have been developed, *e.g.* symplectic tracking [@caprimap], normal form analysis [@monthnf], verified integration [@rdaint], global optimization [@makino2005verified,], fast multipole method for pairwise interactions between particles [@FMMCPO2010].  

# Statement of need

TPSA/DA methods for particle beam dynamic analysis were developed in the 1980s. These tools are available in several popular programs for particle accelerator design and simulations, such as COSY Infinity 9 [@COSYCAP04], MAD-X [@grote2003mad; @MADX], and PTC [@forest2002introduction]. In recent years, the application of TPSA/DA has been extended to other fields, motivating the development of TPSA/DA libraries in popular programming languages. However, the existing programs are not convenient for developers from other fields. For instance, MAD-X is specifically developed for accelerator design and cannot be used as a general programming language. Although COSY Infinity can be used as a general programming languages, it lacks some of the convenient programming features found in modern languages, such as C++ or Python, along with abundant libraries and  a large supporting community. PTC does include a TPSA/DA library in Fortran 90 but it lacks a user-friendly interface. TPSA/DA libraries in C++ are rare. DACE [@massari2018differential; @DACE] is one alternative. The DACE repository on GitHub had been created but no codes had been released when the author began developing cppTPSA[@cppTPSA]. Now DACE is available to the public,  providing fundamental DA operations and some advanced algorithms based on DA. However, it does not support complex DA vectors, which are useful in normal form analysis. To the author's best knowledge, there is no other TPSA/DA library in Python 3. 

# Features

This library consists of a C++ library that performs TPSA/DA calculations and a Python wrapper. Users can compile the source code into a static or shared library or generate a Python library for a Python 3 environment.  The readme file in the repository provides detailed instructions on how to compile  both the C++ library and the Python library, respectively. 



The C++ library is based on Lingyun Yang's TPSA code [@yang9array], which is also incorporated into MAD-X [@MADX]. During development, we tried to make minimal changes from the original code, but had to revise or rewrite some functions for better efficiency and/or consistency.  One big change is the memory management. In Yang’s code, the pointers to all the DA vectors are stored in a vector. Whenever a new DA vector is required, the program searches this vector for the first empty pointer and allocates the memory. Once a DA vector is out of scope, the memory is freed. In contrast, our library initiates a memory pool for all DA vectors (with the number defined by the user) at the very start, during the initialization of the DA environment. The addresses for the slots, each for one DA vector, in the pool are maintained in a linked list. Whenever we need to create a new DA vector, we take out a slot from the beginning of the list. Whenever a DA vector goes out of the scope, its destructor will set all value in the slot to zero and put it back at the end of the list. The memory pool is managed simply by manipulating the two pointers: one pointing to the start and the other to the end of the list. This method eliminates the repetitive searching and allocation/deallocation operations, thereby achieving better efficiency.



Some new features have been added, as follows: 

1. Add a DA vector data type and define the commonly used math operators for it, so that users can use a DA vector as simple as a normal number in calculations. 
2. Support the complex DA vector defined by the C++ complex template. 
3. More math functions are supported. (A list of the overloaded math functions can be found in the readme file of the repository.) 
4. Add new functions that perform the composition of (complex) DA vectors, which can carry out multiple compositions in a call. 
5. A Python wrapper is provided. 



The following C++ code shows an example of a simple TPSA/DA calculation. After initializing an environment that can contain at most 400 three dimensional DA vectors up to the 4th order, two DA vectors x1 and x2 and a complex DA vector y1 are defined, some trigonometric functions are performed on them, and the results are output to the screen. 

```c++
    #include "da.h"
    da_init(4, 3, 400);
    DAVector x1, x2;
    x1 = da[0] + 2*da[1] + 3*da[2];
    x2 = sin(x1);
    x1 = cos(x1);
    auto y1 = x1 + x2*1i;
    std::cout<<x1<<x2<<std::endl;
    std::cout<<sin(y1)<<std::endl;
```



A Python example doing the same calculation is presented as follows. 

```python
    import tpsa
    tpsa.da_init(4, 3, 400)
    da = tpsa.base()
    x1 = da[0] + 2*da[1] + 3*da[2]
    x2 = tpsa.sin(x1)
    x1 = tpsa.cos(x1)
    y1 = tpsa.complex(x1, x2)
    print(x1)
    print(x2)
    print(tpsa.sin(y1))
```

More examples can be found in the repository. 


# Verification

This library has been verified with COSY Infinity 9.0. As an example, the outputs of calculating sin (0.3+da[0]+2×da[1]) up to the fourth order by both programs are presented in \autoref{fig:cosy} and  \autoref{fig:cpptpsa} respectively. \autoref{fig:cosy} shows the result generated by COSY Infinity, while \autoref{fig:cpptpsa} shows the result generated by cppTPSA. The two programs give exactly the same result. In most cases, the two programs agree to the machine's precision. However, one may observe difference in the coefficients at orders of $10^{-15}$ or $10^{-16}$ for some special functions such as arcsin. This is because different numerical recipes are used in calculation. For example, a special function may be approximated by different series.  This small deviation is usually considered acceptable in practice. If higher precision is desired, one could/should consider the Taylor Model (TM) datatype in COSY Infinity. The TM vector calculates a DA vector together with its error band. However, it is outside the scope of this code. Please note the sequence of the terms may be different when outputting a DA vector from cppTPSA and from COSY Infinity. 



![COSY Infinity 9.0 output.\label{fig:cosy}](cosy-output.png)

![cppTPSA output.\label{fig:cpptpsa}](cppTPSA-output.png)

# Acknowledgements

The author would like to thank Dr. Lingyun Yang for providing his source code. 

This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Nuclear Physics under contract DE-AC05-06OR23177.



# References

