---
title: 'cppTPSA/pyTPSA: a C++/python package for truncated power series algebra'
tags:
  - differential algebra
  - truncated power series algebra
  - Python
  - C++
  - accelerator physics
  - astronomy
  - beam dynamics
authors:
  - name: He zhang^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0001-7701-4118
    affiliation: 1 
affiliations:
 - name: Thomas Jefferson National Accelerator Facility, Newport News, VA 23606, USA
   index: 1
date: 15 July 2021
bibliography: paper.bib
---

# Summary

The truncated power series algebra (TPSA), also referred to as differential algebra (DA), is a well established and widely used method in particle accelerator physics and astronomy. The most straightforward usage of TPSA/DA is to calculate the Taylor expansion up to order $n$ of a given function at a specific point, based on which more sophisticated methods have been developed, *e.g.* symplectic tracking, normal form analysis, verified integration, global optimization, *etc*. TPSA/DA can also be used in the fast multipole method for pairwise interactions between particles and in the image reconstruction algorithm.  This package implements the TPSA/DA in C++11 and provides the developers a convenient lib to build the advanced TPSA/DA-based method. A Python 3 lib has also been developed based on the C++ lib and is available in a separate GitHub repository.  

# Background

In the following, we give a very brief introduction from the perspective of computation and practice. Please refer to [] for the theory and more details. 

The fundamental element in DA is the DA vector and first we have to define what a DA vector is. To make the concept easier to understand, we can take a DA vector as the Taylor expansion of a function at a specific point.  

Considering a function $f(\mathbf{x})$ and its Taylor expansion $f_{\mathrm{T}}(\mathbf{x}_0)$ up to the order $n$, we can define a define an equivalence relation between the Taylor expansion and the DA vector as follows

$$f_{\mathrm{T}}(\mathbf{x}_0) = [f]_n = \sum {C_{n_1,n_2, ..., n_v}} \cdot d_1^{n_1} \cdot \dots \cdot d_v^{n_v}, $$ where $\mathbf{x} = (x_1, x_2, \dots, x_v)$, and $n \ge n_1 + n_2 + \dots + n_v$. Here $d_i$ is a special number and it represents a small variance in $x_i$. Generally one can define a DA vector by setting values to respective terms, without define the function $f$. The addition and multiplication of two DA vectors can be defined straightforwardly. To add two DA vectors, we simply add  the coefficients of the like terms in them. To multiply two DA vectors, we calculate the multiplication of each term in the first one with all the terms in the second one and combine like terms ignoring all the terms with an order above $n$. So given two DA vectors $[a]_n$ and $[b]_n$ and a scalar c, we have the following formulas:

\begin{eqnarray}
[a]_{n}+[b]_{n} & := & [a+b]_{n},\nonumber \\
c\cdot[a]_{n} & := & [c\cdot a]_{n},\label{eq:addmlt}\\
{}[a]_{n}\cdot[b]_{n} & := & [a\cdot b]_{n},\nonumber 
\end{eqnarray}

According to the fixed point theorem, the inverse of a DA vector that is not infinitely small can be calculated iteratively in a limit number of iterations. 

The derivation operator $\partial_v$ with respect to the $v^{\mathrm{th}}$ variable can be defined as 

$$ \partial_v[a]_n = \left[ \frac{\partial}{\partial x_v} a \right]_{n-1}, $$

which can be carried out term by term on $[a]_n$. The operator $\partial_v$ satisfies the chain rule:

$$ \partial_v([a]\cdot [b]) = [a]\cdot (\partial_v [b]) + (\partial_v [a])\cdot [b]. $$

The inverse operator $\partial^{-1}_v$ can also be defined and carried out easily in a term-by-term manner. Once the fundamental operators are defined, the DA vector can be used in calculations just as a number. 

# Statement of need

The TPSA/DA methods for particle beam dynamic analysis had been developed in 1980s. The tools are available in a few popular programs for particle accelerator design and simulations, *e.g.* COSY Infinity 9.x, MAD-X, PTC, *etc*. In recent years, the use of TPSA/DA has been extended in other fields, which intrigue the needs for TPSA/DA libraries in popular programming languages. The existing programs are not convenient to developers in other fields. MAD-X is specifically developed for accelerator designs and cannot be used as a general programming languages. Although COSY Infinity can be used as a general programming languages with some limits, it is in lack of some convenient programming features, abundant libraries and a large world-wide community of a modern language as C++ or Python. PTC does include a TPSA/DA library in Fortran 90 but it does not have a user-friendly interface. The TPSA/DA library in C++ is rare. DACE is one alternative. The DACE repository on GitHub had been created but no codes had been released when the author started to develop cppTPSA. Now DACE is available to the public. DACE provides the fundamental DA operations as well as some advanced algorithms based on DA but it has not supported the complex DA vectors, which is useful in normal form analysis. To the best knowledge of the author, there is no other TPSA/DA library in Python 3. 

# Features

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:


# Verification

Figures can be included like this:

# Conclusion

Figures can be included like this:

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.



# References

