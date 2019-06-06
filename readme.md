# C++ TPSA Lib

## Announcement

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## About this code

This code allows users to do computations using Truncated Power Series Algebra (TPSA) and/or Differential Algebra (DA). 



For TPSA and DA, please refer to chapter 8 in [*Lecture Notes on Special Topics in Accelerator Physics*](http://inspirehep.net/record/595287/files/slac-pub-9574.pdf)  by  Prof. Alex Chao  and chapter 2 in [*Modern Map Methods in Particle Beam Physics*](http://bt.pa.msu.edu/cgi-bin/display.pl?name=AIEP108book) by Prof. Martin Berz. 



This code is developed based on Dr. Lingyun Yang's tpsa codes in C++ . His codes (tpsa.cpp and tpsa.h) are included in this repository. They are untouched, except for a few functions that are commented off and replaced by functions in tpsa_extend.cc. Please get permission from Dr. Lingyun Yang before you redistribute tpsa.cpp and tpsa.h.



The major change is the memory management. The current memory management works in the following way. Before any TPSA/DA calculation, one needs to reserve the memory for $n$ TPS vector. A memory pool is allocated in the heap, which can be considered as $n$ slots, each for one TPS vector. A linked list is created to reference the  available (unused) slots. Two pointers point to the beginning and the end of the linked list respectively. When a new TPS vector is created, the first slot will be assigned to it and the beginning pointer points to the  following slot in the linked list. When a TPS vector goes out of its scope, the slot will be reset to empty and liked to the end of the list. The end pointer will be updated and point to the new end.  

![Memory management](doc/tpsa_memory_management.png)



A new data type DAVector is created as a wrapper of the TPS vector. The following mathematical operator and functions are overloaded for DAVector, so that a variable of DAVector type can be used as a intrinsic type in calculations. 



Math operator overloaded

| Left hand | Operator | Right hand |
| :-------: | :------: | :--------: |
| DAVector  |    +     |  DAVector  |
|  double   |    +     |  DAVector  |
| DAVector  |    +     |   double   |
|           |    +     |  DAVector  |
| DAVector  |    -     |  DAVector  |
| DAVector  |    -     |   double   |
|  double   |    -     |  DAVector  |
|           |    -     |  DAVector  |
| DAVector  |    *     |  DAVector  |
| DAVector  |    *     |   double   |
|  double   |    *     |  DAVector  |
| DAVector  |    /     |  DAVector  |
| DAVector  |    /     |   double   |
|  double   |    /     |  DAVector  |
| DAVector  |    =     |  DAVector  |
| DAVector  |    =     |   double   |
| DAVector  |    +=    |  DAVector  |
| DAVector  |    +=    |   double   |
| DAVector  |    -=    |  DAVector  |
| DAVector  |    -=    |   double   |
| DAVector  |    *=    |  DAVector  |
| DAVector  |    *=    |   double   |
| DAVector  |    /=    |  DAVector  |
| DAVector  |    /=    |   double   |

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
- erf


Some test results for efficiency are presented in the following. They are done in a Windows 10 desktop with Intel Xeon (R) E5-1620 processor at 3.60 GHz. Table 1 shows the time cost for composition of one DA/TPS vector of six bases with six DA/TPS vectors.  First column shows the order of the vectors, second column the number of terms in each vector, third column the time using the DA data type with revised memory management, and the fourth column the time using the original code. Table 2 shows the time of composition of six DA vectors, each having six bases, with the other group of six DA vectors. The composition in group cost less time if compared with separate compositions. 



Table 1.  Time (in second) of composition 

| Order | No. of terms | DA                    | TPSA                 |
| ----- | ------------ | --------------------- | -------------------- |
| 2     | 28           | $7.57\times 10^{-6}$  | $6.25\times 10^{-6}$ |
| 4     | 210          | $7.50\times 10^{-4}$  | $1.44\times 10^{-2}$ |
| 6     | 924          | $4.48 \times 10^{-3}$ | $8.39\times 10^{-2}$ |
| 8     | 3003         | $9.90 \times 10^{-1}$ | $2.55$               |
| 10    | 8008         | $15.49$               | $44.60$              |



Table 2. Time (in second) of group composition 

| Order | DA                   |
| ----- | -------------------- |
| 2     | $1.51\times 10^{-5}$ |
| 4     | $1.04\times 10^{-3}$ |
| 6     | $4.42\times 10^{-2}$ |
| 8     | $1.05$               |
| 10    | $16.04$              |



More information is available at doc/doxygen/html/index.html.

## How to compile and use this code

You will need a C++ compiler that supports C++ 11 standard. There are three ways to use the code as follows:

* Download the source files. Include "tpsa_extend.h" and "da.h" in your project and compile. 

* The code is developed using Code::Blocks IDE. There are two C::B profiles under the cbp directory: tpsa_lib.cbp and tpsa_dll.cbp for static library and dynamic library respectively. The cbp files are tested in Windows 10 with gcc compiler. 

* You can also use cmake to compile the code into both a static library and a dynamic library. This is only tested in Ubuntu 16.04. 

  `cmake .` 
  `make`

  

  Here is an example of compiling the code under Ubuntu 16.04. 

  Assume I have cloned the codes to the following folder:

  $HOME/tpsa

  

  Inside the above folder, run:

  `cmake .` 

  The Makefile will be generated. 

  

  Then run:

  `make`

  Both the static lib and the shared lib of tpsa will be generated. You can find the following two files:

  libtpsa.a and libtpsaso.so

  

  Now you can use any of them to compile your file. Here let us compile the examples/examples.cc using libtpsa.a. 

  

  `gcc examples/examples.cc -o tpsa_exp -I ./include/ -L ./lib -ltpsa -lstdc++ -lm -std=c++11`

  

  You can also use libtpsaso.so. 

  `gcc examples/examples.cc -o tpsa_exp -std=c++11 -Iinclude -L. -ltpsaso -lstdc++ -lm`

  The executable file tpsa_exp will be generated. 

  

  To run the tpsa_exp file, tell the OS where to find the libtpsaso.so:

  `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/tpsa`

  Run the executable:

  ` ./tpsa_exp `

  

  The result is shown as follows:
  V [16]              Base  [ 4 / 286 ]

  ---------------------------------------------

   1.000000000000000e+00      0  0  0     0  

  -1.000000000000000e+00      1  0  0     1   

    2.000000000000000e+00      0  1  0     2   

    5.000000000000000e-01       0  0  1     3
  
  ......


## Acknowledgement

Thanks to Dr. Lingyun Yang for providing his tpsa code. 

## Contact

Contact the author by hezhang.AT.jlab.org. 

