# TPSA Lib

## Announcement

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## About this code

This code allows users to do computations using Truncated Power Series Algebra (TPSA) and/or Differential Algebra (DA). 



For TPSA and DA, please refer to chapter 8 in [*Lecture Notes on Special Topics in Accelerator Physics*](http://inspirehep.net/record/595287/files/slac-pub-9574.pdf)  by  Prof. Alex Chao  and chapter 2 in [*Modern Map Methods in Particle Beam Physics*](http://bt.pa.msu.edu/cgi-bin/display.pl?name=AIEP108book) by Prof. Martin Berz. 



This code is developed based on Dr. Lingyun Yang's tpsa codes in C++ . His codes (tpsa.cpp and tpsa.h) are included in this repository. They are untouched, except for a few functions that are commented off and replaced in tpsa_extend.cc. Please get permission from Dr. Lingyun Yang before you redistribute tpsa.cpp and tpsa.h.



## How to compile and use this code



## Math operators and functions

The following math operators have been overloaded for the DA Vectors.

| Left hand | Operator | Right hand |
| :-------: | :------: | :--------: |
| DAVector  |    +     |  DAVector  |
|  double   |    +     |  DAVector  |
| DAVector  |    +     |   double   |
| DAVector  |    -     |  DAVector  |
| DAVector  |    -     |   double   |
|  double   |    -     |  DAVector  |
| DAVector  |    *     |  DAVector  |
| DAVector  |    *     |   double   |
|  double   |    *     |  DAVector  |
| DAVector  |    /     |  DAVector  |
| DAVector  |    /     |   double   |
|  double   |    /     |  DAVector  |

The following math functions support the DA vectors. 

- sqrt
- exp
- log
- sin
- cos
- pow
- abs

## Acknowledgement

Thanks to Dr. Lingyun Yang for providing his tpsa code. 

## Contact

Contact the author by hezhang.AT.jlab.org. 

