/**
 * @file da.h
 * @brief DA wrapper for the extended tpsa code.
 * @details Defines the DAVector class and the Base class.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#ifndef DA_H_INCLUDED
#define DA_H_INCLUDED

#include <memory>
#include <iostream>
#include <vector>

/** \brief Differential Algebra (DA) Vector
 * A TPS vector with methods. Can be used with most numerical operators.
 *
 */
struct DAVector {
  unsigned int da_vector_;
  DAVector();
  DAVector(const DAVector& da_vector);
  DAVector(DAVector&& da_vector);
  DAVector(double x);
  void print() const;
  double con() const;
  unsigned int length() const;
  int n_element() const;
  void element(unsigned int i, unsigned int *c, double& elem) const;
  void element(unsigned int i, std::vector<unsigned int>& c, double& elem) const;
  double element(int i);
  double element(std::vector<int> idx);
  double norm();
  double weighted_norm(double w);
  void set_element(int *c, double elem);
  void set_element(std::vector<int> idx, double elem);
  void reset();
  void reset_const(double x = 0);
  void clean(const double epx);
  static int dim();
  static int order();
  static int full_length();
  DAVector& operator=(const DAVector& da_vector);
  DAVector& operator=(DAVector&& da_vector);
  DAVector& operator=(double x);
  DAVector& operator=(int x);
  DAVector& operator+=(const DAVector& da_vector);
  DAVector& operator+=(DAVector&& da_vector_);
  DAVector& operator+=(double x);
  DAVector& operator+=(int x);
  DAVector& operator-=(const DAVector& da_vector);
  DAVector& operator-=(DAVector&& da_vector_);
  DAVector& operator-=(double x);
  DAVector& operator-=(int x);
  DAVector& operator*=(const DAVector& da_vector);
  DAVector& operator*=(DAVector&& da_vector);
  DAVector& operator*=(double x);
  DAVector& operator*=(int x);
  DAVector& operator/=(const DAVector& da_vector);
  DAVector& operator/=(DAVector&& da_vector);
  DAVector& operator/=(double x);
  DAVector& operator/=(int x);
  ~DAVector();
};

/** \brief Bases for DA calculations.
 * The i-th base can be accessed as base[i], where base is an object of Base.
 *
 */
struct Base {
  std::vector<DAVector> base;
  Base(const unsigned int n);
  Base(){};
  void set_base(const unsigned int n);
  void set_base();
  const DAVector& operator[](unsigned int i) {return base.at(i);}
};

extern Base da; // Bases for DA calculations. The i-th base can be accessed as da[i].

//Initialize the environment for DA computation. Call this function before any DA computation.
int da_init(unsigned int da_order, unsigned int num_da_variables, unsigned int num_da_vectors);
void da_clear(); //Destroy the DA environment and release memory.
int da_change_order(unsigned int new_order);    //Temporary lower the da order.
int da_restore_order();                         //Restore the original da order.
int da_count();                                 //Number of da variable allocated.
int da_remain();                                //Space (number) available for new da variables.
int da_full_length();                           //Full length of the da vector.
//Take derivative w.r.t. a specific base.
void da_der(const DAVector &da_vector, unsigned int base_id, DAVector &da_vector_der);
//Integrate w.r.t. a specific base.
void da_int(const DAVector &da_vector, unsigned int base_id, DAVector &da_vector_int);
DAVector da_der(const DAVector &da_vector, unsigned int base_id);   //Take derivative w.r.t. a specific base.
DAVector da_int(const DAVector &da_vector, unsigned int base_id);   //Integrate w.r.t. a specific base.
void da_substitute_const(const DAVector &iv, unsigned int base_id, double x, DAVector &ov);
void da_substitute(const DAVector &iv, unsigned int base_id, const DAVector &v, DAVector &ov);
void da_substitute(const DAVector &iv, std::vector<unsigned int> &base_id, std::vector<DAVector> &v, DAVector &ov);
void da_substitute(std::vector<DAVector> &ivecs, std::vector<unsigned int> &base_id, std::vector<DAVector> &v,
                  std::vector<DAVector> &ovecs);
void da_composition(std::vector<DAVector> &ivecs, std::vector<DAVector> &v, std::vector<DAVector> &ovecs);
void da_composition(std::vector<DAVector> &ivecs, std::vector<double> &v, std::vector<double> &ovecs);

DAVector operator+(const DAVector &da_vector, double real_number);
DAVector operator+(double real_number, const DAVector &da_vector) ;
DAVector operator+(const DAVector &da_vector_1, const DAVector &da_vector_2);
DAVector operator*(const DAVector &da_vector, double real_number);
DAVector operator*(double real_number, const DAVector &da_vector);
DAVector operator*(const DAVector &da_vector_1, const DAVector &da_vector_2);
DAVector operator-(const DAVector &da_vector, double real_number);
DAVector operator-(double real_number, const DAVector &da_vector);
DAVector operator-(const DAVector &da_vector_1, const DAVector &da_vector_2);
DAVector operator/(const DAVector &da_vector, double real_number);
DAVector operator/(double real_number, const DAVector &da_vector);
DAVector operator/(const DAVector &da_vector_1, const DAVector &da_vector_2);
DAVector operator+(const DAVector &da_vector);
DAVector operator-(const DAVector &da_vector);

DAVector sqrt(const DAVector &da_vector);
DAVector exp(const DAVector &da_vector);
DAVector log(const DAVector &da_vector);
DAVector sin(const DAVector &da_vector);
DAVector cos(const DAVector &da_vector);
DAVector tan(const DAVector &da_vector);
DAVector asin(const DAVector &da_vector);
DAVector acos(const DAVector &da_vector);
DAVector atan(const DAVector &da_vector);
DAVector sinh(const DAVector &da_vector);
DAVector cosh(const DAVector &da_vector);
DAVector tanh(const DAVector &da_vector);
DAVector pow(const DAVector &da_vector, const int order);
DAVector pow(const DAVector &da_vector, const double order);
double abs(const DAVector &da_vector);
DAVector erf(const DAVector& x);

std::ostream& operator<<(std::ostream &os, const DAVector &da_vector);

void inv_map(std::vector<DAVector> &ivecs, int dim, std::vector<DAVector> &ovecs);
#endif // DA_H_INCLUDED
