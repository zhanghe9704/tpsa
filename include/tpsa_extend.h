/**
 * @file tpsa_extend.h
 * @brief Extend the original tpsa.h file by Dr. Lingyun Yang.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#ifndef TPSA_EXTEND
#define TPSA_EXTEND

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "tpsa.h"

////******************** Additional Methods for TSP Vectors ************************////
unsigned int ad_dim(); //Return the TPS base number.
void ad_reset_vector(const TVEC iv); //Reset all the elements to zero, still keep the length of the vector.
void ad_change_order(unsigned int new_order);//temporarily lower the TPS order.
void ad_restore_order();//Restore the original TPS order, which is larger than the current one
void ad_reserve(const unsigned int n);  //Reserve memory for n TPS vectors.
void ad_clear(); //Destroy the TPS environment and release memory.
void ad_assign(unsigned int &i);    //Assign memory to a TPS vector. The length of the vector is zero.
double ad_con(const TVEC iv); //Return the constant part of the TPS vector
void ad_reset_const(const TVEC iv, double x); //Reset the TPS vector constant element as x and all other elements zero
unsigned int ad_remain();       //Space (number) available for new TPS vectors.
unsigned int ad_count();        //Number of TPS vectors allocated.
void ad_mult_c(const TVEC iv, double c, TVEC ov); //Multiple the TPS vector iv with a constant number c, result stored in ov.
void ad_substitute_const(const TVEC iv, unsigned int base_id, double x, TVEC ov);
void ad_substitute(const TVEC iv, unsigned int base_id, const TVEC v, TVEC ov);
void ad_substitute(const TVEC iv, std::vector<unsigned int> &base_id, std::vector<TVEC> &v, TVEC ov);
void ad_substitute(std::vector<TVEC> &ivecs, std::vector<unsigned int> &base_id, std::vector<TVEC> &v, std::vector<TVEC> &ovecs) ;
void ad_composition(std::vector<TVEC> &ivecs, std::vector<TVEC> &v, std::vector<TVEC> &ovecs);
void ad_composition(std::vector<TVEC> &ivecs, std::vector<double> &v, std::vector<double> &ovecs);
void ad_add(const unsigned int idst, const unsigned int jsrc, unsigned int ov);
void ad_sub(const unsigned int idst, const unsigned int jsrc, TVEC ov);
void ad_add_const(const TVEC i, double r, TVEC ov);
double ad_elem(const TVEC &vec, std::vector<int> &idx);
void ad_elem(const TVEC &vec, unsigned int idx, std::vector<unsigned int>& c, double& x);
void ad_pok(const TVEC &vec, std::vector<int> &idx, double x);
int ad_order(); //Return the order of the TPS environment.
void ad_clean(const TVEC &iv, const double eps); //Set the coefficients smaller than eps to zero.
int ad_full_length();   //Return the maximum length of a TPS.
void ad_int(TVEC iv, unsigned int base_id, TVEC ov); //Integrate w.r.t. the specific base.
int ad_n_element(TVEC v);   //Number of non-zero element!
double ad_norm(TVEC v); //Norm of the TPS vector.
double ad_weighted_norm(TVEC v, double w); //Weighted norm of the TPS vector.
unsigned int ad_last_note();    //Index of the last available slot in the TPS vector pool.
unsigned int ad_next_note();    //Index of the next available slot in the TPS vector pool.

////******************** Expose existing methods in tpsa.cc but not in tpsa.h ************************////
void print_vec(unsigned int ii, std::ostream& os);
#endif // TPSA_EXTEND

