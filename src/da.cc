/**
 * @file da.cc
 * @brief DA wrapper for the extended tpsa code.
 * @details Defines methods for DAVectors.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#include "../include/da.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <set>
#include "../include/tpsa_extend.h"

using std::complex;
using std::vector;

Base da; /**< Bases for DA calculations. The i-th base can be accessed as da[i].  */

double DAVector::eps = 1e-16;/**<Global threshold of the abs value of DAVector coefficients. */

DAVector::DAVector() {ad_alloc(&da_vector_);}

//DAVector::DAVector(bool b) {
//    if(b) ad_alloc(&da_vector_);
//    else ad_assign(da_vector_);
//}

DAVector::DAVector(const DAVector& da_vector) {
      ad_alloc(&da_vector_);
      ad_copy(&da_vector.da_vector_, &da_vector_);
}

DAVector::DAVector(DAVector&& da_vector) {
    ad_assign(da_vector_);
    std::swap(da_vector_, da_vector.da_vector_);
    ad_free(&da_vector.da_vector_);
}

DAVector::DAVector(double x) {
    ad_alloc(&da_vector_);
    ad_reset_const(da_vector_, x);
}

DAVector::DAVector(int i) {
    ad_alloc(&da_vector_);
    ad_reset_const(da_vector_, i);
}

DAVector::~DAVector() {
    ad_free(&da_vector_);
}

void DAVector::clear() {
    ad_free(&da_vector_);
}

void DAVector::print() const { ad_print(&da_vector_);}      /**< Print out a DA vector. */
double DAVector::con() const {return ad_con(da_vector_);}   /**< Return the constant element of a DA vector. */
int DAVector::n_element() const {return ad_n_element(da_vector_);} /**< Return the number of non-zero element a DA vector. */

DAVector& DAVector::operator=(DAVector&& da_vector) {
    if (this==&da_vector) return *this;
    std::swap(this->da_vector_, da_vector.da_vector_);
    ad_reset(&da_vector.da_vector_);
    return *this;
}

DAVector& DAVector::operator=(const DAVector &da_vector) { ad_copy(&da_vector.da_vector_, &da_vector_); return *this;}

DAVector& DAVector::operator=(double x) {
    this->reset_const(x);
    return *this;
}

DAVector& DAVector::operator=(int x) {
    this->reset_const(static_cast<double>(x));
    return *this;
}

DAVector& DAVector::operator+=(const DAVector& da_vector) {ad_add(&da_vector_, &da_vector.da_vector_); return *this;}

DAVector& DAVector::operator+=(DAVector&& da_vector) {*this = *this + da_vector; return *this;}

DAVector& DAVector::operator+=(double x) {ad_add_const(&da_vector_, &x); return *this;}

DAVector& DAVector::operator+=(int x) {double xx = static_cast<double>(x); ad_add_const(&da_vector_, &xx); return *this;}

DAVector& DAVector::operator-=(const DAVector& da_vector) {ad_sub(&da_vector_, &da_vector.da_vector_); return *this;}

DAVector& DAVector::operator-=(DAVector&& da_vector) {*this = *this - da_vector; return *this;}

DAVector& DAVector::operator-=(double x) {x*=-1; ad_add_const(&da_vector_, &x); return *this;}

DAVector& DAVector::operator-=(int x) {double xx = static_cast<double>(-x); ad_add_const(&da_vector_, &xx); return *this;}

DAVector& DAVector::operator*=(const DAVector& da_vector) {*this = *this * da_vector; return *this;}

DAVector& DAVector::operator*=(DAVector&& da_vector) {*this = *this * da_vector; return *this;}

DAVector& DAVector::operator*=(double x){ad_mult_const(&da_vector_, &x); return *this;}

DAVector& DAVector::operator*=(int x){double xx = static_cast<double>(x); ad_mult_const(&da_vector_, &xx); return *this;}

DAVector& DAVector::operator/=(const DAVector& da_vector) {*this = *this / da_vector; return *this;}

DAVector& DAVector::operator/=(DAVector&& da_vector) {*this = *this / da_vector; return *this;}

DAVector& DAVector::operator/=(double x) {ad_div_c(&da_vector_, &x); return *this;}

DAVector& DAVector::operator/=(int x){double xx = static_cast<double>(x); ad_div_c(&da_vector_, &xx); return *this;}

void DAVector::reset() { ad_reset_vector(da_vector_);}   /**< Reset all element to zero, vector length unchanged. */
/// Reset the value to the given number. Vector length is set to 1.
void DAVector::reset_const(double x) {ad_reset_const(da_vector_, x);}
int DAVector::dim() { return ad_dim();}          /**< Return the DA base number. */
int DAVector::order() {return ad_order();}

/// Return the length of the DA vector.
unsigned int DAVector::length() const {
  	unsigned int length;
  	ad_length(&da_vector_, &length);
  	return length;
}

int DAVector::full_length() {
    return ad_full_length();
}

/** \brief Return the value and the order pattern of the specific element.
 * Given the ordinal number of an element, return the value and order pattern of the element. Following the c++ tradition,
 * the ordinal number, i, starts from zero. (In ad_elem, the ordinal number starts from one.) The size of vector c should be
 * equal to the number of bases. For example, if i matches the element (x^nx)*(n^ny)*(z^nz), c = {nx, ny, nz}, where
 * x, y, and z are the bases.
 * \param[in] i The ordinal number of the element.
 * \param[out] c The order pattern of the element.
 * \param[out] elem The value of the element.
 * \return
 *
 */
void DAVector::element(unsigned int i, std::vector<unsigned int>& c, double& elem) const {
    unsigned int ii = i+1;
    if (i<0) ii = 0;
  	ad_elem(da_vector_, ii, c, elem);
}

/** \brief Return the value and the order pattern of the specific element.
 * Given the ordinal number of an element, return the value and order pattern of the element. Following the c++ tradition,
 * the ordinal number, i, starts from zero. (In ad_elem, the ordinal number starts from one.) The size of array c should be
 * equal to the number of bases. For example, if i matches the element (x^nx)*(n^ny)*(z^nz), c = {nx, ny, nz}, where
 * x, y, and z are the bases.
 * \param[in] i The ordinal number of the element.
 * \param[out] c The order pattern of the element.
 * \param[out] elem The value of the element.
 * \return
 *
 */
void DAVector::element(unsigned int i, unsigned int *c, double& elem) const {
    unsigned int ii = i+1;
    if (i<0) ii = 0;
  	ad_elem(&da_vector_, &ii, c, &elem);
}

/** \brief Return the value of a specific element.
 * Given the ordinal number of the element, return the value. The ordinal number starts from zero, following c++ tradition.
 * \param i The ordinal number of the element.
 * \return The value of the element.
 *
 */
double DAVector::element(int i) {
    unsigned int ii = i+1;
    if (i<0) ii = 0;
    double elem = 0;
    int d = dim();
    unsigned int *c = new unsigned int[d];
    ad_elem(&da_vector_, &ii, c, &elem);
    delete[] c;
    return elem;
}

/** \brief Return the value of a specific element.
 * Given the order pattern of the element, return the value. The size of the vector idx should be equal to the base number.
 * For example, assuming the base number is three, idx = {nx, ny, nz} referring to the element (x^nx)*(y^ny)*(z^nz),
 * where x, y, and z are the bases.
 * \param idx The order pattern of the element.
 * \return The value of the element.
 *
 */
double DAVector::element(std::vector<int>idx) {
    return ad_elem(da_vector_, idx);
}

std::vector<int>& DAVector::element_orders(int i) {
    return ad_element_orders(i);
}

/** \brief Check if all the abs value of the coefficients in the DAVector are smaller than the cut-off value.
 * \return True or false.
 *
 */
bool DAVector::iszero() const {
    return ad_zero_check(da_vector_, eps);
}

/** \brief Check if all the abs value of the coefficients in the DAVector are smaller than the given value eps.
 * \param eps The value to compare with the abs value of the coefficients.
 * \return True or false.
 *
 */
bool DAVector::iszero(double eps) const{
    return ad_zero_check(da_vector_, eps);
}

/** \brief Return the norm of the DA vector, e.g. the maximum of the absolute value of the DA coefficients.
 * \return The norm of the DA vector.
 *
 */
double DAVector::norm() {
    return ad_norm(da_vector_);
}

/** \brief Return the weighted norm of the DA vector.
 * Calculate the absolute value of coef*w^n for each term, where coef is the coefficient of the term, n is the total
 * order of the term and w is the weight. Return the maximum of the calculation.
 *  \param w The weight.
 *  \return The norm of the DA vector.
 */
double DAVector::weighted_norm(double w) {
    return ad_weighted_norm(da_vector_, w);
}

/** \brief Set the value of a specific element.
 *
 * \param c The order pattern of the element.
 * \param elem The value of the element.
 * \return void.
 *
 */
void DAVector::set_element(int *c, double elem) {
    size_t n = DAVector::dim();
    ad_pok(&da_vector_, c, &n, &elem);
}

/** \brief Set the value of a specific element.
 *
 * \param idx The order pattern of the element.
 * \param elem The value of the element.
 * \return void.
 *
 */
void DAVector::set_element(std::vector<int> idx, double elem) {
    ad_pok(da_vector_, idx, elem);
}

/** \brief Set a coefficient in the DA Vector to be zero if the abs of the
 * coefficient is less than the given eps.
 * \param[in] the threshold for zero.
 * \return void.
 */
void DAVector::clean(const double eps) {
    ad_clean(da_vector_, eps);
}

/** \brief Set a coefficient in the DA Vector to be zero if the abs of the
 * coefficient is less than the global eps.
 * \return void.
 */
void DAVector::clean() {
    ad_clean(da_vector_, DAVector::eps);
}

Base::Base(const unsigned int n) {
    for(unsigned int i = 0; i<n; ++i) {
        DAVector v;
        double zero = 0.0;
		ad_var(&v.da_vector_, &zero, &i);
		base.push_back(v);
    }
}

/** \brief Set the DA base vectors.
 *
 * \param n Number of the bases.
 * \return void.
 *
 */
void Base::set_base(const unsigned int n) {
    base.resize(n);
    for(unsigned int i = 0; i<n; ++i) {
        double zero = 0.0;
		ad_var(&base.at(i).da_vector_, &zero, &i);
    }
}

/// Set the DA base vectors.
void Base::set_base() {
    unsigned int n = DAVector::dim();
    set_base(n);
}

/** \brief Set the cut-off value for DA coefficients.If the abs value of an coefficient is less than eps, it will be
 * set to zero.
 * \param eps The cut-off value for DA coefficients. Should be greater than zero.
 * \return void
 *
 */
void da_set_eps(double eps) {
    if(eps>0) DAVector::eps = eps;
    else std::cout<<"Warning: negative value of eps is ignored in da_set_eps!"<<std::endl;
}

/** \brief Take derivative of a da vector w.r.t. a specific base
 *
 * \param da_vector The input da vector.
 * \param base_id Id of the base.
 * \param da_vector_der The output da vector. Should be different with the input vector.
 * \return void
 *
 */
void da_der(const DAVector &da_vector, unsigned int base_id, DAVector &da_vector_der) {
    assert(base_id>=0&&base_id<DAVector::dim()&&"Base out of limits in DA_DER!");
    ad_derivative(&da_vector.da_vector_, &base_id, &da_vector_der.da_vector_);
}

/** \brief Integrate a da vector w.r.t. a specific base
 *
 * \param da_vector The input da vector.
 * \param base_id Id of the base.
 * \param da_vector_int The output da vector. Should be different with the input vector.
 * \return void
 *
 */
void da_int(const DAVector &da_vector, unsigned int base_id, DAVector &da_vector_int) {
    assert(base_id>=0&&base_id<DAVector::dim()&&"Base out of limits in da_int!");
    ad_int(da_vector.da_vector_, base_id, da_vector_int.da_vector_);
}

/** \brief Take derivative of a da vector w.r.t. a specific base
 *
 * \param da_vector The input da vector.
 * \param base_id Id of the base.
 * \return da vector as the derivation result.
 *
 */
DAVector da_der(const DAVector &da_vector, unsigned int base_id) {
    DAVector res;
    assert(base_id>=0&&base_id<DAVector::dim()&&"Base out of limits in da_der!");
    ad_derivative(&da_vector.da_vector_, &base_id, &res.da_vector_);
    return res;
}

/** \brief Integrate a da vector w.r.t. a specific base
 *
 * \param da_vector The input da vector.
 * \param base_id Id of the base.
 * \return da vector as the integration result.
 *
 */
DAVector da_int(const DAVector &da_vector, unsigned int base_id) {
    DAVector res;
    assert(base_id>=0&&base_id<DAVector::dim()&&"Base out of limits in da_int!");
    ad_int(da_vector.da_vector_, base_id, res.da_vector_);
    return res;
}

/// Return the full length of the da vector.
int da_full_length() {
    return ad_full_length();
}

///Return the orders of each base as a vector for the i-th element.
std::vector<int>& da_element_orders(int i) {
    return ad_element_orders(i);
}

/// Count the number of existing DA Vectors, including the bases.
int da_count() {
    return ad_count();
}

/** \brief Check the remaining capacity in the DA pool.
 *
 * \return The maximum number of DA vectors that can be created in the remaining place of the pool.
 *
 */
int da_remain() {
    return ad_remain();
}
//Initialize the DA environment.
/** \brief Initialize the DA environment.
 * This function initialize the DA environment w.r.t. the given DA order and the number of bases.
 * A memory pool is created for the given number of DA vectors, which should be large enough for perspective calculations.
 * Base vectors are also created. The i-th base can be accessed as da[i].
 * This function should be called before any DA calculation.
 * \param da_order Highest order of a DA Vector.
 * \param num_da_variables Number of bases.
 * \return num_da_vectors Size of the DA memory pool.
 *
 */
int da_init(unsigned int da_order, unsigned int num_da_variables, unsigned int num_da_vectors, bool table) {
    assert(da_order>0&&num_da_variables>0&&num_da_vectors>0&&"Wrong parameters when initialize DA environment!");
	ad_init(&num_da_variables, &da_order);
	ad_reserve(&num_da_vectors);
	da.set_base(num_da_variables);
	if(table) ad_generate_order_table();
	return 0;
}

/** \brief Destroy the DA environment and release memory.
 * This function should only be called when all the DAVectors are out of scope. Otherwise it will cause a segment fault.
 * \return void
 *
 */

void da_clear() {
    da.base.clear();
    ad_clear();
}

///Temporarily change the DA order.
int da_change_order(unsigned int new_order) {ad_change_order(new_order); return 0;}

///Restore the original DA order.
int da_restore_order(){ad_restore_order(); return 0;}

/** \brief Substitute the given value to the specific base in a DA vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal number of the base.
 * \param[in] x The alue to substitute.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_substitute_const(const DAVector &iv, unsigned int base_id, double x, DAVector &ov) {
    ad_substitute_const(iv.da_vector_, base_id, x, ov.da_vector_);
}

/** \brief Substitute a DA vector to a specific base in a given DA Vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal number of the base.
 * \param[in] v The DA vector to substitute.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_substitute(const DAVector &iv, unsigned int base_id, const DAVector &v, DAVector &ov) {
    ad_substitute(iv.da_vector_, base_id, v.da_vector_, ov.da_vector_);
}

/** \brief Substitute DA vectors to the specific bases in a DA vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal numbers of the bases.
 * \param[in] v The DA vectors to substitute. The size of v should be equal to the size of base_id.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_substitute(const DAVector &iv, std::vector<unsigned int> &base_id, std::vector<DAVector> &v, DAVector &ov) {
    std::set<unsigned int> check;
    for(auto& id: base_id) check.insert(id);
    assert(check.size()==base_id.size()&&"Duplicate base id in da_subscribe!");
    std::vector<TVEC> ad_v;
    for(auto&& i : v) ad_v.push_back(i.da_vector_);
    ad_substitute(iv.da_vector_, base_id, ad_v, ov.da_vector_);
}

/** \brief Substitute DA vectors to the specific bases in a group of DA vectors.
 *
 * \param ivecs A group of DA Vectors.
 * \param base_id The ordinal numbers for the bases.
 * \param v The DA vectors to substitute. The size of v should be equal to the size of base_id.
 * \param ovecs The DA vectors that saves the result. The size of ovecs should be equal to the size of ivecs.
 * \return
 *
 */
void da_substitute(std::vector<DAVector> &ivecs, std::vector<unsigned int> &base_id, std::vector<DAVector> &v,
                  std::vector<DAVector> &ovecs) {
    std::set<unsigned int> check;
    for(auto& id: base_id) check.insert(id);
    assert(check.size()==base_id.size()&&"Duplicate base id in da_subscribe!");
    std::vector<TVEC> ad_iv, ad_v, ad_ov;
    for(auto&& i : v) ad_v.push_back(i.da_vector_);
    for(auto&& i : ivecs) ad_iv.push_back(i.da_vector_);
    for(auto&& i : ovecs) ad_ov.push_back(i.da_vector_);
    ad_substitute(ad_iv, base_id, ad_v, ad_ov);
}

/** \brief Composition of a group of DA vectors with another group of DA vectors.
 * Calculate the composition of two groups of DA vectors, and save the result into the third group of DA vectors.
 * Call the first group of DA vectors, ivecs, as f, which includes n DA vectors: f_1, f_2, ..., f_n.
 * The second group of DA vectors, v, should contain the same number of DA vectors as the number of DA bases.
 * Call the third group of DA vectors, ovecs, as g, which also includes n DA vectors: g_1, g_2, ..., g_n.
 * g = f(v), or g_1 = f_1(v), g_2 = f_2(v), ..., g_n = f_n(v).
 * \param[in] ivecs DA vectors.
 * \param[in] v DA vectors.
 * \param[out] ovecs Result DA vectors.
 * \return void.
 *
 */
void da_composition(std::vector<DAVector> &ivecs, std::vector<DAVector> &v, std::vector<DAVector> &ovecs) {
    std::vector<TVEC> ad_iv, ad_v, ad_ov;
    for(auto&& i : v) ad_v.push_back(i.da_vector_);
    for(auto&& i : ivecs) ad_iv.push_back(i.da_vector_);
    for(auto&& i : ovecs) ad_ov.push_back(i.da_vector_);
    ad_composition(ad_iv, ad_v, ad_ov);
}

/** \brief Submitting all the bases for a group of DA vectors with given values.
 *
 * \param[in] ivecs A group of DA Vectors.
 * \param[in] v The values for all the bases.
 * \param[out] ovecs Results.
 * \return void.
 *
 */
void da_composition(std::vector<DAVector> &ivecs, std::vector<double> &v, std::vector<double> &ovecs) {
    std::vector<TVEC> ad_iv;
    for(auto&& i : ivecs) ad_iv.push_back(i.da_vector_);
    ad_composition(ad_iv, v, ovecs);
}

/** \brief Submitting all the bases for a group of DA vectors with given complex values.
 *
 * \param[in] ivecs A group of DA Vectors.
 * \param[in] v The complex values for all the bases.
 * \param[out] ovecs Results.
 * \return void.
 *
 */
void da_composition(std::vector<DAVector> &ivecs, std::vector<std::complex<double>> &v,
                    std::vector<std::complex<double>> &ovecs) {
    std::vector<TVEC> ad_iv;
    for(auto&& i : ivecs) ad_iv.push_back(i.da_vector_);
    ad_composition(ad_iv, v, ovecs);
}

complex<DAVector>& complex_da_pow_int_pos(const complex<DAVector>& iv, std::vector<complex<DAVector>> &power_v,
                                          const int order, int order_rec) {
    if (abs(power_v.at(order*order_rec))>0) {
        return power_v.at(order*order_rec);
    }
    else if (order==0) {
        power_v.at(0) = 1;
        return power_v.at(0);
    }
    else if (order==1) {
        double norm = abs(power_v.at(1*order_rec));
        if(norm<std::numeric_limits<double>::min()) {
            power_v.at(1*order_rec) = iv;
        }
        return power_v.at(1*order_rec);
    }
    else if (order==2) {
        double norm  = abs(power_v.at(2*order_rec));
        if(norm<std::numeric_limits<double>::min()) {
            power_v.at(2*order_rec) = iv*iv;
        }
        return power_v.at(2*order_rec);
    }
    else if (order&1) { //Odd order
        unsigned int order_idx = order_rec*2;

        double norm = abs(power_v.at(order_idx));
        if (norm<std::numeric_limits<double>::min()) power_v.at(order_idx) = iv*iv;
        complex<DAVector>& vres = power_v.at(order_idx);

        vres = complex_da_pow_int_pos(vres, power_v, order/2, order_idx);
        power_v.at(order_rec+(order/2)*order_idx) = vres * iv;
        return power_v.at(order_rec+(order/2)*order_idx);
    }
    else {              //Even order
        unsigned int order_idx = order_rec*2;
        double norm = abs(power_v.at(order_idx));
        if (norm<std::numeric_limits<double>::min()) power_v.at(order_idx) = iv * iv;
        complex<DAVector>& vres = power_v.at(order_idx);
        return complex_da_pow_int_pos(vres, power_v, order/2, order_idx);
    }
}

void da_composition(std::vector<DAVector> &ivecs, std::vector<std::complex<DAVector>> &v,
                    std::vector<std::complex<DAVector>> &ovecs) {
    int gnv = DAVector::dim();
    int gnd = DAVector::order();
    assert(gnv==v.size()&&"Error in da_composition: No. of DA vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size()==ovecs.size()&&"Error in da_composition: No. of input vectors NOT EQUAL to No. of output vectors!");
    vector<vector<complex<DAVector>>>  power_vv(gnv, vector<complex<DAVector>>(gnd+1));
    for(int i=0; i<gnv; ++i) {
        power_vv.at(i).at(0) = 1;
        power_vv.at(i).at(1) = v.at(i);
    }

    int veclen_max = 0;
    for(auto& v: ivecs) {
        if(v.length()>veclen_max) veclen_max = v.length();
    }

    //Copy the constant part of the input vectors to the output vectors.
    int vec_size = ivecs.size();
    for(int iv=0; iv<vec_size; ++iv) {
        ovecs.at(iv) = ivecs.at(iv).con();
    }
    if(!ad_valid_order_table()) ad_generate_order_table();
    for(int i=1; i<veclen_max; ++i) {//loop over all elements, except for the constant element.
        auto orders = da_element_orders(i);
        bool product_flag = true;  //Calculate the product.
        complex<DAVector> product(1,0);
        bool c_flag = true;
        for(int iv=0; iv<vec_size; ++iv) {
            if(i>= ivecs.at(iv).length()) continue;
            if (std::abs(ivecs.at(iv).element(i))<std::numeric_limits<double>::min()) continue;
            if (c_flag) {
                for(int id=0; id<gnv; ++id) {
                    if (orders.at(id)>0) {
                         complex_da_pow_int_pos(v.at(id), power_vv.at(id), orders.at(id), 1);
                    }
                }
                c_flag = false;
            }
            double coef = ivecs.at(iv).element(i);
            if (product_flag) {
                for(int id=0; id<gnv; ++id) {
                    int order = orders.at(id);
                    if(order>0) product *= power_vv.at(id).at(order);
                }
                product_flag = false;
            }
            ovecs.at(iv) += coef*product;
        }
    }
    for(auto& v: ovecs){
        v.real().clean();
        v.imag().clean();
    }
}

void da_composition(std::vector<std::complex<DAVector>> &ivecs, std::vector<std::complex<DAVector>> &v,
                    std::vector<std::complex<DAVector>> &ovecs) {
    assert(DAVector::dim()==v.size()&&"Error in da_composition: No. of DA vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size()==ovecs.size()&&"Error in da_composition: No. of input vectors NOT EQUAL to No. of output vectors!");
    int n = ivecs.size();
    std::vector<DAVector> iv(2*n);
    for(int i=0; i<n; ++i) {
        iv.at(2*i) = ivecs.at(i).real();
        iv.at(2*i+1) = ivecs.at(i).imag();
    }
    std::vector<complex<DAVector>> ov(2*n);
    da_composition(iv, v, ov);
    complex<double> i1(0,1);
    for(int i=0; i<n; ++i) {
        ovecs.at(i) = ov.at(2*i) + i1*ov.at(2*i+1);
    }
}

//Overload the operators for DA
DAVector operator+(const DAVector &da_vector, double real_number) {
	DAVector res;
	ad_add_const(da_vector.da_vector_, real_number, res.da_vector_);
	return res;
}

DAVector operator+(double real_number, const DAVector &da_vector) {
    DAVector res;
	ad_add_const(da_vector.da_vector_, real_number, res.da_vector_);
	return res;
}

complex<DAVector> operator+(const DAVector &da_vector, std::complex<double> complex_number) {
    complex<DAVector> res(complex_number.real()+da_vector, complex_number.imag());
    return res;
}

complex<DAVector> operator+( std::complex<double> complex_number, const DAVector &da_vector) {
    complex<DAVector> res(complex_number.real()+da_vector, complex_number.imag());
    return res;
}

complex<DAVector> operator-(const DAVector &da_vector, std::complex<double> complex_number) {
    complex<DAVector> res(da_vector-complex_number.real(), complex_number.imag());
    return res;
}

complex<DAVector> operator-( std::complex<double> complex_number, const DAVector &da_vector) {
    complex<DAVector> res(complex_number.real()-da_vector, complex_number.imag());
    return res;
}

complex<DAVector> operator*(const DAVector &da_vector, std::complex<double> complex_number) {
    complex<DAVector> res(da_vector*complex_number.real(), da_vector*complex_number.imag());
    return res;
}

complex<DAVector> operator*( std::complex<double> complex_number, const DAVector &da_vector) {
    complex<DAVector> res(complex_number.real()*da_vector, complex_number.imag()*da_vector);
    return res;
}

complex<DAVector> operator/(const DAVector &da_vector, std::complex<double> complex_number) {
    complex<double> c = 1.0/complex_number;
    complex<DAVector> res(da_vector*c.real(), da_vector*c.imag());
    return res;
}

complex<DAVector> operator/( std::complex<double> complex_number, const DAVector &da_vector) {
    DAVector d = 1.0/da_vector;
    complex<DAVector> res(complex_number.real()*d, complex_number.imag()*d);
    return res;
}

std::complex<DAVector>  operator+(const std::complex<DAVector> &cd_vector, double number) {
    complex<DAVector> res(cd_vector.real()+number, cd_vector.imag());
    return res;
}
std::complex<DAVector>  operator+(double number, const std::complex<DAVector> &cd_vector) {
    complex<DAVector> res(cd_vector.real()+number, cd_vector.imag());
    return res;
}
std::complex<DAVector>  operator-(const std::complex<DAVector> &cd_vector, double number) {
    complex<DAVector> res(cd_vector.real()-number, cd_vector.imag());
    return res;
}
std::complex<DAVector>  operator-(double number, const std::complex<DAVector> &cd_vector) {
    complex<DAVector> res(number-cd_vector.real(), cd_vector.imag());
    return res;
}
std::complex<DAVector>  operator*(const std::complex<DAVector> &cd_vector, double number) {
    DAVector rv = cd_vector.real()*number;
    DAVector iv = cd_vector.imag()*number;
    complex<DAVector> res(rv, iv);
    return res;
}
std::complex<DAVector>  operator*(double number, const std::complex<DAVector> &cd_vector) {
    DAVector rv = cd_vector.real()*number;
    DAVector iv = cd_vector.imag()*number;
    complex<DAVector> res(rv, iv);
    return res;
}
std::complex<DAVector>  operator/(const std::complex<DAVector> &cd_vector, double number){
    DAVector rv = cd_vector.real()/number;
    DAVector iv = cd_vector.imag()/number;
    complex<DAVector> res(rv, iv);
    return res;
}
std::complex<DAVector>  operator/(double number, const std::complex<DAVector> &cd_vector){
    DAVector rv = cd_vector.real();
    DAVector iv = cd_vector.imag();
    DAVector v = 1/(rv*rv+iv*iv);
    complex<DAVector> res(number*v*rv, -number*v*iv);
    return res;
}

std::complex<DAVector>  operator+(const std::complex<DAVector> &cd_vector, std::complex<double> complex_number) {
    complex<DAVector> res(cd_vector.real()+complex_number.real(), cd_vector.imag()+complex_number.imag());
    return res;
}
std::complex<DAVector>  operator+(std::complex<double> complex_number, const std::complex<DAVector> &cd_vector) {
    complex<DAVector> res(cd_vector.real()+complex_number.real(), cd_vector.imag()+complex_number.imag());
    return res;
}
std::complex<DAVector>  operator-(const std::complex<DAVector> &cd_vector, std::complex<double> complex_number){
    complex<DAVector> res(cd_vector.real()-complex_number.real(), cd_vector.imag()-complex_number.imag());
    return res;
}
std::complex<DAVector>  operator-(std::complex<double> complex_number, const std::complex<DAVector> &cd_vector){
    complex<DAVector> res(complex_number.real()-cd_vector.real(), complex_number.imag()-cd_vector.imag());
    return res;
}
std::complex<DAVector>  operator*(const std::complex<DAVector> &cd_vector, std::complex<double> complex_number) {
    complex<DAVector> vr = cd_vector*complex_number.real();
    complex<DAVector> vi = cd_vector*complex_number.imag();
    complex<DAVector> res(vr.real()-vi.imag(), vr.imag()+vi.real());
    return res;
}
std::complex<DAVector>  operator*(std::complex<double> complex_number, const std::complex<DAVector> &cd_vector){
    complex<DAVector> vr = cd_vector*complex_number.real();
    complex<DAVector> vi = cd_vector*complex_number.imag();
    complex<DAVector> res(vr.real()-vi.imag(), vr.imag()+vi.real());
    return res;
}
std::complex<DAVector>  operator/(const std::complex<DAVector> &cd_vector, std::complex<double> complex_number) {
    complex<double> c = 1.0/complex_number;
    return c*cd_vector;
}
std::complex<DAVector>  operator/(std::complex<double> complex_number, const std::complex<DAVector> &cd_vector) {
    complex<DAVector> v = 1.0/cd_vector;
    return complex_number*v;
}


bool operator==(const DAVector &da_vector_1, const DAVector &da_vector_2) {
    DAVector c = da_vector_1 - da_vector_2;
    return c.iszero();
}

DAVector operator+(const DAVector &da_vector_1, const DAVector &da_vector_2) {
	DAVector res;
	ad_add(da_vector_1.da_vector_, da_vector_2.da_vector_, res.da_vector_);
	return res;
}

DAVector operator*(const DAVector &da_vector, double real_number) {
	DAVector res;
	ad_mult_c(da_vector.da_vector_, real_number, res.da_vector_);
	return res;
}

DAVector operator*(double real_number, const DAVector &da_vector) {
	DAVector res;
	ad_mult_c(da_vector.da_vector_, real_number, res.da_vector_);
	return res;
}

DAVector operator*(const DAVector &da_vector_1, const DAVector &da_vector_2) {
    DAVector res;
    ad_mult(&da_vector_1.da_vector_, &da_vector_2.da_vector_, &res.da_vector_);
	return res;
}

DAVector operator-(const DAVector &da_vector, double real_number) {
	DAVector res;
	ad_add_const(da_vector.da_vector_, -1*real_number, res.da_vector_);
	return res;
}

DAVector operator-(double real_number, const DAVector &da_vector) {
	DAVector res;
	ad_mult_c(da_vector.da_vector_, -1, res.da_vector_);
	ad_add_const(&res.da_vector_, &real_number);
	return res;
}

DAVector operator-(const DAVector &da_vector_1, const DAVector &da_vector_2) {
	DAVector res;
	ad_sub(da_vector_1.da_vector_, da_vector_2.da_vector_, res.da_vector_);
	return res;
}

DAVector operator/(const DAVector &da_vector, double real_number) {
    if (std::abs(real_number) < std::numeric_limits<double>::min()) {
        std::cerr << "ERROR: divide a too small number! " << real_number << std::endl;
        std::exit(-1);
    }
	DAVector res;
	ad_mult_c(da_vector.da_vector_, 1.0/real_number, res.da_vector_);
	return res;
}

DAVector operator/(double real_number, const DAVector &da_vector) {
    DAVector res;
    ad_c_div(&da_vector.da_vector_, &real_number, &res.da_vector_);
	return res;
}

DAVector operator/(const DAVector &da_vector_1, const DAVector &da_vector_2) {
    DAVector res;
    ad_div(&da_vector_1.da_vector_, &da_vector_2.da_vector_, &res.da_vector_);
	return res;
}

DAVector operator-(const DAVector &da_vector) {
    DAVector res;
    ad_mult_c(da_vector.da_vector_, -1, res.da_vector_);
    return res;
}

DAVector operator+(const DAVector &da_vector) {
    return da_vector;
}

DAVector sqrt(const DAVector &da_vector) {
    DAVector res;
    ad_sqrt(&da_vector.da_vector_, &res.da_vector_);
	return res;
}

DAVector exp(const DAVector &da_vector) {
    DAVector res;
    ad_exp(&da_vector.da_vector_, &res.da_vector_);
	return res;
}

DAVector log(const DAVector &da_vector) {
    DAVector res;
    ad_log(&da_vector.da_vector_, &res.da_vector_);
	return res;
}

DAVector tan(const DAVector &da_vector) {
    DAVector res;
    res = sin(da_vector)/cos(da_vector);
    return res;
}

/** \brief Calculate the arctan of the given DA vector.
 * Use formula 4.4.34 and 4.4.42 in the Handbook of Mathematical Functions
 * by Milton Abramowitz and Irene A. Stegun.
 * Assume z0 is the constant part of the given DA vector z and dz = z - z0.
 * In 4.4.34, let z1 = z0 and z2 = dz/(1+z*z0), then the right hand size of
 * 4.4.34 is atan(z). So we get atan(z) = atan(z0) + atan(z2). Use 4.4.42 to
 * calculate atan(z2).
 * \param[in] DA vector.
 * \return DA vector.
 *
 */
DAVector atan(const DAVector &da_vector) {
    double da_cons = da_vector.con();
    DAVector da_high = (da_vector - da_cons)/(da_vector*da_cons + 1);
    DAVector result = atan(da_cons);
    DAVector da_high2 = da_high*da_high;
    for(double i=1, j=-1; i<da_vector.order()+1; i+=2) {
        j *= -1;
        result += da_high*j/i;
        da_high *= da_high2;
    }
    return result;
}

DAVector asin(const DAVector &da_vector) {
    double cons = da_vector.con();
    if(cons<-1 or cons>1) {
        throw std::domain_error("The constant part of the input DA vector of asin() should be within [-1,1].");
    }
    DAVector x = da_vector/(1+sqrt(1-da_vector*da_vector));
    DAVector result = 2*atan(x);
    return result;
}

DAVector acos(const DAVector &da_vector) {
    double cons = da_vector.con();
    if(cons<-1 or cons>1) {
        throw std::domain_error("The constant part of the input DA vector of acos() should be within [-1,1].");
    }
    double half_pi = 1.57079632679489661923132169163975144209858469968755291048;
    DAVector result = half_pi-asin(da_vector);
    return result;
}

DAVector sin(const DAVector &da_vector) {
    double cons = da_vector.con();
    DAVector pol = da_vector - cons;
    DAVector high_order = pol;
    int order = da_vector.order();
    std::vector<double> coefs(order+1);
    coefs.at(0) = sin(cons);
    coefs.at(1) = cos(cons);
    for(int i=2; i<order+1; ++i) {
        coefs.at(i) = -coefs.at(i-2)/(i*(i-1));
    }
    DAVector result = coefs.at(0);
    for(int i=1; i<order+1; ++i) {
        result += coefs.at(i)*high_order;
        high_order *= pol;
    }
    return result;
}

DAVector cos(const DAVector &da_vector) {
    double cons = da_vector.con();
    DAVector pol = da_vector - cons;
    DAVector high_order = pol;
    int order = da_vector.order();
    std::vector<double> coefs(order+1);
    coefs.at(0) = cos(cons);
    coefs.at(1) = -sin(cons);
    for(int i=2; i<order+1; ++i) {
        coefs.at(i) = -coefs.at(i-2)/(i*(i-1));
    }
    DAVector result = coefs.at(0);
    for(int i=1; i<order+1; ++i) {
        result += coefs.at(i)*high_order;
        high_order *= pol;
    }
    return result;
}

DAVector sinh(const DAVector &da_vector) {
    double cons = da_vector.con();
    DAVector pol = da_vector - cons;
    DAVector high_order = pol;
    int order = da_vector.order();
    std::vector<double> coefs(order+1);
    coefs.at(0) = sinh(cons);
    coefs.at(1) = cosh(cons);
    for(int i=2; i<order+1; ++i) {
        coefs.at(i) = coefs.at(i-2)/(i*(i-1));
    }
    DAVector result = coefs.at(0);
    for(int i=1; i<order+1; ++i) {
        result += coefs.at(i)*high_order;
        high_order *= pol;
    }
    return result;
}

DAVector cosh(const DAVector &da_vector) {
    double cons = da_vector.con();
    DAVector pol = da_vector - cons;
    DAVector high_order = pol;
    int order = da_vector.order();
    std::vector<double> coefs(order+1);
    coefs.at(0) = cosh(cons);
    coefs.at(1) = sinh(cons);
    for(int i=2; i<order+1; ++i) {
        coefs.at(i) = coefs.at(i-2)/(i*(i-1));
    }
    DAVector result = coefs.at(0);
    for(int i=1; i<order+1; ++i) {
        result += coefs.at(i)*high_order;
        high_order *= pol;
    }
    return result;
}

DAVector tanh(const DAVector &da_vector) {
    DAVector res = sinh(da_vector)/cosh(da_vector);
    return res;
}

double abs(const DAVector &da_vector) {
    double r = 0;
    ad_abs(&da_vector.da_vector_, &r);
    return r;
}

double abs(const std::complex<DAVector> &complex_dav) {
    double r1, r2;
    r1 = abs(complex_dav.real());
    r2 = abs(complex_dav.imag());
    return r1>r2?r1:r2;
}


DAVector pow_pos(const DAVector &da_vector, const int order) {
    DAVector res;
    if (order==0) {
        res.reset_const(1);
    }
    else if (order==1) {
        res = da_vector;
    }
    else if (order&1) { //Odd order
        res = da_vector * pow_pos(da_vector*da_vector, order/2);
    }
    else { //Even order
        res = pow_pos(da_vector*da_vector, order/2);
    }
    return res;
}

DAVector pow(const DAVector &da_vector, const int order) {
    DAVector res;
    if (order<0) {
        if (da_vector.con()==0) {
            res.reset_const(INFINITY);
            return res;
        }
        else {
            res = pow_pos(da_vector, -order);
            res = 1/res;
        }
    }
    else {
        res = pow_pos(da_vector, order);
    }
    return res;
}


DAVector pow(const DAVector &da_vector, const double order) {
    DAVector res;
    if (std::floor(order) == order) {//order is integer
        res = pow(da_vector, static_cast<int>(order));
    }
    else {
        double bas = da_vector.con();
        if(bas>0) {
            res = exp(order*log(da_vector));    //Is there a better way?
        }
        else if(bas==0) {
            res.reset_const(INFINITY);
        }
        else {
            res.reset_const(NAN);
        }
    }

    return res;
}

std::ostream& operator<<(std::ostream &os, const DAVector &da_vector) {
    print_vec(da_vector.da_vector_, os);
    return os;
}

std::ostream& operator<<(std::ostream &os, const std::complex<DAVector> &cd_vector) {
    print_vec(cd_vector.real().da_vector_, cd_vector.imag().da_vector_,os);
    return os;
}


void _ludcmp(std::vector<std::vector<double>> &a, const int n, std::vector<int> &idx, int &d){
	int imax;
	double big, dum, sum, temp;
	std::vector<double> vv(n);      //the implicit scaling of each row.
	double tiny = 1.0e-20;

	d = 1; 		//No row interchanges yet.
	for (int i=0; i<n; ++i){
		big = 0.0;
		for (int j=0; j<n; ++j){
			if((temp=fabs(a[i][j]))>big) big = temp;
		}
		if (big<tiny) {			//The largest element is zero
			std::cout<<"Singular matrix in routine LUDcmp"<<std::endl;
			exit(EXIT_FAILURE);
		}

		vv[i] = 1/big;		//Scaling
	}

	for (int j=0; j<n; ++j){		//Loop over columns of Crout's method
		for (int i=0; i<j; ++i){
			sum = a[i][j];
			for (int k=0; k<i; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}

		//Search for largest pivot element
		big = 0.0;
		for (int i=j; i<n; ++i){
			sum = a[i][j];
			for (int k=0; k<j; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ( (dum=vv[i]*fabs(sum))>=big){
				big = dum;
				imax = i;
			}
		}

		if (j != imax){			// Do we need to interchange rows?
			for (int k=0; k<n; ++k){	//Yes
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;				//change the parity of d;
			vv[imax] = vv[j];	//Interchange the scale factor;
		}

		idx[j]=imax;

		if (fabs(a[j][j])<tiny)	a[j][j] = tiny;	//If the pivot element is zero, submitted by a tiny value

		if	(j!=(n-1)){						//Divide by the pivot element
			dum = 1.0/(a[j][j]);
			for (int i=j+1; i<n; ++i) a[i][j] *= dum;
		}
	}										//Go back to the next column in the reduction
}

//Solves the set of n linear equations AX = B. Here the input a[n][n] is not the matrix A, but
//its LU decomposition. And idx[n] is the permutation vector. b[n] is the right-hand size vector
//B, and returns the solution vector X.
//This algorithm can be found in Numerical Recipes in C, 2nd ed.
void _lubksb(std::vector<std::vector<double> >&a, const int n, std::vector<int> &idx, std::vector<double> &b){

	int ii=-1, ip;
	double sum;
	for(int i=0; i<n; ++i){
		ip = idx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii+1)
			for(int j=ii; j<=i-1; ++j)	sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i] = sum;
	}

	for (int i=n-1; i>=0; --i){
		sum = b[i];
		for(int j=i+1;j<n; ++j)	{
			sum -= a[i][j]*b[j];
		}

		b[i] = sum/a[i][i];
	}
}

/** \brief Calculate the inverse of a square matrix by LU decomposition. The input matrix is destroyed after calculation.
 *
 * \param a Stores the input matrix.
 * \param n Dimension of the matrix is n*n.
 * \return y Stores the output matrix.
 *
 */
void _inv_matrix(std::vector<std::vector<double> > &a, const int n, std::vector<std::vector<double> > &y) {
    int d = 0;
    std::vector<double> col(n);
    std::vector<int> idx(n);
    _ludcmp(a, n, idx, d);
    for(int j=0; j<n; ++j) {
        for(int i=0; i<n; ++i) col[i] = 0;
        col[j] = 1;
        _lubksb(a, n, idx, col);
        for(int i=0; i<n; ++i) y[i][j] = col[i];
    }
}

/** \brief Calculate the inverse map of a given map
 *
 * \param ivecs Saves the input map with dim da vectors.
 * \param dim   Number of da vectors and number of bases used in the map.
 * \param ovecs Saves the output map.
 * \return void.
 *
 */
void inv_map(std::vector<DAVector> &ivecs, int dim, std::vector<DAVector> &ovecs) {
    assert(dim<=DAVector::dim()&&"Wrong dimension of map in inv_map!");
    for(auto& v : ivecs) assert(v.con()<std::numeric_limits<double>::min()&&"Constant part of the input map is NOT zero in inv_map!");

    std::vector<std::vector<double> >lin_matrix(dim, std::vector<double>(dim)); //Save the linear coefficients.

    //fetch linear map coefficient;
    std::vector<int> c(dim, 0);
    for(int i=0; i<dim; ++i) {
        for(int j=0; j<dim; ++j) {
            c.at(j) = 1;
            lin_matrix[i][j] = ivecs.at(i).element(c);
            c.at(j) = 0;
        }
    }

    //nonlinear part of the given map;
    std::vector<DAVector> nlin_map;
    unsigned int da_order = DAVector::order();
    for(auto& v : ivecs) {
        da_change_order(1);
        DAVector t = v;
        da_change_order(da_order);
        t = v - t;
        nlin_map.push_back(t);
    }

    std::vector<std::vector<double> >inv_lin_matrix(dim, std::vector<double>(dim)); //Save the inverse of the linear map.
    _inv_matrix(lin_matrix, dim, inv_lin_matrix);

    std::vector<DAVector> inv_lin_map;
    for(int i=0; i<dim; ++i) {
        DAVector t;
        t.reset_const(0);
        for(int j=0; j<dim; ++j) {
            if (fabs(inv_lin_matrix[i][j])>1e-16)
                t = t + inv_lin_matrix[i][j]*da[j];
        }
        inv_lin_map.push_back(t);
    }

    std::vector<DAVector> tmp_ovecs;
    for(auto t : inv_lin_map) tmp_ovecs.push_back(t);
    std::vector<DAVector> tmp(dim);

    for(int i=0; i<DAVector::order()+1; ++i) {
        da_composition(nlin_map, tmp_ovecs, tmp);
        for(int j=0; j<dim; ++j)
            tmp_ovecs.at(j) = da[j] - tmp.at(j);
        da_composition(inv_lin_map, tmp_ovecs, tmp);
        for(int j=0; j<dim; ++j)
            tmp_ovecs.at(j) = tmp.at(j);
    }
    for(int i=0; i<dim; ++i) {
        ovecs.at(i) = tmp_ovecs.at(i);
    }
}


/** \brief error function for a DA vector
 * Calculate the error function value of the input DA vector. The erf() in math.h is used to calculate the cosntant
 * apart. 2/sqrt(pi)*exp(-x^2) is the derivative of the error function, by which the high order terms could be calculated.
 *
 * \param x The input DA vector.
 * \return erf(x) as a DA vector.
 *
 */
DAVector erf(const DAVector& x) {
  double coef = 1.1283791670955125585607;  //coef = 2/sqrt(pi);
  double cc = x.con();
  DAVector dal = cc + da[0];
  dal = exp(-1*dal*dal);
  dal = da_int(dal, 0);
  DAVector ada = x - cc;
  DAVector da_erf;
  da_substitute(dal, 0, ada, da_erf);
  return da_erf*coef + erf(cc);
}
