/**
 * @file da.cc
 * @brief DA wrapper for the extended tpsa code.
 * @details Defines methods for DAVectors.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#include "../include/da.h"
#include <cmath>
#include "../include/tpsa_extend.h"

Base da; /**< Bases for DA calculations. The i-th base can be accessed as da[i].  */

DAVector::DAVector() {ad_alloc(&da_vector_);}

DAVector::DAVector(const DAVector& da_vector) {
      ad_alloc(&da_vector_);
      ad_copy(&da_vector.da_vector_, &da_vector_);
}

DAVector::DAVector(DAVector&& da_vector) {
    ad_assign(da_vector_);
    std::swap(da_vector_, da_vector.da_vector_);
    ad_free(&da_vector.da_vector_);
}

DAVector::~DAVector() {
    ad_free(&da_vector_);
}

void DAVector::print() const { ad_print(&da_vector_);}      /**< Print out a DA vector. */
double DAVector::con() const {return ad_con(da_vector_);}   /**< Return the constant element of a DA vector. */

DAVector& DAVector::operator=(DAVector&& da_vector) {
    if (this==&da_vector) return *this;
    std::swap(this->da_vector_, da_vector.da_vector_);
    ad_reset(&da_vector.da_vector_);
    return *this;
}

void DAVector::reset() { ad_reset_vector(da_vector_);}   /**< Reset all element to zero, vector length unchanged. */
/// Reset the value to the given number. Vector length is set to 1.
void DAVector::reset_const(double x) {ad_reset_const(da_vector_, x);}
int DAVector::dim() { return ad_dim();}          /**< Return the DA base number. */
DAVector& DAVector::operator=(const DAVector &da_vector) { ad_copy(&da_vector.da_vector_, &da_vector_); return *this;}

/// Return the length of the DA vector.
unsigned int DAVector::length() const {
  	unsigned int length;
  	ad_length(&da_vector_, &length);
  	return length;
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
void DAVector::element(unsigned int &i, unsigned int *c, double &elem) const {
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
    base.clear();
    for(unsigned int i = 0; i<n; ++i) {
        DAVector v;
        double zero = 0.0;
		ad_var(&v.da_vector_, &zero, &i);
		base.push_back(v);
    }
}

/// Set the DA base vectors.
void Base::set_base() {
    unsigned int n = DAVector::dim();
    set_base(n);
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
int da_init(unsigned int& da_order, unsigned int& num_da_variables, unsigned int& num_da_vectors) {
	ad_init(&num_da_variables, &da_order);
	ad_reserve(&num_da_vectors);
	da.set_base(num_da_variables);
	return 0;
}

///Temporarily reduce the DA order.
int da_reduce_order(unsigned int new_order) {ad_change_order(new_order); return 0;}

///Restore the original DA order.
int da_restore_order(){ad_restore_order(); return 0;}

/** \brief Subscribe the given value to the specific base in a DA vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal number of the base.
 * \param[in] x The alue to subscribe.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_subscribe_const(const DAVector &iv, unsigned int base_id, double x, DAVector &ov) {
    ad_subscribe(iv.da_vector_, base_id, x, ov.da_vector_);
}

/** \brief Subscribe a DA vector to a specific base in a given DA Vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal number of the base.
 * \param[in] v The DA vector to subscribe.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_subscribe(const DAVector &iv, unsigned int base_id, const DAVector &v, DAVector &ov) {
    ad_subscribe(iv.da_vector_, base_id, v.da_vector_, ov.da_vector_);
}

/** \brief Subscribe DA vectors to the specific bases in a DA vector.
 *
 * \param[in] iv A DA vector.
 * \param[in] base_id The ordinal numbers of the bases.
 * \param[in] v The DA vectors to subscribe. The size of v should be equal to the size of base_id.
 * \param[out] ov The DA vector that saves the result.
 * \return void.
 *
 */
void da_subscribe(const DAVector &iv, std::vector<unsigned int> &base_id, std::vector<DAVector> &v, DAVector &ov) {
    std::vector<TVEC> ad_v;
    for(auto&& i : v) ad_v.push_back(i.da_vector_);
    ad_subscribe(iv.da_vector_, base_id, ad_v, ov.da_vector_);
}

/** \brief Subscribe DA vectors to the specific bases in a group of DA vectors.
 *
 * \param ivecs A group of DA Vectors.
 * \param base_id The ordinal numbers for the bases.
 * \param v The DA vectors to subscribe. The size of v should be equal to the size of base_id.
 * \param ovecs The DA vectors that saves the result. The size of ovecs should be equal to the size of ivecs.
 * \return
 *
 */
void da_subscribe(std::vector<DAVector> &ivecs, std::vector<unsigned int> &base_id, std::vector<DAVector> &v,
                  std::vector<DAVector> &ovecs) {
    std::vector<TVEC> ad_iv, ad_v, ad_ov;
    for(auto&& i : v) ad_v.push_back(i.da_vector_);
    for(auto&& i : ivecs) ad_iv.push_back(i.da_vector_);
    for(auto&& i : ovecs) ad_ov.push_back(i.da_vector_);
    ad_subscribe(ad_iv, base_id, ad_v, ad_ov);
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

DAVector sin(const DAVector &da_vector) {
    DAVector res;
    ad_sin(&da_vector.da_vector_, &res.da_vector_);
	return res;
}

DAVector cos(const DAVector &da_vector) {
    DAVector res;
    ad_cos(&da_vector.da_vector_, &res.da_vector_);
	return res;
}

double abs(const DAVector &da_vector) {
    double r = 0;
    ad_abs(&da_vector.da_vector_, &r);
    return r;
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

