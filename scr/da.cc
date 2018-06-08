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
int DAVector::n_element() const {return ad_n_element(da_vector_);} /**< Return the number of non-zero element a DA vector. */

DAVector& DAVector::operator=(DAVector&& da_vector) {
    if (this==&da_vector) return *this;
    std::swap(this->da_vector_, da_vector.da_vector_);
    ad_reset(&da_vector.da_vector_);
    return *this;
}

DAVector& DAVector::operator=(const DAVector &da_vector) { ad_copy(&da_vector.da_vector_, &da_vector_); return *this;}

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
int da_init(unsigned int da_order, unsigned int num_da_variables, unsigned int num_da_vectors) {
    assert(da_order>0&&num_da_variables>0&&num_da_vectors>0&&"Wrong parameters when initialize DA environment!");
	ad_init(&num_da_variables, &da_order);
	ad_reserve(&num_da_vectors);
	da.set_base(num_da_variables);
	return 0;
}

///Temporarily change the DA order.
int da_change_order(unsigned int new_order) {ad_change_order(new_order); return 0;}

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

void _ludcmp(std::vector<std::vector<double>> &a, const int n, std::vector<int> &idx, int &d){
	int imax;
	double big, dum, sum, temp;
	std::vector<double> vv(n);      //the implicit scaling of each row.
//	double *vv = new double[n];		//the implicit scaling of each row.
	double tiny = 1.0e-20;

	d = 1; 		//No row interchanges yet.
	for (int i=0; i<n; ++i){
		big = 0.0;
		for (int j=0; j<n; ++j){
			if((temp=fabs(a[i][j]))>big) big = temp;
		}
		if (big==0) {			//The largest element is zero
			std::cout<<"Singular matrix in routine LUDcmp"<<std::endl;
//			delete[] vv;
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

		if (a[j][j]==0)	a[j][j] = tiny;	//If the pivot element is zero, submitted by a tiny value

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
    int da_order = DAVector::order();
    for(auto& v : ivecs) {
        da_change_order(1);
        DAVector t = v;
        da_change_order(da_order);
//        da_restore_order();
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
