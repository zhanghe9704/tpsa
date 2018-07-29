/**
 * @file tpsa_extend.cc
 * @brief Extend the original tpsa.cpp file by Dr. Lingyun Yang.
 * @details Put the tpsa.cpp from Dr. Lingyun Yang in the save directory with this file.
 * But do NOT included in the project. Otherwise there will be multiple definition of many functions.
 * Some static variables defined in tpsa.cpp are used in this file.
 * Because one cannot reference a static variable outside the file, tpsa.cpp is included, instead of tpsa.h.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#include "../include/tpsa_extend.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include "tpsa.cpp"


static unsigned int ad_flag = 0; ///< The index of the next available TPS vector
static unsigned int ad_end = 0;  ///< The index of the last available TPS vector.
///Linked list for TPSA memory management. When ad_flag == adlist[ad_end], memory runs out!
static std::vector <unsigned int> adlist;
static TNVND gnd_record = 0; ///< Temporarily record the TSP order, used only when reducing and restoring the TSP order.

/** \brief Return the order of the TPS environment.
 *
 * \return order of the TPS environment.
 *
 */
int ad_order() {
    return gnd;
}

/** \brief Return the full length of the TPS with the current order and dimension
 *
 * \return full length of the TPS
 *
 */
int ad_full_length() {
    return FULL_VEC_LEN;
}

/** \brief Reduce the TPS order to the specific value.
 *
 * \param new_order New TPS order. The new order should be lower than the current order.
 * \return void
 *
 */
void ad_change_order(unsigned int new_order) {
    if (0==gnd_record) gnd_record = gnd;
    if (new_order<gnd_record) {
        gnd = new_order;
        FULL_VEC_LEN = comb_num(gnv+gnd, gnd);
    }
}

/** \brief Return the number of TPS bases.
 *
 * \return TPS base number.
 *
 */

unsigned int ad_dim() {
    return gnv;
}

/** \brief Restore the TPS order to the original value.
 *
 * \return void
 *
 */
void ad_restore_order() {
    if (gnd_record>gnd) {
        gnd = gnd_record;
        FULL_VEC_LEN = comb_num(gnv+gnd, gnd);
    }
}

void ad_list_order(const TVEC iv) {
    TNVND* p = base;
    for (size_t i=0; i<adveclen[iv]; ++i) {
        if (std::abs(advec[iv][i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        for (size_t j = 0; j < gnv-1; ++j) {
            std::cout << (unsigned int) (*p-*(p+1)) <<' ';
            ++p;
        }
        std::cout << (unsigned int)*p++ << ' ' << i << std::endl;
    }
}

TVEC ad_pow_int_pos(const TVEC iv, std::vector<unsigned int> &power_v, const int order, unsigned int order_rec) {
    TVEC res;
    if (adveclen[power_v.at(order*order_rec)]>0) {
        res = power_v.at(order*order_rec);
    }
    else if (order==0) {
        advec[power_v.at(0)][0] = 1;
        adveclen[power_v.at(0)] = 1;
        res = power_v.at(0);
    }
    else if (order==1) {
        res = power_v.at(1*order_rec);
        if(adveclen[res]==0) ad_copy(&iv, &res);
    }
    else if (order==2) {
        res = power_v.at(2*order_rec);
        if(adveclen[res]==0) ad_mult(&iv, &iv, &res);
    }
    else if (order&1) { //Odd order
        TVEC vres;
        unsigned int order_idx = order_rec*2;
        vres = power_v.at(order_idx);
        if (adveclen[vres]==0) ad_mult(&iv, &iv, &vres);

        vres = ad_pow_int_pos(vres, power_v, order/2, order_idx);
        res = power_v.at(order_rec+(order/2)*order_idx);
        ad_mult(&iv, &vres, &res);

    }
    else {              //Even order
        TVEC vres;
        unsigned int order_idx = order_rec*2;
        vres = power_v.at(order_idx);
        if (adveclen[vres]==0) ad_mult(&iv, &iv, &vres);
        res = ad_pow_int_pos(vres, power_v, order/2, order_idx);
    }
    return res;
}

/** \brief Check the remaining capacity in the TPS pool.
 *
 * \return The maximum number of TPS vectors that can be created in the remaining place of the pool.
 *
 */
unsigned int ad_remain() {
    unsigned int i = 0;
    unsigned int tmp_flag = ad_flag;
//    while(tmp_flag != ad_end) {
//        ++i;
//        tmp_flag = adlist.at(tmp_flag);
//    }
    while(tmp_flag != adlist.at(ad_end)) {
        ++i;
        tmp_flag = adlist.at(tmp_flag);
    }
    return i;
}



/** \brief Find the value of the constant part of a given TPS vector.
 *
 * \param iv A TPS vector.
 * \return Value of the constant part of the given TPS vector.
 *
 */
double ad_con(const TVEC iv) {
    return advec[iv][0];
}

/** \brief Reset all the elements to zero and keep the length unchanged for a given TPS vector.
 *
 * \param iv A TPS vector.
 * \return void.
 *
 */
void ad_reset_vector(const TVEC iv)
{
    memset(advec[iv], 0, FULL_VEC_LEN*sizeof(double));
}

/** \brief Reset the given TPS vector to a given constant.
 * The constant part of the given TPS vector will be equal to the given constant. All higher order elements are zeros.
 * \param iv A TPS vector.
 * \param x A number.
 * \return void.
 *
 */
void ad_reset_const(const TVEC iv, double x) {
    ad_reset_vector(iv);
    advec[iv][0] = x;
    adveclen[iv] = 1;
}

unsigned int get_order_index(unsigned int i) {
    return order_index[i];
}

/** \brief Integrate a TPS w.r.t. a specific base
 *
 * \param iv The TPS to integrate.
 * \param base_id The base.
 * \param ov The output TPS.
 * \return void.
 *
 */
void ad_int(TVEC iv, unsigned int base_id, TVEC ov) {
    TVEC v;
    double x = 0;
    ad_alloc(&v);
    ad_var(&v, &x, &base_id);
    ad_mult(&iv, &v, &ov);

    TNVND* p = base;
    std::vector<unsigned int> c(gnv);

    for (size_t i=0; i<adveclen[ov]; ++i) {
        if (std::abs(advec[ov][i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        for (size_t j = 0; j < gnv-1; ++j) {
            c[j] = (unsigned int) (*p-*(p+1));
            ++p;
        }
        c[gnv-1] = (unsigned int)*p++ ;
        advec[ov][i] /= c[base_id];
    }
    ad_free(&v);
}

/** \brief Return the number of non-zero elements in a TPS.
 *
 * \param v The TPS.
 * \return Number of non-zero elements in the TPS.
 *
 */
int ad_n_element(TVEC v) {
    int n = 0;
    for(size_t i=0; i<adveclen[v]; ++i)
        if (fabs(advec[v][i])>std::numeric_limits<double>::min())
            ++n;
    return n;
}

// ***** The following functions provide alternative ones instead of the original ones in tpsa.cpp. *****

/** \brief Count the number of existing TPS vectors.
 * This is an alternative function for the original one in tpsa.cpp.
 * \return The number of existing TPS vectors.
 *
 */
unsigned int ad_count() {
    return adlist.size() - ad_remain();
}

/** \brief Multiply a constant number to a TPS vector.
 * Multiply a constant number c to the given TPS vector iv, save the result in ov.
 * This function does not change the value of iv, as in ad_mult_const.
 * This is an alternative function for the original one (ad_mult_const) in tpsa.cpp.
 * \param[in] iv A TPS vector.
 * \param[in] c A number.
 * \param[out] ov The result TPS vector as ov = c * iv.
 * \return void.
 *
 */
void ad_mult_c(const TVEC iv, double c, TVEC ov)
{
    adveclen[ov] = adveclen[iv];
    for (size_t i = 0; i < std::min(adveclen[iv], order_index[gnd+1]); ++i)
//    for (size_t i = 0; i < adveclen[iv]; ++i)
        advec[ov][i] = c*advec[iv][i];
}

//Reserve space for n TPS vectors.
/** \brief Reserve memory for n TPS vectors.
 * Allocate the memory pool for n TPS vectors and initialize the linked list (adlist, ad_flag and ad_end) for memory management.
 * This is an alternative function for the original one in tpsa.cpp.
 * \param n Number of TPS vectors.
 * \return void.
 *
 */
void ad_reserve(const unsigned int n)
{
    if (n <= 0) return;

    if (advecpool!=nullptr) {   //If the memory has already allocated before, release it first!
        delete [] advecpool[0];
        delete [] advecpool;
    }
    advecpool = new double*[n];
    advecpool[0] = new double[FULL_VEC_LEN*n];
    memset(advecpool[0], 0, FULL_VEC_LEN*n*sizeof(double));
    ad_flag = 0;
    advec.clear();
    advec.push_back(advecpool[0]);
    adlist.clear();
    for(size_t i=1; i<n; ++i) {
        adlist.push_back(i);
        advecpool[i] = &advecpool[0][FULL_VEC_LEN*i];
        advec.push_back((advecpool[i]));
    }
    adveclen.clear();
    adveclen.resize(n,0);
    adlist.push_back(n);
    ad_end = n-1;
}

/** \brief Assign a memory slot for a TPS vector.
 * Assign memory slot for a TPS vector. The length of the vector is zero.
 * This is an alternative function for the original one (ad_alloc) in tpsa.cpp.
 * \param i Return a new integer representing the assigned TPS vector.
 * \return void.
 * \note
 *     - The whole TPS vector pool was allocated by ad_reserve(n). This command does not actually allocate memory.
 * It only assign a slot to the specific TPS vector.
 *     -Use ad_alloc to assign memory slot for an TPS vector with value zero and length 1.
 *
 */
void ad_assign(unsigned int &i) {
    if (ad_flag == adlist.at(ad_end)) {
        std::cerr << "Run out of vectors" << std::endl;
        exit(-1);
    }
    i = ad_flag;
    adveclen[i] = 0; //zero length. Use ad_alloc to set the value 0 and the length 1.
    ad_flag = adlist.at(ad_flag);
}

/** \brief Composition of a group of TPS vectors with a group of numbers.
 * Submitting the numbers into the bases of a group of TPS vectors, and save the result into another group of TPS vectors.
 * Call the first group of TPS vectors, ivecs, as f, which includes n TPS vectors: f_1, f_2, ..., f_n.
 * A group of m numbers are saved in v, where m is equal to the number of TPS bases.
 * Call the second group of TPS vectors, ovecs, as g, which also includes n TPS vectors: g_1, g_2, ..., g_n.
 * g = f(v), or g_1 = f_1(v), g_2 = f_2(v), ..., g_n = f_n(v).
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] ivecs A group of TPS vectors, saved as an std::vector.
 * \param[in] v A group of TPS vectors, saved as an std::vector.
 * \param[out] ovecs A group of TPS vectors, saved as an std::vector.
 * \return void.
 *
 */
void ad_composition(std::vector<TVEC> &ivecs, std::vector<double> &v, std::vector<double> &ovecs) {
    assert(gnv==v.size()&&"Error in ad_composition: No. of TPS vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size()==ovecs.size()&&"Error in ad_composition: No. of input vectors NOT EQUAL to No. of output vectors!");

    TNVND* p = base;
    std::vector<unsigned int> c(gnv);
    std::vector<unsigned int> bv(gnv);
    std::vector< std::vector<double> > power_vv (gnv, std::vector<double>(gnd+1));
    for(auto& power_v : power_vv) {
        power_v.at(0) = 1;
    }
    for (unsigned int idx = 0; idx<gnv; ++idx) {
        double val = v.at(idx);
        for(unsigned int i = 1; i<gnd+1; ++i) {
            power_vv.at(idx).at(i) = val*power_vv.at(idx).at(i-1);
        }
    }
    //Find the length of the longest TPS vector.
    unsigned int veclen_max = 0;
    for(auto iv : ivecs)
        if (adveclen[iv]>veclen_max) veclen_max = adveclen[iv];

    //Copy the constant element to output vectors
    for(unsigned int iv = 0; iv<ivecs.size(); ++iv)
        ovecs.at(iv) = advec[ivecs.at(iv)][0];

    for (size_t j = 0; j < gnv-1; ++j) {
        c.at(j) = (unsigned int) (*p-*(p+1));
        ++p;
    }
    c.at(gnv-1) = (unsigned int)*p++ ;
    unsigned int vec_size = ivecs.size();
    for (size_t i=1; i<veclen_max; ++i) { //loop over all elements, except for the constant element.
        bool product_flag = true; //Calculate the product.
        bool c_flag = true;
        double product = 1;
        unsigned int zero_coef = 0;
        for(unsigned int iv = 0; iv<vec_size; ++iv) {
            if (i>=adveclen[ivecs.at(iv)]) {
                ++zero_coef;
                continue;
            }
            double coef = advec[ivecs.at(iv)][i];
            if (std::abs(coef) < std::numeric_limits<double>::min()) {
                ++zero_coef;
                continue;
            }

            if(c_flag) {
                for (size_t j = 0; j < gnv-1; ++j) {
                    c.at(j) = (unsigned int) (*p-*(p+1));
                    ++p;
                }
                c.at(gnv-1) = (unsigned int)*p++ ;
                c_flag = false;
            }

            if (product_flag) {
                for(unsigned int id=0; id<gnv; ++id) product *= power_vv.at(id).at(c.at(id));
                product_flag = false;
            }
            ovecs.at(iv) += product*coef;
        }
        if (zero_coef == vec_size) p += gnv;
    }
}

/** \brief Composition of two groups of TPS vectors.
 * Calculate the composition of two groups of TPS vectors, and save the result into the third group of TPS vectors.
 * Call the first group of TPS vectors, ivecs, as f, which includes n TPS vectors: f_1, f_2, ..., f_n.
 * The second group of TPS vectors, v, should contain the same number of TPS vectors as the number of TPS bases.
 * Call the third group of TPS vectors, ovecs, as g, which also includes n TPS vectors: g_1, g_2, ..., g_n.
 * g = f(v), or g_1 = f_1(v), g_2 = f_2(v), ..., g_n = f_n(v).
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] ivecs A group of TPS vectors, saved as an std::vector.
 * \param[in] v A group of TPS vectors, saved as an std::vector.
 * \param[out] ovecs A group of TPS vectors, saved as an std::vector.
 * \return void.
 *
 */
void ad_composition(std::vector<TVEC> &ivecs, std::vector<TVEC> &v, std::vector<TVEC> &ovecs) {
    assert(gnv==v.size()&&"Error in ad_composition: No. of TPS vectors NOT EQUAL to No. of bases!");
    assert(ivecs.size()==ovecs.size()&&"Error in ad_composition: No. of input vectors NOT EQUAL to No. of output vectors!");

    TNVND* p = base;
    unsigned int* c = new unsigned int[gnv];
    unsigned int* bv = new unsigned int[gnv];

    std::vector< std::vector<TVEC> > power_vv(gnv, std::vector<TVEC>(gnd+1));
    for(auto& power_v : power_vv) {
        for(auto& vec : power_v) {
            TVEC idx;
            ad_assign(idx);
            vec = idx;
        }
        TVEC idx = power_v.at(0);
        advec[idx][0] = 1;
        adveclen[idx] = 1;
    }

    for(unsigned int i=0; i<gnv; ++i) {
        ad_free(&power_vv.at(i).at(1));
        power_vv.at(i).at(1) = v.at(i);
    }

    TVEC tmp, product;
    ad_alloc(&tmp);
    ad_alloc(&product);

    unsigned int veclen_max = 0;
    for(auto iv : ivecs)
        if (adveclen[iv]>veclen_max) veclen_max = adveclen[iv];

    //Copy the constant element to output vectors
    for(auto ov : ovecs) ad_reset(&ov);
    for(unsigned int iv = 0; iv<ivecs.size(); ++iv)
        advec[ovecs.at(iv)][0] = advec[ivecs.at(iv)][0];

    for (size_t j = 0; j < gnv-1; ++j) {
        c[j] = (unsigned int) (*p-*(p+1));
        ++p;
    }
    c[gnv-1] = (unsigned int)*p++ ;

    unsigned int vec_size = ivecs.size();
    for (size_t i=1; i<veclen_max; ++i) { //loop over all elements, except for the constant element.
        bool c_flag = true; //Calculate c;
        bool product_flag = true; //Calculate the product.
        unsigned int zero_coef = 0;
        for(unsigned int iv = 0; iv<vec_size; ++iv) {
            if (i>=adveclen[ivecs.at(iv)]) {
                ++zero_coef;
                continue;
            }
            if (std::abs(advec[ivecs.at(iv)][i]) < std::numeric_limits<double>::min()) {
                ++zero_coef;
                continue;
            }
            if (c_flag) {
                for (size_t j = 0; j < gnv-1; ++j) {
                    c[j] = (unsigned int) (*p-*(p+1));
                    ++p;
                }
                c[gnv-1] = (unsigned int)*p++ ;

                for(unsigned int id=0; id<gnv; ++id) {
                    if (c[id]>0) {
                        ad_pow_int_pos(v.at(id), power_vv.at(id), c[id], 1);
                    }
                }
                c_flag = false;  //Already calculated c.
            }

            double coef = advec[ivecs.at(iv)][i];
            if (product_flag) {
                ad_reset_vector(product);
                advec[product][0] = 1;
                adveclen[product] = 1;
                for(unsigned int id=0; id<gnv; ++id) {
                    if(c[id]>0) {
                        ad_mult(&product, &power_vv.at(id).at(c[id]), &tmp);
                        std::swap(product, tmp);
                    }
                }
                product_flag = false;
            }
            ad_mult_c(product, coef, tmp);
            for(unsigned int idx=0; idx<adveclen[tmp]; ++idx)
                advec[ovecs.at(iv)][idx] += advec[tmp][idx];
        }
        if (zero_coef==vec_size) p += gnv;
    }

    for(auto ov : ovecs) {
        adveclen[ov] = 1;
        for(int i=order_index[gnd+1]-1; i>=0; --i) {
            if(std::abs(advec[ov][i])>std::numeric_limits<double>::min()) {
                adveclen[ov] = i+1;
                break;
            }
        }
//        std::cout<<ov<<' '<<adveclen[ov]<<std::endl;
    }

    //Release the temporary TPS vectors.
    ad_free(&tmp);
    ad_free(&product);

    for(auto power_v : power_vv) {
        ad_free(&power_v.at(0));
        for(unsigned int i=2; i<(gnd+1); ++i) ad_free(&power_v.at(i));
    }
    delete[] c;
    delete[] bv;
}

/** \brief Subscribe a group of TPS vectors to the specific bases of another group of TPS vectors. Result saved into the third group of TPS vectors.
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] ivecs The TPS vectors into which other TPS vectors are subscribed.
 * \param[in] base_id The indexes of the bases.
 * \param[in] v The TPS vectors to be subscribed. The number of the TPS vectors should be equal to the number of indexes.
 * \param[out] ovecs The TPS vectors that save the result.
 * \return void.
 *
 */
void ad_subscribe(std::vector<TVEC> &ivecs, std::vector<unsigned int> &base_id, std::vector<TVEC> &v, std::vector<TVEC> &ovecs) {
    assert(base_id.size()==v.size()&&"Error in ad_subscribe: No. of TPS vectors NOT EQUAL to No. of base ids!");
    assert(ivecs.size()==ovecs.size()&&"Error in ad_subscribe: No. of input vectors NOT EQUAL to No. of output vectors!");
    TNVND* p = base;
    unsigned int* c = new unsigned int[gnv];
    unsigned int* bv = new unsigned int[gnv];
    unsigned int nv = v.size();

    std::vector< std::vector<TVEC> > power_vv(nv, std::vector<TVEC>(gnd+1));
    for(auto& power_v : power_vv) {
        for(auto& vec : power_v) {
            TVEC idx;
            ad_assign(idx);
            vec = idx;
        }
        TVEC idx = power_v.at(0);
        advec[idx][0] = 1;
        adveclen[idx] = 1;
    }

    for(unsigned int i=0; i<nv; ++i) {
        ad_free(&power_vv.at(i).at(1));
        power_vv.at(i).at(1) = v.at(i);
    }
    for(auto ov : ovecs) ad_reset(&ov);
    TVEC tmp, product;
    ad_alloc(&tmp);
    ad_alloc(&product);

    std::vector<unsigned int> rc(nv);

    unsigned int veclen_max = 0;
    for(auto iv : ivecs)
        if (adveclen[iv]>veclen_max) veclen_max = adveclen[iv];

    for (size_t i=0; i<veclen_max; ++i) { //loop over all elements.
        bool c_flag = true; //Calculate c;
        bool k_flag = true; //Calculate k;
        bool sub_flag = false; //No subscribing.
        bool product_flag = true; //Calculate the product.
        std::fill(rc.begin(), rc.end(), 0);
        unsigned int k = 0;
        unsigned int order; //Highest order to be multiplied in the TPS vector v^n.
        unsigned int idx_limit; //Limit of the index of the terms to be multiplied in v^n.
        for(unsigned int iv = 0; iv<ivecs.size(); ++iv) {
            if (i>=adveclen[iv]) continue;
            if (std::abs(advec[ivecs.at(iv)][i]) < std::numeric_limits<double>::min()) {
                p += gnv;
                continue;
            }
            if (c_flag) {
                for (size_t j = 0; j < gnv-1; ++j) {
                    c[j] = (unsigned int) (*p-*(p+1));
                    ++p;
                }
                c[gnv-1] = (unsigned int)*p++ ;

                for(unsigned int id=0; id<nv; ++id) {
                    if (c[base_id.at(id)]>0) {
                        ad_pow_int_pos(v.at(id), power_vv.at(id), c[base_id.at(id)], 1);
                        rc.at(id) = c[base_id.at(id)];
                        c[base_id.at(id)] = 0;
                        sub_flag = true;
                    }
                }
                c_flag = false;  //Already calculated c.
            }

            if (sub_flag) {
                if (k_flag) {
                    unsigned int d = 0;
                    for (unsigned int i = 0; i < gnv; ++i) {
                        d += c[i];
                    }                           //d - total order of the term after eliminating the specific base.
                    order = gnd-d; //Highest order to be multiplied in the TPS vector v^n.
                    idx_limit = order_index[order+1]; //Limit of the index of the terms to be multiplied in v^n.

//                    unsigned int k = 0;
                    for (unsigned int i = 0; i < gnv; ++i){
                        bv[i] = d;
                        d -= c[i];
                        k += H[gnv-i][bv[i]];
                    }                           //k - index of the term after eliminating the specific base.
                    k_flag = false;  //Already calculated k.
                }
                double coef = advec[ivecs.at(iv)][i];
                if (k>0) {
                    if(product_flag) {
                        ad_reset_vector(product);
                        advec[product][0] = 1;
                        adveclen[product] = 1;
                        ad_change_order(order);
                        for(unsigned int id=0; id<nv; ++id) {
                            if(rc.at(id)>0) {
                                ad_mult(&product, &power_vv.at(id).at(rc.at(id)), &tmp);
                                std::swap(product, tmp);
                            }
                        }
                        ad_restore_order();
                        product_flag = false;
                    }

                    ad_mult_c(product, coef, tmp);
                    if (std::abs(advec[tmp][0]) > std::numeric_limits<double>::min())
                        advec[ovecs.at(iv)][k] += advec[tmp][0];
                    for (unsigned int idx = 1; idx<adveclen[tmp] && idx<idx_limit; ++idx) {
                        advec[ovecs.at(iv)][prdidx[k][idx]] += advec[tmp][idx];
                    }
                }
                else {
                    if (product_flag) {
                        ad_reset_vector(product);
                        advec[product][0] = 1;
                        adveclen[product] = 1;
                        for(unsigned int id=0; id<nv; ++id) {
                            if(rc.at(id)>0) {
                                ad_mult(&product, &power_vv.at(id).at(rc.at(id)), &tmp);
                                std::swap(product, tmp);
                            }
                        }
                        product_flag = false;
                    }
                    ad_mult_c(product, coef, tmp);
                    for(unsigned int idx=0; idx<adveclen[tmp]; ++idx)
                        advec[ovecs.at(iv)][idx] += advec[tmp][idx];
                }
            }
            else {              //No subscribing to this term.
                advec[ovecs.at(iv)][i] += advec[ivecs.at(iv)][i];
            }
        }
    }

    for(auto ov : ovecs) {
        adveclen[ov] = 1;
        for(int i=order_index[gnd+1]-1; i>=0; --i) {
            if(std::abs(advec[ov][i])>std::numeric_limits<double>::min()) {
                adveclen[ov] = i+1;
                break;
            }
        }
    }

    //Release the temporary TPS vectors.
    ad_free(&tmp);
    ad_free(&product);

    for(auto power_v : power_vv) {
        ad_free(&power_v.at(0));
        for(unsigned int i=2; i<(gnd+1); ++i) ad_free(&power_v.at(i));
    }
    delete[] c;
    delete[] bv;
}

/** \brief Subscribe a group of TPS vectors to the specific bases of a given TPS vector. Save the result into another TPS vector.
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] ivecs The TPS vector into which other TPS vectors are subscribed.
 * \param[in] base_id The indexes of the bases.
 * \param[in] v The TPS vectors to be subscribed. The number of the TPS vectors should be equal to the number of indexes.
 * \param[out] ovecs The TPS vector that saves the result.
 * \return void.
 *
 */
void ad_subscribe(const TVEC iv, std::vector<unsigned int> &base_id, std::vector<TVEC> &v, TVEC ov) {

    assert(base_id.size()==v.size()&&"Error in ad_subscribe: No. of TPS vectors NOT EQUAL to No. of base ids!");
    TNVND* p = base;
    unsigned int* c = new unsigned int[gnv];
    unsigned int* bv = new unsigned int[gnv];
    unsigned int nv = v.size();

    std::vector< std::vector<TVEC> > power_vv(nv, std::vector<TVEC>(gnd+1));
    for(auto& power_v : power_vv) {
        for(auto& vec : power_v) {
            TVEC idx;
            ad_assign(idx);
            vec = idx;
        }
        TVEC idx = power_v.at(0);
        advec[idx][0] = 1;
        adveclen[idx] = 1;
    }

    for(unsigned int i=0; i<nv; ++i) {
        ad_free(&power_vv.at(i).at(1));
        power_vv.at(i).at(1) = v.at(i);
    }

    ad_reset(&ov);
    TVEC tmp, product;
    ad_alloc(&tmp);
    ad_alloc(&product);

    std::vector<unsigned int> rc(nv);

    for (size_t i=0; i<adveclen[iv]; ++i) { //loop over all elements.
        if (std::abs(advec[iv][i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        for (size_t j = 0; j < gnv-1; ++j) {
            c[j] = (unsigned int) (*p-*(p+1));
            ++p;
        }
        c[gnv-1] = (unsigned int)*p++ ;

        bool flag = false;
        for(unsigned int id=0; id<nv; ++id) {
            if (c[base_id.at(id)]>0) {
                ad_pow_int_pos(v.at(id), power_vv.at(id), c[base_id.at(id)], 1);
                rc.at(id) = c[base_id.at(id)];
                c[base_id.at(id)] = 0;
                flag = true;
            }
        }
        if (flag) {
            double coef = advec[iv][i];
            unsigned int d = 0;
            for (unsigned int i = 0; i < gnv; ++i) {
                d += c[i];
            }                           //d - total order of the term after eliminating the specific base.
            unsigned int order = gnd-d; //Highest order to be multiplied in the TPS vector v^n.
            unsigned int idx_limit = order_index[order+1]; //Limit of the index of the terms to be multiplied in v^n.

            unsigned int k = 0;
            for (unsigned int i = 0; i < gnv; ++i){
                bv[i] = d;
                d -= c[i];
                k += H[gnv-i][bv[i]];
            }                           //k - index of the term after eliminating the specific base.

            if (k>0) {
                ad_reset_vector(product);
                advec[product][0] = 1;
                adveclen[product] = 1;
                ad_change_order(order);
                for(unsigned int id=0; id<nv; ++id) {
                    if(rc.at(id)>0) {
                        ad_mult(&product, &power_vv.at(id).at(rc.at(id)), &tmp);
                        std::swap(product, tmp);
//                        ad_copy(&tmp, &product);
                    }
                }
                ad_mult_const(&product, &coef);
                ad_restore_order();
                if (std::abs(advec[product][0]) > std::numeric_limits<double>::min()) advec[ov][k] += advec[product][0];
                for (unsigned int idx = 1; idx<adveclen[product] && idx<idx_limit; ++idx) {
                    advec[ov][prdidx[k][idx]] += advec[product][idx];
                }
            }
            else {
                ad_reset_vector(product);
                advec[product][0] = 1;
                adveclen[product] = 1;
                for(unsigned int id=0; id<nv; ++id) {
                    if(rc.at(id)>0) {
                        ad_mult(&product, &power_vv.at(id).at(rc.at(id)), &tmp);
                        std::swap(product, tmp);
//                        ad_copy(&tmp, &product);
                    }
                }
                ad_mult_const(&product, &coef);
                for(unsigned int idx=0; idx<adveclen[product]; ++idx)
                    advec[ov][idx] += advec[product][idx];
            }

            std::fill(rc.begin(), rc.end(), 0);

        }
        else {              //No subscribing to this term.
            advec[ov][i] += advec[iv][i];
        }
    }

    adveclen[ov] = 1;
    for(int i=order_index[gnd+1]-1; i>=0; --i) {
        if(std::abs(advec[ov][i])>std::numeric_limits<double>::min()) {
            adveclen[ov] = i+1;
            break;
        }
    }




    //Release the temporary TPS vectors.
    ad_free(&tmp);
    ad_free(&product);

    for(auto power_v : power_vv) {
        ad_free(&power_v.at(0));
        for(unsigned int i=2; i<(gnd+1); ++i) ad_free(&power_v.at(i));
    }
    delete[] c;
    delete[] bv;
}

/** \brief Subscribe a TPS vector to the specific base of another TPS vector. Save the result into the third TPS vector.
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] iv The TPS vector into which another TPS vector is subscribed.
 * \param[in] base_id The index of the base.
 * \param[in] v The TPS vector to be subscribed.
 * \param[out] ov The TPS vector that saves the result.
 * \return void.
 *
 */
void ad_subscribe(const TVEC iv, unsigned int base_id, const TVEC v, TVEC ov) {
    TNVND* p = base;
    unsigned int* c = new unsigned int[gnv];
    unsigned int* bv = new unsigned int[gnv];

    std::vector<unsigned int> power_v(gnd+1, 0);
    for(auto&& i : power_v) {
        unsigned int idx;
        ad_assign(idx);
        i = idx;
 	}

 	advec[power_v.at(0)][0] = 1;
 	adveclen[power_v.at(0)] = 1;
    ad_copy(&v, &power_v[1]);

    ad_reset(&ov);
    TVEC tmp;

    for (size_t i=0; i<adveclen[iv]; ++i) { //loop over all element.
        if (std::abs(advec[iv][i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        for (size_t j = 0; j < gnv-1; ++j) {
            c[j] = (unsigned int) (*p-*(p+1));
            ++p;
        }
        c[gnv-1] = (unsigned int)*p++ ;
        if (c[base_id]>0) {
            tmp = ad_pow_int_pos(v, power_v, c[base_id], 1);
            double coef = advec[iv][i];
            c[base_id] = 0;

            unsigned int d = 0;
            for (unsigned int i = 0; i < gnv; ++i) {
                d += c[i];
            }                           //d - total order of the term after eliminating the specific base.
            unsigned int order = gnd-d; //Highest order to be multiplied in the TPS vector v^n.
            unsigned int idx_limit = order_index[order+1]; //Limit of the index of the terms to be multiplied in v^n.

            unsigned int k = 0;
            for (unsigned int i = 0; i < gnv; ++i){
                bv[i] = d;
                d -= c[i];
                k += H[gnv-i][bv[i]];
            }                           //k - index of the term after eliminating the specific base.

            if (k>0) {
                if (std::abs(advec[tmp][0]) > std::numeric_limits<double>::min()) advec[ov][k] += coef*advec[tmp][0];
                for (unsigned int idx = 1; idx<adveclen[tmp] && idx<idx_limit; ++idx) {
                    advec[ov][prdidx[k][idx]] += coef*advec[tmp][idx];
                }

            }
            else {
                for(unsigned int idx = 0; idx<adveclen[tmp]; ++idx) {
                    if (std::abs(advec[tmp][idx]) > std::numeric_limits<double>::min()) {
                        advec[ov][idx] += coef*advec[tmp][idx];
                    }
                }
            }
        }
        else {
            advec[ov][i] += advec[iv][i];
        }
    }

    adveclen[ov] = 1;
    for(int i=order_index[gnd+1]-1; i>=0; --i) {
        if(std::abs(advec[ov][i])>std::numeric_limits<double>::min()) {
            adveclen[ov] = i+1;
            break;
        }
    }

    //Release the temporary TPS vectors.
    for(auto i : power_v) {
        ad_free(&i);
 	}
 	delete[] c;
 	delete[] bv;
}

/** \brief Subscribe a number to the specific base of a given TPS vector. Save the result into another TPS vector.
 * This is an alternative function for the original one (ad_subst) in tpsa.cpp.
 * \param[in] iv The TPS vector into which the number is subscribed.
 * \param[in] base_id The index of the base.
 * \param[in] x The number to be subscribed.
 * \param[out] ov The TPS vector that saves the result.
 * \return void.
 *
 */
void ad_subscribe_const(const TVEC iv, unsigned int base_id, double x, TVEC ov) {
    TNVND* p = base;
    unsigned int* c = new unsigned int[gnv];
    double* power_x = new double[gnd+1];
    power_x[0] = 1;
    for(unsigned int i = 1; i<gnd+1; ++i) power_x[i] = power_x[i-1]*x;
    unsigned int* bv = new unsigned int[gnv];

    ad_copy(&iv, &ov);
    if (std::abs(x)>= std::numeric_limits<double>::min()) { // x is not zero.
        for (size_t i=0; i<adveclen[iv]; ++i) { //loop over all element.
            if (std::abs(advec[iv][i]) < std::numeric_limits<double>::min()) {
                p += gnv;
                continue;
            }
            for (size_t j = 0; j < gnv-1; ++j) {
                c[j] = (unsigned int) (*p-*(p+1));
                ++p;
            }
            c[gnv-1] = (unsigned int)*p++ ;
            if (c[base_id]>0) {
                double coef = advec[iv][i]*power_x[c[base_id]];
                c[base_id] = 0;
                advec[ov][i] = 0;

                unsigned int d = 0;
                for (unsigned int i = 0; i < gnv; ++i) {
                    d += c[i];
                }
                //Find the index of the new element eliminating the given base.
                unsigned int k = 0;
                for (unsigned int i = 0; i < gnv; ++i){
                    bv[i] = d;
                    d -= c[i];
                    k += H[gnv-i][bv[i]];
                }

                advec[ov][k] += coef;
            }
        }
    }
    else {      //Consider x as zero! All the respective terms are set to zero.
        for (size_t i=0; i<adveclen[ov]; ++i) { //loop over all elements.
            if (std::abs(advec[ov][i]) < std::numeric_limits<double>::min()) {
                p += gnv;
                continue;
            }
            for (size_t j = 0; j < gnv-1; ++j) {
                c[j] = (unsigned int) (*p-*(p+1));
                ++p;
            }
            c[gnv-1] = (unsigned int)*p++ ;
            if (c[base_id]>0) advec[ov][i] = 0;
        }
    }

    for(unsigned int i=adveclen[ov]-1; i>=0; --i) {
        if(std::abs(advec[ov][i])<std::numeric_limits<double>::min()) adveclen[ov] -= 1;
        else break;
    }
    if (adveclen[ov]==0) adveclen[ov] = 1;

    delete[] c;
    delete[] power_x;
}

/** \brief Summation of  two TPS vectors and save the result to the third one.
 * Calculate ov = idst + jsrc.
 * This is an alternative function for the original one in tpsa.cpp.
 * \param[in] idst TPS vector for summation.
 * \param[in] jsrc TPS vector for summation.
 * \param[out] ov TPS vector that saves the summation result.
 * \return void.
 *
 */
void ad_add(const unsigned int idst, const unsigned int jsrc, unsigned int ov) {
    double *v = advec[idst];
    double *rhsv = advec[jsrc];
    double *resv = advec[ov];

    unsigned int n_copy, n_add;
    n_add = adveclen[idst];
    n_copy = adveclen[jsrc] - n_add;
    if (n_add > adveclen[jsrc]) {
        n_add = adveclen[jsrc];
        n_copy *= -1;
    }

    for (size_t ii = 0; ii<n_add; ++ii) resv[ii] = v[ii] + rhsv[ii];
    if (n_copy>0) {
        if (n_add==adveclen[jsrc]) std::copy(&advec[idst][n_add], &advec[idst][adveclen[idst]], &advec[ov][n_add]);
        else std::copy(&advec[jsrc][n_add], &advec[jsrc][adveclen[jsrc]], &advec[ov][n_add]);
    }
    adveclen[ov] = n_add + n_copy;
}

/** \brief Add a number to a TSP vector and save the result to another TSP vector.
 * Calculate ov = i + r.
 * This is an alternative function for the original one in tpsa.cpp.
 * \param[in] i TSP vector.
 * \param[in] r A number.
 * \param[out] ov TSP vector that saves the result.
 * \return void.
 *
 */
void ad_add_const(const TVEC i, double r, TVEC ov) {
    ad_copy(&i, &ov);
    advec[ov][0] += r;
}

/** \brief Subscribe a TPS vector from another TPS vector and saves the result to the third TPS vector.
 * This is an alternative function for the original one in tpsa.cpp.
 * \param[in] idst The minuend TPS vector.
 * \param[in] jsrc The subtrahend TPS vector.
 * \param[out] ov The TPS vector that saves the results.
 * \return void.
 *
 */
void ad_sub(const unsigned int idst, const unsigned int jsrc, TVEC ov) {
    double *v = advec[idst];
    double *rhsv = advec[jsrc];
    double *resv = advec[ov];

    if (adveclen[idst] == 0 || adveclen[jsrc] == 0) {
        std::cerr << "ERROR: AD sub zero length vector" << std::endl;
        return;
    }

    if (adveclen[idst] < adveclen[jsrc]) {
        // two blocks for overlap part and non-overlap part.
        for (size_t ii = 0; ii < adveclen[idst]; ++ii) resv[ii] = v[ii] - rhsv[ii];
        for (size_t ii = adveclen[idst]; ii < adveclen[jsrc]; ++ii) resv[ii] = -rhsv[ii];
        adveclen[ov] = adveclen[jsrc];
    }else {
        for (size_t ii = 0; ii < adveclen[jsrc]; ++ii) resv[ii] = v[ii] -rhsv[ii];
        std::copy(&v[adveclen[jsrc]], &v[adveclen[idst]], &resv[adveclen[jsrc]]);
        adveclen[ov] = adveclen[idst];
    }
}

/** \brief Return value of a specific element in a TPS vector
 * Take a TPS vector with three bases as an example. Given the order pattern idx = {nx,ny,nz}, the function returns the
 * value of the element (x^nx)*(n^ny)*(z^nz) in vec, the TPS vector.
 * This is an alternative function for the original one in tpsa.cpp.
 * \param vec A TPS vector.
 * \param idx Order pattern of the element.
 * \return Value of the element.
 *
 */
double ad_elem(const TVEC &vec, std::vector<int> &idx) {
    assert(gnv==idx.size()&&"Error in ad_elem: No. of indexes NOT EQUAL to No. of bases!");
    //Find the index of the element
    unsigned int d = 0;
    for (unsigned int i = 0; i < gnv; ++i) {
        d += idx.at(i);
    }

    unsigned int k = 0;
    for (unsigned int i = 0; i < gnv; ++i){
        auto b = d;
        d -= idx.at(i);
        k += H[gnv-i][b];
    }
    return advec[vec][k];
}

// ***** The following functions replace the original ones in tpsa.cpp. *****


/** \brief Reserve space for n TPS vectors.
 * Allocate the memory for n TPS vectors.
 * This function replaces the original one in the tpsa.cpp.
 * \param n Number of TPS vectors.
 * \return void.
 *
 */
void ad_reserve(const unsigned int* n) {
    ad_reserve(*n);
}


/** \brief Assign a memory slot for a TPS vector.
 * Assign a memory slot for a TPS vector with value zero and length one.
 * \param i return a new integer representing the allocated TPS vector.
 * \note This function do not allocate memory. The memory pool was allocated by ad_reserve(n).
 */
void ad_alloc(TVEC* i) {
    ad_assign(*i);
    adveclen[*i] = 1;
}

/** \brief Copy TPS vector
 * This function replaces the original one in tpsa.cpp.
 * \param[in] isrc source
 * \param[out] idst destination
 */
void ad_copy(const TVEC* isrc, const TVEC* idst) {
    unsigned int i = *isrc;
    unsigned int j = *idst;
    if (i == j) return;
    memcpy(advec[j], advec[i], FULL_VEC_LEN*sizeof(double));
    adveclen[j] = adveclen[i];
}

/** \brief Reset a TPS vector to a constant 0
 * This function replaces the original one in tpsa.cpp.
 * \param iv TPS vector.
 */
void ad_reset(const TVEC* iv) {
    ad_reset_vector(*iv);
    adveclen[*iv] = 0;
}

/** \brief Free the memory slot occupied by the TPS vector.
 * Free the memory slot so that it is available for future use by updating adlist and ad_end. The freed memory slot is
 * added to the end of the linked list.
 * This function replaces the original one in tpsa.cpp.
 * \param i TPS vector
 * \note This function does not delete any memory.
 *
 */
void ad_free(const TVEC* i) {
    ad_reset(i);
    adlist.at(*i) = adlist.at(ad_end);
    adlist.at(ad_end) = *i;
    ad_end = *i;
}


/** \brief Internal multiplication of two TPS vectors. Result is saved into the third TPS vector.
 * Calculate idst = ilhs * irhs, where idst should be different from ilhs and irhs.
 * This function replaces the original one in tpsa.cpp.
 * \param[in] ilhs TPS vector
 * \param[in] irhs TPS vector
 * \param[out] idst TPS vector that saves the result.
 * \return void
 *
 */
void ad_mult(const TVEC* ilhs, const TVEC* irhs, TVEC* idst) {
    unsigned int lhs = *ilhs;
    unsigned int rhs = *irhs;
    unsigned int dst = *idst;

    memset(advec[dst], 0, FULL_VEC_LEN*sizeof(double));
//    memset(advec[dst], 0, adveclen[dst]*sizeof(double));   //This will cause a bug of extra terms in some cases!

    adveclen[dst] = adveclen[lhs];

    advec[dst][0] = advec[lhs][0] * advec[rhs][0];
    if (advec[lhs][0]!=0) {
        for (size_t i = 1; i < std::min(adveclen[rhs], order_index[gnd+1]); ++i) {
//        for (size_t i = 1; i <adveclen[rhs]; ++i) {
            advec[dst][i] += advec[lhs][0]*advec[rhs][i];
        }
    }

    if (advec[rhs][0]!=0) {
        for (size_t i = 1; i < std::min(adveclen[lhs], order_index[gnd+1]); ++i) {
//        for (size_t i = 1; i < adveclen[lhs]; ++i) {
            advec[dst][i] += advec[lhs][i]*advec[rhs][0];
        }
    }

    // It is symmetric, but not when lhs rhs have different size, .....
    TNVND ord = 1;
    size_t L = std::max(adveclen[lhs], adveclen[rhs]);

    for (size_t i = 1; i < std::min(adveclen[lhs], order_index[gnd]); ++i) {
        if (order_index[ord+1] <= i) ++ord;
        unsigned int M = order_index[gnd-ord+1];
        if (M > adveclen[rhs]) M = adveclen[rhs];
        if(advec[lhs][i]!=0) {
            for (size_t j = 1; j < M; ++j) {
                advec[dst][prdidx[i][j]] += advec[lhs][i]*advec[rhs][j];
            }
        }

        if (prdidx[i][M-1] >= L) {
            L = prdidx[i][M-1] + 1;
        }
    }

    if (L > FULL_VEC_LEN) L = FULL_VEC_LEN;
    adveclen[dst] = L;

    while(std::abs(advec[dst][adveclen[dst]-1]) < std::numeric_limits<double>::min()  && adveclen[dst] > 1) --adveclen[dst];

}


/** \brief A TPS vector divided by a number.
 * Calculate iv = iv / c.
 * This function replaces the original one in tpsa.cpp.
 * \param iv A TPS vector.
 * \param c a real number.
 * \return void.
 *
 */
void ad_div_c(const TVEC* iv, const double* c) {
    if (std::abs(*c) < std::numeric_limits<double>::min()) {
        std::cerr << "ERROR: divide a two small number! " << *c << std::endl;
        std::exit(-1);
    }

    double c_inv = 1 / *c;
    for (size_t i = 0; i < std::min(adveclen[*iv], order_index[gnd+1]); ++i)
//    for (size_t i = 0; i < adveclen[*iv]; ++i)
        advec[*iv][i] *= c_inv;
}

/** \brief A real number is divided by a TPS vector. Result is saved to another TPS vector.
 * Calculate ivret = c / iv.
 * This function replaces the original one in tpsa.cpp.
 * \param[in] iv A TPS vector.
 * \param[in] c A real number.
 * \param[out] ivret A TPS vector that saves the result.
 * \return void.
 *
 */
void ad_c_div(const TVEC* iv, const double* c, TVEC* ivret) {
    // If c is zero, set ivret zero and return.
    if (std::abs(*c) < std::numeric_limits<double>::min()) {
        ad_reset(ivret);
        return;
    }

    TVEC ipn, ip, itmp;
    ad_alloc(&ipn);
    ad_alloc(&ip);
    ad_alloc(&itmp);

    TVEC iret = *ivret;

    double *pn = advec[ipn];
    double *p = advec[ip];
    double *ret = advec[iret];
    double *v = advec[*iv];

    ad_copy(iv, &ip);
    ad_copy(iv, &ipn);

    double x0 = v[0];

    double inv_x0 = 1/x0;
    inv_x0 *= -1;
    for (size_t i = 1; i < adveclen[ip]; ++i) {
        p[i] *= inv_x0;
    }
    pn[0] = p[0] = 0;
    inv_x0 *= -1;
    memcpy(pn, p, adveclen[ip]*sizeof(double));

    ret[0] = 1;
    adveclen[iret] = 1;

    ad_add(&iret, &ip);

    for (TNVND nd = 2; nd < gnd+1; ++nd) {
        ad_mult(&ipn, &ip, &itmp);
        ad_copy(&itmp, &ipn);
        for (size_t i = 0; i < adveclen[ipn]; ++i)
            ret[i] += pn[i];
    }
    adveclen[iret] = FULL_VEC_LEN;
    double ret_coef = inv_x0*(*c);

    for (size_t i = 0; i < adveclen[iret]; ++i) {
        ret[i] *= ret_coef;
    }
    ad_free(&itmp);
    ad_free(&ip);
    ad_free(&ipn);
}

/** \brief Set the small TPS coefficients as 0 if they are less than eps.
 * This function replaces the original one in tpsa.cpp.
 * \param[in] iv A TPS vector.
 * \param[in] eps A real number.
 * \return void.
 */
void ad_clean(const TVEC& iv, const double eps)
{
    unsigned int N = 0;
    for(size_t i = 0; i < adveclen[iv]; ++i) {
        if (abs(advec[iv][i]) < abs(eps)) advec[iv][i] = 0;
        else N = i;
    }
    if (adveclen[iv] > N+1) {
        adveclen[iv] = N+1;
    }
}
