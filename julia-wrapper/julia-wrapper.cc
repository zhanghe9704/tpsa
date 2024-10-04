#include <jlcxx/jlcxx.hpp>
#include "da.h"

Base& get_base() {return da;}
std::complex<DAVector> sin(const std::complex<DAVector>& c) {return sin(c);}
std::complex<DAVector> cos(const std::complex<DAVector>& c) {return cos(c);}
std::complex<DAVector> tan(const std::complex<DAVector>& c) {return tan(c);}

// Define the module and expose DAVector class and non-member functions separately
JLCXX_MODULE juliaTPSA(jlcxx::Module& mod) {
    // Register the DAVector class
    mod.add_type<DAVector>("DAVector")
        // Constructors
        .constructor<>()
        .constructor<const DAVector&>()
        // .constructor<DAVector&&>()
        .constructor<double>()
        .constructor<int>()
        .constructor<std::vector<double>&>()

        // Member Functions
        .method("print", &DAVector::print)
        .method("con", &DAVector::con)
        .method("length", &DAVector::length)
        .method("n_element", &DAVector::n_element)
        .method("element", static_cast<double (DAVector::*)(int)>(&DAVector::element))
        .method("element", static_cast<double (DAVector::*)(std::vector<int>)>(&DAVector::element))
        .method("element", static_cast<void (DAVector::*)(unsigned int, std::vector<unsigned int>&, double&) const>(&DAVector::element))
        .method("derivative", static_cast<double (DAVector::*)(std::vector<int>)>(&DAVector::derivative))
        .method("norm", &DAVector::norm)
        .method("weighted_norm", &DAVector::weighted_norm)
        .method("reset", &DAVector::reset)
        .method("reset", &DAVector::reset_const)
        .method("clean", static_cast<void (DAVector::*)()>(&DAVector::clean))
        .method("clean", static_cast<void (DAVector::*)(const double)>(&DAVector::clean))
        .method("iszero", static_cast<bool (DAVector::*)() const>(&DAVector::iszero))
        .method("iszero", static_cast<bool (DAVector::*)(double) const>(&DAVector::iszero))
        .method("clear", &DAVector::clear)

        // Assignment operators
        .method("assign", static_cast<DAVector& (DAVector::*)(const DAVector&)>(&DAVector::operator=))
        // .method("assign_move", static_cast<DAVector& (DAVector::*)(DAVector&&)>(&DAVector::operator=))
        .method("assign", static_cast<DAVector& (DAVector::*)(double)>(&DAVector::operator=))
        .method("assign", static_cast<DAVector& (DAVector::*)(int)>(&DAVector::operator=))

        // Compound assignment operators (+=)
        .method("add_assign", static_cast<DAVector& (DAVector::*)(const DAVector&)>(&DAVector::operator+=))
        // .method("+=_move", static_cast<DAVector& (DAVector::*)(DAVector&&)>(&DAVector::operator+=))
        .method("add_assign", static_cast<DAVector& (DAVector::*)(double)>(&DAVector::operator+=))
        .method("add_assign", static_cast<DAVector& (DAVector::*)(int)>(&DAVector::operator+=))

        // Compound assignment operators (-=)
        .method("sub_assign", static_cast<DAVector& (DAVector::*)(const DAVector&)>(&DAVector::operator-=))
        // .method("-=_move", static_cast<DAVector& (DAVector::*)(DAVector&&)>(&DAVector::operator-=))
        .method("sub_assign", static_cast<DAVector& (DAVector::*)(double)>(&DAVector::operator-=))
        .method("sub_assign", static_cast<DAVector& (DAVector::*)(int)>(&DAVector::operator-=))

        // Compound assignment operators (*=)
        .method("mul_assign", static_cast<DAVector& (DAVector::*)(const DAVector&)>(&DAVector::operator*=))
        // .method("*=_move", static_cast<DAVector& (DAVector::*)(DAVector&&)>(&DAVector::operator*=))
        .method("mul_assign", static_cast<DAVector& (DAVector::*)(double)>(&DAVector::operator*=))
        .method("mul_assign", static_cast<DAVector& (DAVector::*)(int)>(&DAVector::operator*=))

        // Compound assignment operators (/=)
        .method("div_assign", static_cast<DAVector& (DAVector::*)(const DAVector&)>(&DAVector::operator/=))
        // .method("/=_move", static_cast<DAVector& (DAVector::*)(DAVector&&)>(&DAVector::operator/=))
        .method("div_assign", static_cast<DAVector& (DAVector::*)(double)>(&DAVector::operator/=))
        .method("div_assign", static_cast<DAVector& (DAVector::*)(int)>(&DAVector::operator/=));

    // Initialization and cleanup functions
    mod.method("da_init", &da_init);
    mod.method("da_clear", &da_clear);
    mod.method("da_change_order", &da_change_order);
    mod.method("da_restore_order", &da_restore_order);
    mod.method("da_count", &da_count);
    mod.method("da_remain", &da_remain);
    mod.method("da_full_length", &da_full_length);

    // Functions returning a std::vector
    mod.method("da_element_orders", &da_element_orders);
    
    // Mathematical functions for differentiation and integration
    mod.method("da_der", static_cast<void(*)(const DAVector&, unsigned int, DAVector&)>(&da_der));
    mod.method("da_int", static_cast<void(*)(const DAVector&, unsigned int, DAVector&)>(&da_int));

    mod.method("da_der", static_cast<DAVector(*)(const DAVector&, unsigned int)>(&da_der));
    mod.method("da_int", static_cast<DAVector(*)(const DAVector&, unsigned int)>(&da_int));

    // Substitution functions
    mod.method("da_substitute_const", &da_substitute_const);
    mod.method("da_substitute", static_cast<void(*)(const DAVector&, unsigned int, const DAVector&, DAVector&)>(&da_substitute));
    mod.method("da_substitute", static_cast<void(*)(const DAVector&, std::vector<unsigned int>&, std::vector<DAVector>&, DAVector&)>(&da_substitute));
    mod.method("da_substitute", static_cast<void(*)(std::vector<DAVector>&, std::vector<unsigned int>&, std::vector<DAVector>&, std::vector<DAVector>&)>(&da_substitute));

    // Composition functions
    mod.method("da_composition", static_cast<void(*)(std::vector<DAVector>&, std::vector<DAVector>&, std::vector<DAVector>&)>(&da_composition));
    mod.method("da_composition", static_cast<void(*)(std::vector<DAVector>&, std::vector<double>&, std::vector<double>&)>(&da_composition));
    mod.method("da_composition", static_cast<void(*)(std::vector<DAVector>&, std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>(&da_composition));

    // Complex composition functions
    mod.method("cd_composition", static_cast<void(*)(std::vector<DAVector>&, std::vector<std::complex<DAVector>>&, std::vector<std::complex<DAVector>>&)>(&cd_composition));
    mod.method("cd_composition", static_cast<void(*)(std::vector<std::complex<DAVector>>&, std::vector<std::complex<DAVector>>&, std::vector<std::complex<DAVector>>&)>(&cd_composition));
    mod.method("cd_composition", static_cast<void(*)(std::vector<std::complex<DAVector>>&, std::vector<DAVector>&, std::vector<std::complex<DAVector>>&)>(&cd_composition));

    // Set epsilon
    mod.method("da_set_eps", &da_set_eps);

    // Overloaded + operators
    mod.method("add_dvector_double", static_cast<DAVector (*)(const DAVector&, double)>(&operator+));
    mod.method("add_double_dvector", static_cast<DAVector (*)(double, const DAVector&)>(&operator+));
    mod.method("add_dvector_dvector", static_cast<DAVector (*)(const DAVector&, const DAVector&)>(&operator+));

    // Overloaded * operators
    mod.method("mul_dvector_double", static_cast<DAVector (*)(const DAVector&, double)>(&operator*));
    mod.method("mul_double_dvector", static_cast<DAVector (*)(double, const DAVector&)>(&operator*));
    mod.method("mul_dvector_dvector", static_cast<DAVector (*)(const DAVector&, const DAVector&)>(&operator*));

    // Overloaded - operators
    mod.method("sub_dvector_double", static_cast<DAVector (*)(const DAVector&, double)>(&operator-));
    mod.method("sub_double_dvector", static_cast<DAVector (*)(double, const DAVector&)>(&operator-));
    mod.method("sub_dvector_dvector", static_cast<DAVector (*)(const DAVector&, const DAVector&)>(&operator-));

    // Overloaded / operators
    mod.method("div_dvector_double", static_cast<DAVector (*)(const DAVector&, double)>(&operator/));
    mod.method("div_double_dvector", static_cast<DAVector (*)(double, const DAVector&)>(&operator/));
    mod.method("div_dvector_dvector", static_cast<DAVector (*)(const DAVector&, const DAVector&)>(&operator/));

    // Unary + and - operators
    mod.method("unary_plus", static_cast<DAVector (*)(const DAVector&)>(&operator+));
    mod.method("unary_minus", static_cast<DAVector (*)(const DAVector&)>(&operator-));

    // Equality operator
    mod.method("eq_dvector_dvector", static_cast<bool (*)(const DAVector&, const DAVector&)>(&operator==));

    // Overloaded + with DAVector and std::complex<double>
    mod.method("add_dvector_complex", static_cast<std::complex<DAVector> (*)(const DAVector&, std::complex<double>)>(&operator+));
    mod.method("add_complex_dvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const DAVector&)>(&operator+));
    
    // Overloaded - with DAVector and std::complex<double>
    mod.method("sub_dvector_complex", static_cast<std::complex<DAVector> (*)(const DAVector&, std::complex<double>)>(&operator-));
    mod.method("sub_complex_dvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const DAVector&)>(&operator-));
    
    // Overloaded * with DAVector and std::complex<double>
    mod.method("mul_dvector_complex", static_cast<std::complex<DAVector> (*)(const DAVector&, std::complex<double>)>(&operator*));
    mod.method("mul_complex_dvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const DAVector&)>(&operator*));
    
    // Overloaded / with DAVector and std::complex<double>
    mod.method("div_dvector_complex", static_cast<std::complex<DAVector> (*)(const DAVector&, std::complex<double>)>(&operator/));
    mod.method("div_complex_dvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const DAVector&)>(&operator/));

    // Overloaded operators for std::complex<DAVector> and double
    mod.method("add_cdvector_double", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, double)>(&operator+));
    mod.method("add_double_cdvector", static_cast<std::complex<DAVector> (*)(double, const std::complex<DAVector>&)>(&operator+));
    
    mod.method("sub_cdvector_double", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, double)>(&operator-));
    mod.method("sub_double_cdvector", static_cast<std::complex<DAVector> (*)(double, const std::complex<DAVector>&)>(&operator-));
    
    mod.method("mul_cdvector_double", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, double)>(&operator*));
    mod.method("mul_double_cdvector", static_cast<std::complex<DAVector> (*)(double, const std::complex<DAVector>&)>(&operator*));
    
    mod.method("div_cdvector_double", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, double)>(&operator/));
    mod.method("div_double_cdvector", static_cast<std::complex<DAVector> (*)(double, const std::complex<DAVector>&)>(&operator/));

    // Overloaded operators for std::complex<DAVector> and std::complex<double>
    mod.method("add_cdvector_cdouble", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, std::complex<double>)>(&operator+));
    mod.method("add_cdouble_cdvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const std::complex<DAVector>&)>(&operator+));
    
    mod.method("sub_cdvector_cdouble", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, std::complex<double>)>(&operator-));
    mod.method("sub_cdouble_cdvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const std::complex<DAVector>&)>(&operator-));
    
    mod.method("mul_cdvector_cdouble", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, std::complex<double>)>(&operator*));
    mod.method("mul_cdouble_cdvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const std::complex<DAVector>&)>(&operator*));
    
    mod.method("div_cdvector_cdouble", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, std::complex<double>)>(&operator/));
    mod.method("div_cdouble_cdvector", static_cast<std::complex<DAVector> (*)(std::complex<double>, const std::complex<DAVector>&)>(&operator/));

    // Overloaded +, -, *, / operators for std::complex<DAVector> objects
    mod.method("add_cdvector_cdvector", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const std::complex<DAVector>&)>(&operator+));
    mod.method("sub_cdvector_cdvector", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const std::complex<DAVector>&)>(&operator-));
    mod.method("mul_cdvector_cdvector", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const std::complex<DAVector>&)>(&operator*));
    mod.method("div_cdvector_cdvector", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const std::complex<DAVector>&)>(&operator/));



    // Expose Non-Member Mathematical Functions separately
    mod.method("base", &get_base);
    mod.method("sqrt", static_cast<DAVector (*)(const DAVector&)>(&sqrt));
    mod.method("exp", static_cast<DAVector (*)(const DAVector&)>(&exp));
    mod.method("log", static_cast<DAVector (*)(const DAVector&)>(&log));
    mod.method("sin", static_cast<DAVector (*)(const DAVector&)>(&sin));
    mod.method("cos", static_cast<DAVector (*)(const DAVector&)>(&cos));
    mod.method("tan", static_cast<DAVector (*)(const DAVector&)>(&tan));
    mod.method("asin", static_cast<DAVector (*)(const DAVector&)>(&asin));
    mod.method("acos", static_cast<DAVector (*)(const DAVector&)>(&acos));
    mod.method("atan", static_cast<DAVector (*)(const DAVector&)>(&atan));
    mod.method("sinh", static_cast<DAVector (*)(const DAVector&)>(&sinh));
    mod.method("cosh", static_cast<DAVector (*)(const DAVector&)>(&cosh));
    mod.method("tanh", static_cast<DAVector (*)(const DAVector&)>(&tanh));
    mod.method("asinh", static_cast<DAVector (*)(const DAVector&)>(&asinh));
    mod.method("acosh", static_cast<DAVector (*)(const DAVector&)>(&acosh));
    mod.method("atanh", static_cast<DAVector (*)(const DAVector&)>(&atanh));

    // Power functions
    mod.method("pow", static_cast<DAVector (*)(const DAVector&, const int)>(&pow));
    mod.method("pow", static_cast<DAVector (*)(const DAVector&, const double)>(&pow));
    mod.method("pow", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const int)>(&pow));
    mod.method("pow", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&, const double)>(&pow));

    // Absolute value
    mod.method("abs", static_cast<double (*)(const DAVector&)>(&abs));
    mod.method("abs", static_cast<double (*)(const std::complex<DAVector>&)>(&abs));

    // Error function and atan2
    mod.method("erf", static_cast<DAVector (*)(const DAVector&)>(&erf));
    mod.method("atan2", static_cast<DAVector (*)(const DAVector&, const DAVector&)>(&atan2));


    // Complex DAVector math functions
    mod.method("exp", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&exp));
    mod.method("sqrt", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&sqrt));
    mod.method("log", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&log));
    mod.method("sin", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&sin));
    mod.method("cos", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&cos));
    mod.method("tan", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&tan));
    mod.method("asin", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&asin));
    mod.method("acos", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&acos));
    mod.method("atan", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&atan));
    mod.method("asinh", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&asinh));
    mod.method("acosh", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&acosh));
    mod.method("atanh", static_cast<std::complex<DAVector> (*)(const std::complex<DAVector>&)>(&atanh));


    // Get real and imaginary parts of std::complex<DAVector>
    mod.method("get_real", static_cast<DAVector& (*)(std::complex<DAVector>&)>(&get_real));
    mod.method("get_imag", static_cast<DAVector& (*)(std::complex<DAVector>&)>(&get_imag));
    mod.method("get_real_const", static_cast<const DAVector& (*)(const std::complex<DAVector>&)>(&get_real));
    mod.method("get_imag_const", static_cast<const DAVector& (*)(const std::complex<DAVector>&)>(&get_imag));

    // Copy functions for std::complex<DAVector>
    mod.method("cd_copy_cdvector", static_cast<void (*)(std::complex<DAVector>&, std::complex<DAVector>&)>(&cd_copy));
    mod.method("cd_copy_cdouble", static_cast<void (*)(std::complex<double>, std::complex<DAVector>&)>(&cd_copy));
    mod.method("cd_copy_double", static_cast<void (*)(double, std::complex<DAVector>&)>(&cd_copy));

    // String trimming function
    mod.method("trim_whitespace", &trim_whitespace);

    // Reading DA and CD from files
    mod.method("read_da_from_file", &read_da_from_file);
    mod.method("read_cd_from_file", &read_cd_from_file);

    // DA vector division
    mod.method("devide_by_element", &devide_by_element);

    // DA and CD vector comparison functions
    mod.method("compare_da_vectors", &compare_da_vectors);
    mod.method("compare_da_with_file", &compare_da_with_file);
    mod.method("compare_cd_vectors", &compare_cd_vectors);
    mod.method("compare_cd_with_file", &compare_cd_with_file);


    // Expose the '<<' operator for DAVector
    mod.method("ostream_da", [](std::ostream& os, const DAVector& da_vector) -> std::ostream& {
        return os << da_vector;
    });

    // Expose the '<<' operator for std::complex<DAVector>
    mod.method("ostream_cd", [](std::ostream& os, const std::complex<DAVector>& cd_vector) -> std::ostream& {
        return os << cd_vector;
    });
}
