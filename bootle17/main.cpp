// #define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include "structs.tcc"
#include "sum.tcc"
#include "eq.tcc"
#include "prod.tcc"
#include "shift.tcc"
#include "same_prod.tcc"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/common/profiling.hpp>

int main() {
    typedef libff::Fr<libff::mnt6_pp> FieldT;
    libff::start_profiling();
    libff::mnt6_pp::init_public_params();
    // size_t_test();

    // std::cout << "FieldT::one() = " << FieldT::one() << std::endl;

    // group_test<FieldT>();

    // field_test<FieldT>();

    // matrix_test<FieldT>();

    eq_test<FieldT>();

    sum_test<FieldT>();

    prod_test<FieldT>();

    shift_test<FieldT>();

    same_prod_test<FieldT>();

    // Create an Fr element
    FieldT fr = FieldT::random_element();
    // std::string msg = fr_to_string<FieldT>(fr);  // 待哈希的字符串
    std::string msg = frs_to_string<FieldT>(row_vector<FieldT>::random(12).get_all_items());  // 待哈希的字符串
    std::string hash_hex = sha256_to_hex_string(msg);
    FieldT fr_x = sha256_hex_digest_to_Fr<FieldT>(hash_hex);
    std::cout << "Fr_x: " << fr_x << std::endl;

    return 0;
}