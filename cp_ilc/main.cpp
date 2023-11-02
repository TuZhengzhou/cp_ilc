// #define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include "structs.tcc"
#include "sum.tcc"
#include "eq.tcc"
#include "prod.tcc"
#include "shift.tcc"
#include "same_prod.tcc"
#include "permutation.tcc"
#include "comEq.tcc"
#include "ILC.tcc"
// #include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include "../libsnark/depends/libff/libff/common/profiling.hpp"
#include "../libsnark/depends/libff/libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "../libsnark/depends/libff/libff/algebra/curves/bn128/bn128_pp.hpp"

int main() {
    libff::start_profiling();
    // libff::mnt6_pp::init_public_params();
    // typedef libff::mnt6_pp ppT;

    // std::cout << libff::bn128_pp::G1_type::G1_one << std::endl;
    // std::cout << "***" << libff::mnt6_pp::G1_type::G1_one << "***" << std::endl;
    // std::cout << libff::Fr<ppT>::one() << std::endl;
    // return 0;

    libff::bn128_pp::init_public_params();
    typedef libff::bn128_pp ppT;
    typedef ppT::Fp_type FieldT;
    typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT> > CommitT;

    // eq_test<FieldT>();
    // sum_test<FieldT>();
    // prod_test<FieldT>();
    // shift_test<FieldT>();
    // same_prod_test<FieldT>();
    // permutation_test<FieldT>();
    // comEq_test<FieldT, ppT>();

    libff::inhibit_profiling_counters = true;
    typedef libff::G1<ppT> HType1;
    typedef libff::G2<ppT> HType2;

    // comEq_tests<FieldT, ppT, HType1>();
    // comEq_tests<FieldT, ppT, HType2>();
    // comEq_filling_tests<FieldT, ppT, HType1>();
    // comEq_filling_tests<FieldT, ppT, HType2>();

    libff::inhibit_profiling_counters = true;
    ILC_test_compare_1<FieldT, ppT, HType1>();
    ILC_test_compare_2<FieldT, ppT, HType1>();
    ILC_test_compare_3<FieldT, ppT, HType1>();

    // std::cout << "FieldT::size_in_bits() = " << FieldT::size_in_bits() << std::endl;
    // std::cout << "libff::G1<ppT>::size_in_bits() = " << libff::G1<ppT>::size_in_bits() << std::endl;
    // std::cout << "libff::G2<ppT>::size_in_bits() = " << libff::G2<ppT>::size_in_bits() << std::endl;
    // // Create an Fr element
    // FieldT fr = FieldT::random_element();
    // // std::string msg = fr_to_string<FieldT>(fr);  // 待哈希的字符串
    // std::string msg = frs_to_string<FieldT>(row_vector<FieldT>::random(12).get_all_items());  // 待哈希的字符串
    // std::string hash_hex = sha256_to_hex_string(msg);
    // FieldT fr_x = sha256_hex_digest_to_Fr<FieldT>(hash_hex);
    // std::cout << "Fr_x: " << fr_x << std::endl;

    return 0;
}