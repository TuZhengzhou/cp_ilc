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

#include "../libsnark/depends/libff/libff/common/profiling.hpp"
#include "../libsnark/depends/libff/libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "../libsnark/depends/libff/libff/algebra/curves/bn128/bn128_pp.hpp"

int main() {
    libff::start_profiling();




    libff::bn128_pp::init_public_params();
    typedef libff::bn128_pp ppT;
    typedef ppT::Fp_type FieldT;
    typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT> > CommitT;

    std::cout << "FieldT::size_in_bits(): " << FieldT::size_in_bits() << std::endl;
    std::cout << "FieldT::num_bits: " << FieldT::num_bits << std::endl;


    typedef libff::G1<ppT> HType1;
    typedef libff::G2<ppT> HType2;



    libff::inhibit_profiling_info = true;
    ILC_test_compare_1<FieldT, ppT, HType1>();

    // ILC_test_compare_2<FieldT, ppT, HType1>();

    // ILC_test_compare_3<FieldT, ppT, HType1>();

    libff::print_cumulative_times();

    libff::print_cumulative_op_counts();

    return 0;
}