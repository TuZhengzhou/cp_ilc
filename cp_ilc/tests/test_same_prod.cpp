// #define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include "structs.tcc"
#include "same_prod.tcc"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

int main() {
    libff::start_profiling();
    libff::mnt6_pp::init_public_params();
    typedef libff::mnt6_pp ppT;

    typedef ppT::Fp_type FieldT;

    same_prod_test<FieldT>();

    return 0;
}