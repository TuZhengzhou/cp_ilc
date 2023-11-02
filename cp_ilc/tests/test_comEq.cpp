// #define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include "structs.tcc"
#include "comEq.tcc"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

int main() {
    libff::start_profiling();
    libff::mnt6_pp::init_public_params();
    typedef libff::mnt6_pp ppT;

    typedef ppT::Fp_type FieldT;
    typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT> > CommitT;

    comEq_test<FieldT, ppT>();

    return 0;
}