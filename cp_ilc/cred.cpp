#include<iostream>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "libff/algebra/curves/bn128/bn128_pp.hpp"
#include "libff/algebra/curves/public_params.hpp"
using namespace std;

typedef libff::mnt6_pp ppT;
int main() {
    libff::start_profiling();
    libff::mnt6_pp::init_public_params();
    libff::bn128_pp::init_public_params();
    std::cout << libff::bn128_pp::G1_type::G1_one << std::endl;
    std::cout << "***" << libff::mnt6_pp::G1_type::G1_one << "***" << std::endl;
    // std::cout << "***" << libff::bn128_pp::G1_type::random_element() << "***" << std::endl;
    // std::cout 
    std::cout << libff::Fr<ppT>::one() << std::endl;
    // // cout << libff::G1<ppT>::G1_one << endl;
    // cout << libff::Fr<ppT>::one() << endl;
    cout << "output OK:>" << endl;
}

