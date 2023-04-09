#define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <../libsnark/libsnark/relations/arithmetic_programs/qap/qap.hpp>
#include <../libsnark/libsnark/reductions/r1cs_to_qap/r1cs_to_qap.hpp>


namespace libsnark {
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_single_var_polynomial(std::vector<int32_t> coefs, size_t degree) {
    // 形如 y = c0 * 1 + c1 * x + c2 * x^2 + ... c_(degree) * x ^ (degree) 形式的多项式
    // coefs 表示多项式的系数
    // degree 表示多项式的阶数
    std::cout << coefs[0] << std::endl;

    size_t num_constraints = degree;

    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = 1;
    cs.auxiliary_input_size = num_constraints;
    
    r1cs_variable_assignment<FieldT> full_variable_assignment;

    std::cout << coefs[1] << std::endl;
    FieldT x = FieldT::random_element();
    std::cout << coefs[2] << std::endl;
    FieldT a = x;
    FieldT b = x;
    full_variable_assignment.push_back(x);
    for (size_t i = 0; i < num_constraints-1; ++i)
    {
        // std::cout << "i = " << i << std::endl;
        linear_combination<FieldT> A, B, C;
        A.add_term(1, 1);
        B.add_term(i+1, 1);
        C.add_term(i+2, 1);
        FieldT tmp = a * b;
        full_variable_assignment.push_back(tmp);
        b = tmp;

        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        // std::cout << "cs.constraints[" << i << "] = " << std::endl;
        // std::cout << "A = " << A << std::endl;
        // std::cout << "B = " << B << std::endl;
        // std::cout << "C = " << C << std::endl;
    }

    linear_combination<FieldT> A, B, C;
    FieldT y = FieldT::zero();
    for (size_t i = 0; i <= degree; i++) {
        A.add_term(i, coefs[i]);
        if ( i == 0 ) {
            y = y + FieldT::one() * coefs[0];
        } else {
            y = y + full_variable_assignment[i-1] * coefs[i];
        }
    }
    B.add_term(0, 1);
    C.add_term(degree + 1, 1);
    full_variable_assignment.push_back(y);
    cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    // std::cout << "cs.constraints[" << 2 << "] = " << cs.constraints[2] << std::endl;
    std::cout << cs.constraints << std::endl;

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + 1);
    r1cs_auxiliary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + 1, full_variable_assignment.end());

    /* sanity checks */
    // std::cout << "cs.num_variables() = " << cs.num_variables() << std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= 1);
    assert(cs.num_inputs() == 1);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}
}

using namespace libsnark;

int main() {
    // libff::start_profiling();

    libff::mnt6_pp::init_public_params();

    std::cout << "I an here" << std::endl;

    /*
        x * x = x ^ 2
        x * x^2 = x ^ 3
        15 + 3 * x + 4 * x^2 + 7 * x^3 = y
    */
    // std::vector<int32_t> coefs = {15, 3, 4, 7};  // y = 15 + 3 * x + 4 * x ^2 + 7 * x ^ 3
    std::vector<int32_t> coefs = {15, 3, 4, 7, 139};  // y = 15 + 3 * x + 4 * x^2 + 7 * x^3 + 139 * x^4
    size_t degree = coefs.size()-1;

    r1cs_example<libff::Fr<libff::mnt6_pp> > example;
    example = generate_r1cs_example_from_single_var_polynomial<libff::Fr<libff::mnt6_pp> >(coefs, degree);

    const libff::Fr<libff::mnt6_pp> t = libff::Fr<libff::mnt6_pp>::random_element(),
    d1 = libff::Fr<libff::mnt6_pp>::random_element(),
    d2 = libff::Fr<libff::mnt6_pp>::random_element(),
    d3 = libff::Fr<libff::mnt6_pp>::random_element();
    qap_instance<libff::Fr<libff::mnt6_pp>> qap_inst_1 = r1cs_to_qap_instance_map(example.constraint_system);
    qap_instance_evaluation<libff::Fr<libff::mnt6_pp>> qap_inst_2 = r1cs_to_qap_instance_map_with_evaluation(example.constraint_system, t);
    qap_witness<libff::Fr<libff::mnt6_pp>> qap_wit = r1cs_to_qap_witness_map(example.constraint_system, example.primary_input, example.auxiliary_input, d1, d2, d3);
    assert(qap_inst_1.is_satisfied(qap_wit));
    assert(qap_inst_2.is_satisfied(qap_wit));

    // std::cout << qap_inst_1.A_in_Lagrange_basis[0] << std::endl;
    std::cout << qap_inst_1.A_in_Lagrange_basis << std::endl;
    std::cout << qap_inst_1.B_in_Lagrange_basis << std::endl;
    std::cout << qap_inst_1.C_in_Lagrange_basis << std::endl;
    return 0;
}