#ifndef BOOTLE_PERMUTATION_TCC
#define BOOTLE_PERMUTATION_TCC
#include <iostream>
#include <assert.h>
#include <vector>
#include "permutation.hpp"

template <typename FieldT>
pp_perm<FieldT>::pp_perm(size_t k, size_t mu, size_t n, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, \
                            const permutation<tuple_dim2_t, 2>& perm) {
    this->k_ = k;
    this->mu_ = mu;
    this->m_ = std::pow(2ul, this->mu_);
    this->n_ = n;
    this->N_ = this->k_ * this->m_ * this->n_;
    printf("k = \033[32m%ld\033[37m, mu = \033[32m%ld\033[37m, m = \033[32m%ld\033[37m, n = \033[32m%ld\033[37m, N = \033[32m%ld\033[37m\n", this->k_, this->mu_, this->m_, this->n_, this->N_);

    this->A_ = row_vector_matrix<FieldT>(A);
    this->B_ = row_vector_matrix<FieldT>(B);
    this->PI_ = permutation<tuple_dim2_t, 2>(perm);
    
    this->V_ = row_vector_matrix<FieldT>::linear_grow(this->m_ * this->n_, this->k_, 1, 1);
    this->V_PI_ = apply_permutation<FieldT, 2>(this->V_, perm);
    this->J_ = row_vector_matrix<FieldT>::all_one(this->m_ * this->n_, this->k_);
    // this->V_.print();
    // this->J_.print();
}

    
template <typename FieldT>
bool pp_perm<FieldT>::prove(const FieldT& x, const FieldT& y, const std::map<std::string, FieldT> &mid_challenges) {

    this->U_    = this->V_ * y - this->J_ * x;
    this->U_PI_ = this->V_PI_ * y - this->J_ * x;
    this->A_PRIME_ = this->A_ + this->U_;
    this->B_PRIME_ = this->B_ + this->U_PI_;

    this->eq_U_ = pp_eq<FieldT>(this->U_.get_column_num(), this->U_);
    this->eq_U_PI_ = pp_eq<FieldT>(this->U_PI_.get_column_num(), this->U_PI_);

    this->sum_A_PRIME_ = pp_sum<FieldT>(this->U_.get_column_num(), this->A_, this->U_, this->A_PRIME_);
    this->sum_B_PRIME_ = pp_sum<FieldT>(this->U_.get_column_num(), this->B_, this->U_PI_, this->B_PRIME_);

    this->same_prod_AP_BP_ = pp_same_prod<FieldT>(this->mu_, this->n_, this->k_, this->A_PRIME_, this->B_PRIME_);
    this->same_prod_AP_BP_.prove(mid_challenges);
    return true;
}

template <typename FieldT>
bool pp_perm<FieldT>::verify(const FieldT& x, const FieldT& y, const std::map<std::string, FieldT> &mid_challenges, const std::map<std::string, FieldT> &challenges, const bool output) const {
    row_vector_matrix<FieldT> U_V_    = this->V_ * y - this->J_ * x;
    row_vector_matrix<FieldT> U_PI_V_ = this->V_PI_ * y - this->J_ * x;

    bool check_result = true;

    FieldT x_eq_U_ = FieldT::random_element();
    FieldT x_eq_U_PI_ = FieldT::random_element();
    FieldT x_sum_A_PRIME = FieldT::random_element();
    FieldT x_sum_B_PRIME = FieldT::random_element();

    check_result &= this->eq_U_.verify(x_eq_U_, U_V_);
    check_result &= this->eq_U_PI_.verify(x_eq_U_PI_, U_PI_V_); 
    check_result &= this->sum_A_PRIME_.verify(x_sum_A_PRIME, row_vector<FieldT>::all_zero(this->k_));
    check_result &= this->sum_B_PRIME_.verify(x_sum_B_PRIME, row_vector<FieldT>::all_zero(this->k_));
    check_result &= this->same_prod_AP_BP_.verify(mid_challenges, challenges);

    if (check_result == true) {
        if (output) printf("pp_perm<FieldT>::verify() \033[32mpass\033[37m\n\n");
        return true;
    }
    if (output) printf("pp_perm<FieldT>::verify() \033[31mfail\033[37m\n\n");
    return false;
}

template <typename FieldT>
bool permutation_test() {
    
    size_t mu, n, col_num, row_num;
    mu = std::rand() % 4 + 3;
    n = std::rand() % 4 + 3;
    // mu = 2;
    // n = 2;
    col_num = std::pow(2UL, mu) * n;
    row_num = std::pow(2UL, mu) * n;

    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);

    size_t cycle_len = col_num / 2;
    size_t cycle_num = row_num / 2;
    std::vector<size_t> dim_limits = {row_num, col_num};

    permutation<tuple_dim2_t, 2> PI = permutation<tuple_dim2_t, 2>::random_permutation(cycle_len, cycle_num, 2, dim_limits);

    // 将 A 矩阵中在同一 cycle 位置的元素设置为相同值
    size_t row, col;
    for (size_t i = 0; i < cycle_num; i++) {
        FieldT val = FieldT::random_element();
        const std::vector<tuple_dim2_t> cyc_vec = PI.get_cycle(i).get_content_vector();
        for (size_t j = 0; j < cycle_len; j++) {
            std::tie(row, col) = cyc_vec[j];
            A.set_item(row, col, val);
        }
    }

    pp_perm<FieldT> perm(col_num, mu, n, A, A, PI);
    // same_prod.is_satisfy();
    if (A == apply_permutation(A, PI)) {
        printf("permutation_test.is_satisfy() \033[32m pass\033[37m\n");
    } else {
        printf("permutation_test.is_satisfy() \033[31m fail\033[37m\n");
    }

    /*
        prove part
    */
    std::map<std::string, FieldT> mid_challenges;
    mid_challenges["pp_prod_A_y"]   = FieldT::random_element();
    mid_challenges["pp_prod_A_x0"]  = FieldT::random_element();
    mid_challenges["pp_prod_B_y"]   = FieldT::random_element();
    mid_challenges["pp_prod_B_x0"]  = FieldT::random_element();
    mid_challenges["pp_shift_AB_y"] = FieldT::random_element();
    mid_challenges["pp_shift_AB_x0"] = FieldT::random_element();

    FieldT x = FieldT::random_element();
    FieldT y = FieldT::random_element();
    perm.prove(x, y, mid_challenges);

    /*
        verify part
    */
    std::map<std::string, FieldT> challenges;
    challenges["pp_prod_A_x"]   = FieldT::random_element();
    challenges["pp_prod_B_x"]   = FieldT::random_element();
    challenges["pp_shift_AB_x"] = FieldT::random_element();

    perm.verify(x, y, mid_challenges, challenges, true);
    return true;
    
}

#endif