#ifndef BOOTLE_SUM_TCC
#define BOOTLE_SUM_TCC
#include <iostream>
#include <cassert>
#include "structs.hpp"
#include "sum.hpp"

template<typename FieldT>
pp_sum<FieldT>::pp_sum(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C) {
    assert(A.get_row_num() == row_num && B.get_row_num() == row_num && C.get_row_num() == row_num);
    assert(A.get_column_num() == col_num && B.get_column_num() == col_num && C.get_column_num() == col_num);

    // this->col_num = col_num;
    this->row_num_ = row_num;
    this->col_num_ = col_num;
    this->A = row_vector_matrix<FieldT>(A);
    this->B = row_vector_matrix<FieldT>(B);
    this->C = row_vector_matrix<FieldT>(C);
}

template<typename FieldT>
bool pp_sum<FieldT>::is_satisfy() const {
    assert(this->col_num_ > 0);

    if (this->C != this->A + this->B) {
        printf("pp_sum<FieldT>::is_satisfy() \033[31mfail\033[37m\n");
        return false;
    }
    printf("pp_sum<FieldT>::is_satisfy() \033[32mpass\033[37m\n");
    return true;
}

template<typename FieldT>
bool pp_sum<FieldT>::verify(const row_vector<FieldT> expected_result, const bool output) const {
    assert(this->col_num_ > 0);

    /* 随机挑战 */
    const FieldT challenge = FieldT::random_element();

    /* 生成查询所需的线性组合 */
    std::vector<FieldT> x_exps = get_exps(challenge, this->row_num_);
    std::vector<FieldT> x_exps_neg = ( - row_vector<FieldT>(x_exps)).get_all_items();
    std::vector<FieldT> linear_combination = std::vector<FieldT>(x_exps.begin()+1, x_exps.end());
    std::vector<FieldT> linear_combination_neg = std::vector<FieldT>(x_exps_neg.begin()+1, x_exps_neg.end());

    /* 向 ILC　通道查询 open 结果 */
    row_vector<FieldT> A_open = this->A.open(linear_combination);
    row_vector<FieldT> B_open = this->B.open(linear_combination);
    row_vector<FieldT> C_open = this->C.open(linear_combination_neg);
    row_vector<FieldT> open = A_open + B_open + C_open;

    /* 查询结果 与 期望结果对比 */
    if (open != expected_result) {
        if (output) printf("pp_sum<FieldT>::verify() \033[31mfail\033[37m\n\n");
        return false;
    }
    if (output) printf("pp_sum<FieldT>::verify() \033[32mpass\033[37m\n\n");
    return true;
}

template<typename FieldT>
void sum_test() {
    size_t col_num = std::rand() % 256;
    size_t row_num = std::rand() % 256;
    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
    row_vector_matrix<FieldT> B = row_vector_matrix<FieldT>::random(row_num, col_num);
    row_vector_matrix<FieldT> C = A + B;

    pp_sum<FieldT> sum = pp_sum<FieldT>(row_num, col_num, A, B, C);

    bool satisfy_result = sum.is_satisfy();
    assert(satisfy_result == true);
    
    /* 验证 */
    row_vector<FieldT> expected_result = row_vector<FieldT>::all_zero(col_num);
    bool verify_result = sum.verify(expected_result, true);
    assert(verify_result == true);

    return;
}

#endif