#ifndef BOOTLE_EQ_TCC
#define BOOTLE_EQ_TCC
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
#include "eq.hpp"

template<typename FieldT>
pp_eq<FieldT>::pp_eq(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A) {
    assert(A.get_row_num() == row_num);
    assert(A.get_column_num() == col_num);
    this->row_num_ = row_num;
    this->col_num_ = col_num;
    this->A_ = row_vector_matrix<FieldT>(A);
}

template<typename FieldT>
bool pp_eq<FieldT>::is_satisfy(row_vector_matrix<FieldT>& other) const {
    assert(this->col_num_ > 0UL);
    assert(this->col_num_ == other.get_column_num());
    assert(this->A_.get_row_num() == other.get_row_num());

    if (this->A_ != other) {
        printf("pp_eq<FieldT>::is_satisfy() \033[31mfail\033[37m\n");
        return false;
    }
    printf("pp_eq<FieldT>::is_satisfy() \033[32mpass\033[37m\n");
    return true;
}

template<typename FieldT>
bool pp_eq<FieldT>::verify(const row_vector_matrix<FieldT>& pub_matrix, const bool output) const {
    assert(this->col_num_ > 0UL);
    assert(this->col_num_ == pub_matrix.get_column_num());

    /* 随机挑战 */
    const FieldT challenge = FieldT::random_element();

    /* 向 ILC　通道查询 open 结果 */
    std::vector<FieldT> x_exps = get_exps<FieldT>(challenge, this->A_.get_row_num());
    row_vector<FieldT>  open = this->A_.open(std::vector<FieldT>(x_exps.begin()+1, x_exps.end()));

    /* 根据公共矩阵自行计算 期望计算结果 */
    row_vector<FieldT> expected_result = row_vector<FieldT>::all_zero(pub_matrix.get_column_num());
    for (size_t i = 0; i < pub_matrix.get_row_num(); i++) {
        expected_result += pub_matrix.get_row(i) * x_exps[i+1];
    }

    /* 比较 open 结果与 期望计算结果 */
    if (open != expected_result) {
        if (output) printf("pp_eq<FieldT>::verify() \033[31mfail\033[37m\n\n");
        return false;
    }
    if (output) printf("pp_eq<FieldT>::verify() \033[32mpass\033[37m\n\n");
    return true;
}

template<typename FieldT>
void eq_test() {
    size_t col_num = std::rand() % 256;
    size_t row_num = std::rand() % 256;
    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);

    pp_eq<FieldT> eq = pp_eq<FieldT>(row_num, col_num, A);

    bool satisfy_result = eq.is_satisfy(A);
    assert(satisfy_result == true);

    /* 验证 */
    bool verify_result = eq.verify(A, true);
    assert(verify_result == true);
}
#endif