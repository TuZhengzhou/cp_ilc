#ifndef BOOTLE_SUM_TCC
#define BOOTLE_SUM_TCC
#include <iostream>
#include <cassert>
#include "structs.hpp"
#include "sum.hpp"

template<typename FieldT>
pp_sum<FieldT>::pp_sum(const size_t num_column, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C) {
    assert(A.get_column_num() == num_column);
    assert(B.get_column_num() == num_column);
    assert(C.get_column_num() == num_column);

    this->num_column = num_column;
    this->A = row_vector_matrix<FieldT>(A);
    this->B = row_vector_matrix<FieldT>(B);
    this->C = row_vector_matrix<FieldT>(C);
}

template<typename FieldT>
void pp_sum<FieldT>::set_num_column(const size_t num_column) {
    this->num_column = num_column;
}

template<typename FieldT>
void pp_sum<FieldT>::submit_part(const char Choice, const row_vector_matrix<FieldT>& matrix) {
    assert(Choice == 'A' || Choice == 'B' || Choice == 'C');

    assert(this->num_column > 0);
    assert(this->num_column == matrix.get_column_num());

    switch (Choice)
    {
    case 'A':
        this->A = row_vector_matrix<FieldT>(matrix);
        break;
    case 'B':
        this->B = row_vector_matrix<FieldT>(matrix);
        break;
    default:
        this->C = row_vector_matrix<FieldT>(matrix);
        break;
    }
    return;
}

template<typename FieldT>
void pp_sum<FieldT>::submit_all(const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C) {
    assert(this->num_column > 0);
    assert(this->num_column == A.get_column_num() && A.get_column_num() == B.get_column_num() && B.get_column_num() == C.get_column_num());

    this->A = row_vector_matrix<FieldT>(A);
    this->B = row_vector_matrix<FieldT>(B);
    this->C = row_vector_matrix<FieldT>(C);
    // printf("submit all \033[32msuccess\033[37m\n");
    return;
}

template<typename FieldT>
bool pp_sum<FieldT>::is_satisfy() const {
    
    assert(this->num_column > 0);
    assert(this->num_column == this->A.get_column_num() \
            && this->A.get_column_num() == this->B.get_column_num() \
            && this->B.get_column_num() == this->C.get_column_num());
    assert(this->A.get_row_num() == this->B.get_row_num() \
            && this->B.get_row_num() == this->C.get_row_num());

    if (this->C != this->A + this->B) {
        printf("pp_sum<FieldT>::is_satisfy() \033[31mfail\033[37m\n");
        return false;
    }
    printf("pp_sum<FieldT>::is_satisfy() \033[32mpass\033[37m\n");
    return true;
}

template<typename FieldT>
bool pp_sum<FieldT>::verify(const FieldT& challenge, const row_vector<FieldT> row_vec, const bool output) const {
    assert(this->num_column > 0);
    assert(this->num_column == this->A.get_column_num() \
            && this->num_column == this->B.get_column_num() \
            && this->num_column == this->C.get_column_num() \
            && this->num_column == row_vec.size());
    assert(this->A.get_row_num() == this->B.get_row_num() \
            && this->B.get_row_num() == this->C.get_row_num());

    std::vector<FieldT> linear_combination = {challenge};
    std::vector<FieldT> C_linear_combination = { - challenge};
    
    size_t upper_bound = this->A.get_row_num();
    for (size_t i = 1; i < upper_bound; i++) {
        linear_combination.emplace_back(linear_combination[i-1] * challenge);
        C_linear_combination.emplace_back( - linear_combination[i]);
    }

    row_vector<FieldT> A_open = this->A.open(linear_combination);
    row_vector<FieldT> B_open = this->B.open(linear_combination);
    row_vector<FieldT> C_open = this->C.open(C_linear_combination);

    // std::cout << A_open.get_all_items() << std::endl;
    // std::cout << B_open.get_all_items() << std::endl;
    // std::cout << C_open.get_all_items() << std::endl;
    // std::cout << row_vec.get_all_items() << std::endl;

    if (A_open + B_open + C_open != row_vec) {
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

    // 生成方法 1
    // pp_sum<FieldT> sum = pp_sum<FieldT>(col_num, A, B, C);

    // 生成方法 2
    // pp_sum<FieldT> sum = pp_sum<FieldT>(col_num);
    // sum.submit_all(A, B, C);

    // 生成方法 3
    pp_sum<FieldT> sum;
    sum.set_num_column(col_num);
    sum.submit_all(A, B, C);

    FieldT challenage = FieldT::one() * (std::rand() % 10000000000000001UL);
    row_vector<FieldT> result = row_vector<FieldT>::all_zero(col_num);

    bool satisfy_result = sum.is_satisfy();
    assert(satisfy_result == true);
    
    bool verify_result = sum.verify(challenage, result, true);
    assert(verify_result == true);

    return;
}

#endif