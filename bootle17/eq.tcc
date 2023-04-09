#ifndef BOOTLE_EQ_TCC
#define BOOTLE_EQ_TCC
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
#include "eq.hpp"

template<typename FieldT>
pp_eq<FieldT>::pp_eq(const size_t num_column, const row_vector_matrix<FieldT>& A) {
    assert(A.get_column_num() == num_column);
    this->num_column = num_column;
    this->A = row_vector_matrix<FieldT>(A);
}

template<typename FieldT>
void pp_eq<FieldT>::set_num_column(const size_t num_column) {
    this->num_column = num_column;
}

template<typename FieldT>
void pp_eq<FieldT>::submit(const row_vector_matrix<FieldT>& matrix) {
    assert(this->num_column > 0);
    assert(this->num_column == matrix.get_column_num());

    this->A = row_vector_matrix<FieldT>(matrix);
    return;
}

template<typename FieldT>
bool pp_eq<FieldT>::is_satisfy(row_vector_matrix<FieldT>& other) const {
    assert(this->num_column > 0UL);
    assert(this->num_column == this->A.get_column_num() && this->num_column == other.get_column_num());
    assert(this->A.get_row_num() == other.get_row_num());

    if (this->A != other) {
        printf("pp_eq<FieldT>::is_satisfy() \033[31mfail\033[37m\n");
        return false;
    }
    printf("pp_eq<FieldT>::is_satisfy() \033[32mpass\033[37m\n");
    return true;
}

template<typename FieldT>
bool pp_eq<FieldT>::verify(const FieldT& challenge, const row_vector<FieldT> row_vec) const {
    assert(this->num_column > 0UL);
    assert(this->num_column == this->A.get_column_num() && this->num_column == row_vec.size());
    std::vector<FieldT> linear_combination = {challenge};
    size_t upper_bound = this->A.get_row_num();
    for (size_t i = 1; i < upper_bound; i++) {
        linear_combination.emplace_back(linear_combination[i-1] * challenge);
    }
    row_vector<FieldT> A_open = this->A.open(linear_combination);

    if (A_open != row_vec) {
        printf("pp_eq<FieldT>::verify() \033[31mfail\033[37m\n\n");
        return false;
    }
    printf("pp_eq<FieldT>::verify() \033[32mpass\033[37m\n\n");
    return true;
}

template<typename FieldT>
void eq_test() {
    size_t col_num = std::rand() % 256;
    size_t row_num = std::rand() % 256;
    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);

    pp_eq<FieldT> eq;
    eq.set_num_column(col_num);
    eq.submit(A);

    bool satisfy_result = eq.is_satisfy(A);
    assert(satisfy_result == true);

    FieldT challenage = FieldT::one() * (std::rand() % 10000000000000001UL);
    row_vector<FieldT> result = row_vector<FieldT>::all_zero(col_num);
    FieldT coeff = challenage;
    for (size_t i = 0; i < row_num; i++) {
        result += A.get_row(i) * coeff;
        coeff *= challenage;
    }
    bool verify_result = eq.verify(challenage, result);
    assert(verify_result == true);
}
#endif