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

// template<typename FieldT>
// bool pp_eq<FieldT>::verify(const FieldT& challenge, const row_vector<FieldT>& row_vec, const bool output) const {
//     assert(this->num_column > 0UL);
//     assert(this->num_column == this->A.get_column_num() && this->num_column == row_vec.size());

//     std::vector<FieldT> x_exps = get_exps<FieldT>(challenge, this->A.get_row_num());
//     row_vector<FieldT> A_open = this->A.open(std::vector<FieldT>(x_exps.begin()+1, x_exps.end()));

//     if (A_open != row_vec) {
//         if (output) printf("pp_eq<FieldT>::verify() \033[31mfail\033[37m\n\n");
//         return false;
//     }
//     if (output) printf("pp_eq<FieldT>::verify() \033[32mpass\033[37m\n\n");
//     return true;
// }

template<typename FieldT>
bool pp_eq<FieldT>::verify(const FieldT& challenge, const row_vector_matrix<FieldT>& pub_matrix, const bool output) const {
    assert(this->num_column > 0UL);
    assert(this->num_column == this->A.get_column_num() && this->num_column == pub_matrix.get_column_num());

    std::vector<FieldT> x_exps = get_exps<FieldT>(challenge, this->A.get_row_num());
    row_vector<FieldT> A_open = this->A.open(std::vector<FieldT>(x_exps.begin()+1, x_exps.end()));

    row_vector<FieldT> result = row_vector<FieldT>::all_zero(pub_matrix.get_column_num());
    for (size_t i = 0; i < pub_matrix.get_row_num(); i++) {
        result += pub_matrix.get_row(i) * x_exps[i+1];
    }

    if (A_open != result) {
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

    pp_eq<FieldT> eq;
    eq.set_num_column(col_num);
    eq.submit(A);

    bool satisfy_result = eq.is_satisfy(A);
    assert(satisfy_result == true);

    FieldT challenge = FieldT::one() * (std::rand() % 10000000000000001UL);
    // row_vector<FieldT> result = row_vector<FieldT>::all_zero(col_num);
    // FieldT coeff = challenge;
    // for (size_t i = 0; i < row_num; i++) {
    //     result += A.get_row(i) * coeff;
    //     coeff *= challenge;
    // }
    // bool verify_result = eq.verify(challenge, result, true);
    bool verify_result = eq.verify(challenge, A, true);
    assert(verify_result == true);
}
#endif