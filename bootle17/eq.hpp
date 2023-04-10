#ifndef BOOTLE_EQ_HPP
#define BOOTLE_EQ_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
/*
    bootle17 中的 eq 协议, prove 阶段证明者无工作, 故无 prove 函数
*/
template<typename FieldT>
class pp_eq {
private:
    size_t num_column = 0;
    row_vector_matrix<FieldT> A;

public:
    pp_eq() {};
    pp_eq(const size_t num_column): num_column(num_column) {};
    pp_eq(const size_t num_column, const row_vector_matrix<FieldT>& A);

    void set_num_column(const size_t num_column);
    void submit(const row_vector_matrix<FieldT>& matrix);
    bool is_satisfy(row_vector_matrix<FieldT>& other) const;

    bool verify(const FieldT& challenge, const row_vector<FieldT> row_vec, const bool output = false) const;
};

template<typename FieldT>
void eq_test();
#endif