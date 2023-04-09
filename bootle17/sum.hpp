#ifndef BOOTLE_SUM_HPP
#define BOOTLE_SUM_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
/*
    bootle17 中的 sum 协议, prove 阶段证明者无工作, 故无 prove 函数
*/
template<typename FieldT>
class pp_sum {
private:
    size_t num_column = 0;
    row_vector_matrix<FieldT> A;
    row_vector_matrix<FieldT> B;
    row_vector_matrix<FieldT> C;

public:
    pp_sum() {};
    pp_sum(const size_t num_column): num_column(num_column) {};
    pp_sum(const size_t num_column, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);

    void set_num_column(const size_t num_column);
    void submit_part(const char Choice, const row_vector_matrix<FieldT>& matrix);
    void submit_all(const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);
    bool is_satisfy() const;

    bool verify(const FieldT& challenge, const row_vector<FieldT> row_vec) const;
};

template<typename FieldT>
void sum_test();
#endif