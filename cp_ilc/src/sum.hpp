#ifndef BOOTLE_SUM_HPP
#define BOOTLE_SUM_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
/*
  bootle17 中的 sum 协议, prove 阶段证明者无工作, 故无 prove 函数

  初始化函数为: 
  
  pp_sum(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A, \
    const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);
*/
template<typename FieldT>
class pp_sum {
public:
  pp_sum() {};
  /* 向 ILC 通道承诺 A, B, C */
  pp_sum(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A, \
        const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);
  
  bool is_satisfy() const;
  bool verify(const row_vector<FieldT> row_vec, const bool output = false) const;

private:
  size_t row_num_;
  size_t col_num_;
  row_vector_matrix<FieldT> A;
  row_vector_matrix<FieldT> B;
  row_vector_matrix<FieldT> C;
};

template<typename FieldT>
void sum_test();
#endif