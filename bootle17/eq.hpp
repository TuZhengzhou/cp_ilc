#ifndef BOOTLE_EQ_HPP
#define BOOTLE_EQ_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"
/*
  bootle17 中的 eq 协议, prove 阶段证明者无工作, 故无 prove 函数

  初始化函数为: 

  pp_eq(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A)
*/
template<typename FieldT>
class pp_eq {
public:
  pp_eq() {};
  /* 向 ILC 通道承诺矩阵 A */
  pp_eq(const size_t row_num, const size_t col_num, const row_vector_matrix<FieldT>& A);

  bool is_satisfy(row_vector_matrix<FieldT>& other) const;
  bool verify(const row_vector_matrix<FieldT>& pub_matrix, const bool output = false) const;

private:
  size_t row_num_;
  size_t col_num_;
  row_vector_matrix<FieldT> A_;
};

template<typename FieldT>
void eq_test();
#endif