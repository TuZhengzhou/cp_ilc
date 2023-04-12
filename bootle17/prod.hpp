#ifndef BOOTLE_PROD_HPP
#define BOOTLE_PROD_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"

/*
  prove hadmard product relation
  for matrix A, B, C, which's row is 2^mu * n, column is col_num prove A * B = C, where * means hadmard product

  初始化: 
  向 ILC 通道承诺 A, B, C
  pp_prod(const size_t mu, const size_t n, const size_t col_num, const row_vector_matrix<FieldT>& A, \
          const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);
*/
template<typename FieldT>
class pp_prod {

typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;

public:
  pp_prod() {}
  /* 向 ILC 通道承诺 A, B, C */
  pp_prod(const size_t mu, const size_t n, const size_t col_num, const row_vector_matrix<FieldT>& A, \
          const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);

  /* 检查明文是否符合 */
  bool is_satisfy() const;
  /*
    选取 a0, b0, c0
    计算 d_plus_s, d_sub_s, error_s
    并向 ILC 通道承诺 a0, b0, c0, d_plus_s, d_sub_s, error_s
  */
  bool prove();
  // 验证
  bool verify(const bool output = false) const;

private:
  size_t mu_, m_, n_, row_num_, col_num_;
  row_vector_matrix<FieldT> A_, B_, C_;
  
  FieldT y_;
  std::vector<FieldT> compress_xs_;

  row_vector<FieldT> a0, b0, c0;
  row_vector_matrix<FieldT> d_plus_s;
  row_vector_matrix<FieldT> d_sub_s;
  row_vector_matrix<FieldT> error_s;


  // 返回 y^0, y^1,..., y^(m-1), y^m, y^(m+1), ..., y^(m*n+m-1)
  std::vector<FieldT> get_y_exps(const FieldT& y) const;
  // 输入 x, 返回 x^0, x^1, x^2, ..., x^n
  std::vector<FieldT> get_x_exps(const FieldT& x) const;

  /*
    使用动态规划算法计算压缩得到的 每一个 a_j 中的 每一项 x_i 前的系数
    例如，将一个长为 8 的 向量压缩为 a_j
    层 0 (8 项): a_j_0  a_j_1  a_j_2  a_j_3  a_j_4  a_j_5  a_j_6  a_j_7

    层 1 (4 项): a_j_0 + a_j_1 * x_0, 
              a_j_2 + a_j_3 * x_0, 
              a_j_4 + a_j_5 * x_0, 
              a_j_6 + a_j_7 * x_0 

    层 2 (2 项): (a_j_0 + a_j_1 * x_0) + (a_j_2 + a_j_3 * x_0) * x_1, 
              (a_j_4 + a_j_5 * x_0) + (a_j_6 + a_j_7 * x_0) * x_1

    层 3 (1 项): ((a_j_0 + a_j_1 * x_0) + (a_j_2 + a_j_3 * x_0) * x_1) + ((a_j_4 + a_j_5 * x_0) + (a_j_6 + a_j_7 * x_0) * x_1) * x_2

    x_0 项前的系数   (对应层0, 4 项): a_j_1, a_j_3, a_j_5, a_j_7 (由 A 里 x_0^(1) * B 里 x_0^(0) 得到 x_0^(1) 次方)
    x_0^(-1) 前的系数(对应层0, 4 项): a_j_0, a_j_2, a_j_4, a_j_6 (由 A 里 x_0^(0) * B 里 x_0^(1) 得到 x_0^(-1) 次方)

    x_1 项前的系数   (对应层1, 2 项): a_j_2 + a_j_3 * x_0, a_j_6 + a_j_7 * x_0 (由 A 里 x_1^(1) * B 里 x_1^(0) 得到 x_1^(1) 次方)
    x_1^(-1) 前的系数(对应层1, 2 项): a_j_0 + a_j_1 * x_0, a_j_4 + a_j_5 * x_0 (由 A 里 x_1^(0) * B 里 x_1^(1) 得到 x_1^(-1) 次方)


    x_2 项前的系数   (对应层2, 1 项): (a_j_4 + a_j_5 * x_0) + (a_j_6 + a_j_7 * x_0) * x_1 (由 A 里 x_2^(1) * B 里 x_2^(0) 得到 x_2^(1) 次方)
    x_2^(-1) 前的系数(对应层2, 1 项): (a_j_0 + a_j_1 * x_0) + (a_j_2 + a_j_3 * x_0) * x_1 (由 A 里 x_2^(0) * B 里 x_2^(1) 得到 x_2^(-1) 次方)

    x 项前的系数     (对应层3, 1 项): ((a_j_0 + a_j_1 * x_0) + (a_j_2 + a_j_3 * x_0) * x_1) + ((a_j_4 + a_j_5 * x_0) + (a_j_6 + a_j_7 * x_0) * x_1) * x_2
  */ 
  void get_ab_related(relate_t& a_related, relate_t& b_related, const std::vector<FieldT>& y_exps) const;
  // 检查 a_related, b_realted 计算是否正确
  bool check_ab_related(const relate_t& a_related, const relate_t& b_related, const std::vector<FieldT>& y_exps) const;

  // 检查 error_s 是否计算正确
  bool check_error_s( const relate_t& a_related, \
                      const relate_t& b_related) const;

  row_vector<FieldT> open_A(const std::vector<FieldT>& y_exps, const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_exps) const ;
  row_vector<FieldT> open_B(const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_exps) const ;
  row_vector<FieldT> open_C(const std::vector<FieldT>& y_exps, const std::vector<FieldT>& compress_xs, const std::vector<FieldT>& x_exps) const ;
  row_vector<FieldT> open_d_plus_s(const std::vector<FieldT>& compress_xs) const ;
  row_vector<FieldT> open_d_sub_s (const std::vector<FieldT>& compress_xs) const ;
  row_vector<FieldT> open_error_s (const std::vector<FieldT>& x_exps) const ;
};

template<typename FieldT>
void prod_test();

#endif