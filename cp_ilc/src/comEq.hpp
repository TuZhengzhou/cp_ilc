#ifndef BOOTLE_COMEQ_HPP
#define BOOTLE_COMEQ_HPP
#include <iostream>
#include <libff/algebra/curves/public_params.hpp>
#include <libsnark/knowledge_commitment/knowledge_commitment.hpp>
#include "structs.hpp"

/*
  为使 bootle17 支持承诺输入的证明而设计的 comEq 子协议

  初始化函数为: 

  pp_comEq(const size_t mu, const size_t n, const size_t col_num, const size_t com_num\
          const row_vector_matrix<FieldT> &A, const std::vector<FieldT> &rs,\
          const libff::G1<ppT> &g_base, const HType &h_base, const std::vector<CommitT> &Com);

  其中 com_num 表示实际公开承诺的数量。因为填充的关系，实际承诺数量可能比给定矩阵 A 中元素个数少
*/

template<typename FieldT, typename ppT, typename HType>
class pp_comEq {

typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
typedef libsnark::knowledge_commitment<libff::G1<ppT>, HType> CommitT;

public:
  pp_comEq() {};
  /* 向 ILC 通道承诺矩阵 A */
  pp_comEq(const size_t mu, const size_t n, const size_t col_num, const size_t com_num, \
          const row_vector_matrix<FieldT> &A, const std::vector<FieldT> &rs,\
          const libff::G1<ppT> &g_base, const HType &h_base, const std::vector<CommitT> &Com);

  bool is_satisfy() const;
  bool prove();
  bool verify(const bool output = false) const;

private:
  // 矩阵尺寸参数: (m_ * n_) * k_ = (2^mu_ * n_) * k_. k_ 列数, (m_ * n_) 行数
  size_t mu_, m_, n_, row_num_, col_num_;
  size_t N_;  // N_ = (m_ * n_) * k_

  /* 实际承诺个数 */
  size_t com_num_;

  /* 公开参数：证据矩阵 A_ 中元素的承诺向量 Coms_, 生成承诺使用的两个底数 */
  std::vector<CommitT> Coms_;
  libff::G1<ppT>       g_base_;
  HType       h_base_;

  /* 证据：输入矩阵 A_, 和生成承诺 Coms_ 时对应使用的随机数向量 rs_ */
  row_vector_matrix<FieldT> A_;
  std::vector<FieldT>       rs_;

  /* V -> ILC -> P  随机挑战 */
  FieldT y_;
  std::vector<FieldT> compress_xs_;
  /* 随机挑战派生值 */
  std::vector<FieldT> y_exps_;    // y^0, y^1, y^2, ...
  std::vector<FieldT> xs_exps_;   // x0^i0 * x1^i1 * ... * x_(mu-1)^i_(mu-1), for i in [0, 2^(mu)-1], i_t (t in [0, mu-1]) 是 i 的比特拆分

  /* 盲化因子 */
  row_vector<FieldT> a0_;    // blinding factor
  row_vector<FieldT> e0_;

  /* P -> V */
  libff::G1<ppT> E_; /* P -> V */
  HType H_; /* P -> V */

  /* P -> ILC */
  row_vector<FieldT> f_plus_s_, f_sub_s_; /* P -> ILC */
  row_vector<FieldT> gr_s_;               /* P -> ILC */

  bool check_related(const relate_t& a_rel, const relate_t& w_a_rel) const;

  void w_hat(const FieldT& x, row_vector<FieldT>& w_hat) const;
  
  void open(row_vector<FieldT>& a_hat, FieldT& e_hat, const FieldT& x) const;

  void get_part_sum(relate_t& hat_r, relate_t& w_hat_r) const;
};

template<typename FieldT, typename ppT, typename HType>
void comEq_test();

template<typename FieldT, typename ppT, typename HType>
void comEq_tests();

template<typename FieldT, typename ppT, typename HType>
void comEq_filling_tests();
#endif