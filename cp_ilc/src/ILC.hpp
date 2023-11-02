#ifndef BOOTLE_ILC_HPP
#define BOOTLE_ILC_HPP
#include "structs.hpp"
#include "eq.hpp"
#include "sum.hpp"
#include "prod.hpp"
#include "permutation.hpp"
#include "comEq.hpp"

/*
  pos_marker(const tuple_dim2_t& pos, const bool is_copy, const tuple_dim2_t& copy_from)
*/
struct pos_marker {
  tuple_dim2_t pos_;
  
  bool         is_copy_ = false;
  tuple_dim2_t copy_from_;
  std::vector<tuple_dim2_t> copies_;

  pos_marker(){};
  pos_marker(const tuple_dim2_t& pos) {
    this->pos_ = pos;
    this->is_copy_ = false;
  };
  pos_marker(const tuple_dim2_t& pos, const bool is_copy, const tuple_dim2_t& copy_from) {
    this->pos_ = pos;
    this->is_copy_ = is_copy;
    this->copy_from_ = copy_from;
  };

  void add_copy(const tuple_dim2_t& copy_pos) {
    this->copies_.push_back(copy_pos);
  }

  friend std::ostream& operator<<(std::ostream& out, pos_marker& marker) {
    out << "marker: ";
    out << "\tpos: " << std::get<0>(marker.pos_) << ", " << std::get<1>(marker.pos_) << std::endl;
    out << "is_copy: " << marker.is_copy_ << std::endl;
    if (marker.is_copy_) {
      out << "copy_from: " << std::get<0>(marker.copy_from_) << ", " << std::get<1>(marker.copy_from_) << std::endl;
    }
    for (size_t i = 0; i < marker.copies_.size(); i++) {
      out << "copies[" << i << "]: " << std::get<0>(marker.copies_[i]) << ", " << std::get<1>(marker.copies_[i]) << std::endl;
    }
    return out;
  }
};

/*
pp_ILC(const size_t mu_A, const size_t mu_M, const size_t mu_V, const size_t n_A, const size_t n_M, const size_t n_V, \
       const size_t col_num, const row_vector_matrix<FieldT>& V, const permutation<tuple_dim2_t, 2>& PI, const std::vector<size_t>& S);
*/
template <typename FieldT, typename ppT, typename HType>
class pp_ILC {

typedef libsnark::knowledge_commitment<libff::G1<ppT>, HType > CommitT;

public:

  pp_ILC() {};
  pp_ILC(const size_t mu_A, const size_t mu_M, const size_t mu_V, const size_t n_A, const size_t n_M, const size_t n_V, \
         const size_t col_num, const row_vector_matrix<FieldT>& V, const permutation<tuple_dim2_t, 2>& PI, \
         const std::vector<size_t>& S, const bool com_input);

  bool prove();
  bool verify(const bool output = false);

  /* 根据尺寸参数生成随机可通过的 ILC 矩阵和 排列 */
  static pp_ILC<FieldT, ppT, HType> random(const size_t mu_A, const size_t mu_M, const size_t mu_V, \
                                    const size_t n_A, const size_t n_M, const size_t n_V, \
                                    const size_t col_num, const bool com_input);
private:

  size_t mu_A_, mu_M_, mu_V_, n_A_, n_M_, n_V_, col_num_;
  size_t m_A_, m_M_;  
  size_t row_num_A_, row_num_M_, row_num_V_;  /* 加法行向量个数, 乘法行向量个数, 行向量总数 */
  
  permutation<tuple_dim2_t, 2> PI_;
  std::vector<size_t> S_;
  row_vector_matrix<FieldT> U_;
  
  row_vector_matrix<FieldT> V_;
  row_vector_matrix<FieldT> A_, B_, C_, D_, E_, F_;

  pp_eq<FieldT> eq_;
  pp_sum<FieldT> sum_;
  pp_prod<FieldT> prod_;
  pp_perm<FieldT> perm_;

  bool com_input_;
  pp_comEq<FieldT, ppT, HType> comEq_add_;
  pp_comEq<FieldT, ppT, HType> comEq_mult_;


  /* 
    给定已有的加法和乘法数量, 随机选择已有值作为源输入 
  */
  bool get_random_source(const size_t add_num, const size_t mutiply_num, FieldT& source_val, size_t& source_idx_in_V);
  
  /*
    给出目标位置和源位置
    如果源输入B 不是 copy 自其他源, 在源输入 B 中做标记
    如果源输入B copy 自其他源 A, 在源输入 A 中做标记
  */
  bool mark_copy(std::vector<std::vector<pos_marker> >& markers, const size_t target_idx, const size_t source_idx);
};



/* 将坐标转换为向量里的对应下标 */
size_t to_vector_index(const tuple_dim2_t& pos, const size_t col_num) {
  size_t row, col;
  std::tie(row, col) = pos;
  return row * col_num + col;
}

/* 将向量里的下表转换为坐标 */
tuple_dim2_t to_tuple_pos(const size_t vector_index, const size_t col_num) {
  size_t row, col;
  row = vector_index / col_num;
  col = vector_index % col_num;
  return std::make_tuple(row, col);
}

// template <typename FieldT, typename ppT, typename HType>
// bool ILC_test_origin(const size_t mu, double& prove_time, double& verifiy_time);

// template <typename FieldT, typename ppT, typename HType>
// bool ILC_test_with_com_input(const size_t mu, double& prove_time, double& verifiy_time);

/* 
  承诺数量不变，电路变大（行数固定，增加列数）
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_1();

/* 
  电路不变，承诺数量倍数增加
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_2();

/*
  电路和承诺数量比值不变（行数不变，承诺数量与列数相同）
*/
void ILC_test_compare_3();

#endif