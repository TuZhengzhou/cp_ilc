#ifndef BOOTLE_PERMUTATION_HPP
#define BOOTLE_PERMUTATION_HPP
#include <vector>
#include "structs.tcc"
#include "eq.tcc"
#include "sum.tcc"
#include "same_prod.tcc"

typedef std::tuple<size_t, size_t> tuple_dim2_t;

/* pp_perm(
  const size_t mu, const size_t n, const size_t col_num, 
  const row_vector_matrix<FieldT>& A, 
  const row_vector_matrix<FieldT>& B, 
  const permutation<tuple_dim2_t, 2>& PI; 
*/
template <typename FieldT>
class pp_perm {
public:
  pp_perm() {};
  pp_perm(const size_t mu, const size_t n, const size_t col_num, const row_vector_matrix<FieldT>& A, \
          const row_vector_matrix<FieldT>& B, const permutation<tuple_dim2_t, 2>& PI);

  bool is_satisfy();
  bool prove();
  bool verify(const bool output = false) const;

private:
  size_t k_, mu_, m_, n_;
  size_t N_;  // N_ = (m_ * n_) * k_

  row_vector_matrix<FieldT> A_, B_;
  permutation<tuple_dim2_t, 2> PI_;

  FieldT x_, y_;

  row_vector_matrix<FieldT> V_, V_PI_;
  row_vector_matrix<FieldT> J_;
  row_vector_matrix<FieldT> U_, U_PI_;
  row_vector_matrix<FieldT> A_PRIME_, B_PRIME_;

  pp_eq<FieldT>  eq_U_;
  pp_eq<FieldT>  eq_U_PI_;
  pp_sum<FieldT> sum_A_PRIME_;
  pp_sum<FieldT> sum_B_PRIME_;
  pp_same_prod<FieldT> same_prod_AP_BP_;
};

template <typename FieldT, size_t N>
bool permutation_test();

#endif