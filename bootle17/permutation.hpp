#ifndef BOOTLE_PERMUTATION_HPP
#define BOOTLE_PERMUTATION_HPP
#include <vector>
#include "structs.hpp"
#include "eq.hpp"
#include "sum.hpp"
#include "same_prod.hpp"

typedef std::tuple<size_t, size_t> tuple_dim2_t;

/* pp_perm(
  const size_t k, const size_t mu, const size_t n, 
  const row_vector_matrix<FieldT>& A, 
  const row_vector_matrix<FieldT>& B, 
  const permutation<tuple_dim2_t, 2>& perm; 
*/
template <typename FieldT>
class pp_perm {
public:
  pp_perm() {};
  pp_perm(const size_t k, const size_t mu, const size_t n, const row_vector_matrix<FieldT>& A, \
          const row_vector_matrix<FieldT>& B, const permutation<tuple_dim2_t, 2>& perm);

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