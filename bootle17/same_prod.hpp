#ifndef BOOTLE_SAME_PROD_HPP
#define BOOTLE_SAME_PROD_HPP
#include "structs.hpp"
#include "prod.hpp"
#include "shift.hpp"

/*

initialize:
    pp_same_prod(const size_t mu, const size_t n, const size_t k,\
      const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B);


Now that we have an proof for the double-shift condition, and a Hadamardproduct proof, 
it is easy to construct an proof which shows that the product of all entries in a matrix A 
is the same of the product of all entries of a matrix B.

This is achieved by computing the partial products of entries of the matrix A, beginning 
with 1, and storing them in a matrix A1 with the same dimensions as A. The partial products, 
ending with the product of all elements of A, are stored in another matrix A2, and similarly 
for B, B1, B2. Now, A2 = A ◦ A1, and B2 = B ◦ B1 by design. Note that the product of all 
entries in A is the same as the product of the entries in B if and only if A1, A2, B1 and B2 
satisfy the double shift condition.
*/
template<typename FieldT>
class pp_same_prod {
public:
  pp_same_prod() {};
  pp_same_prod(const size_t mu, const size_t n, const size_t k, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B);

  bool is_satisfy() const;
  bool prove();
  bool verify(const bool output = false) const;

private:
  size_t mu_, m_, n_, col_num_;
  row_vector_matrix<FieldT> A_,  B_;
  row_vector_matrix<FieldT> A1_, B1_;
  row_vector_matrix<FieldT> A2_, B2_;

  pp_prod<FieldT> pp_prod_A_;
  pp_prod<FieldT> pp_prod_B_;
  pp_shift<FieldT> pp_shift_AB_;

};

template<typename FieldT>
bool same_prod_test();

#endif