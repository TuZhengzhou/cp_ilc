#ifndef BOOTLE_SHIFT_HPP
#define BOOTLE_SHIFT_HPP
#include "structs.hpp"


/*
prove of Double Shift condition

initialize:
    pp_shift(size_t k, size_t mu, size_t n, const row_vector_matrix<FieldT>& A, \
            const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C, \
            const row_vector_matrix<FieldT>& D);

Shift relation:
    Consider the matrices A and B, which have mn rows, given respectively by
    vectors ai,j , bi,j ∈ F^k, with 0 ≤ i ≤ m−1, 1 ≤ j ≤ n. The top-left element of A
    is a 1. Columns 2 up to k of A are equal to columns 1 up to k − 1 of B. Further,
    we can obtain the first column of A from the last column of B by deleting the
    last entry c and appending a 1. In this case, A is said to be the shift of B.

    B                                       A
    b1,1    b1,2    ...     b1,n            1       b1,1    b1,2    ....    b1,n-1
    b2,1    ...             b2,n            b1,n    b2,1    b2,2    ....    b2,n-1
    ...                                     b2,n    ....
    bm-1,1  ...             bm-1,n          ....
    bm,1    ...     bm,n-1  c               bm-1,n  bm,1    ....    ....    bm,n-1

Double Shift condition:
    for committed matrices A, B, C and D, we have A the shift
    of B, C the shift of D, and B and D have the same bottom-right-most entry
    bmn,k = dmn,k. This is referred to as the double-shift condition.
*/
template<typename FieldT>
class pp_shift {
    typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
public:
    
    pp_shift() {};
    pp_shift(size_t k, size_t mu, size_t n, const row_vector_matrix<FieldT>& A, \
            const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C, \
            const row_vector_matrix<FieldT>& D);

    bool is_satisfy();

    bool prove(const FieldT& y, const FieldT& x0);

    bool verify(const FieldT& y, const FieldT& x0, const FieldT& x, const bool output = false) const;

private:
    // 矩阵尺寸参数: (m_ * n_) * k_ = (2^mu_ * n_) * k_. k_ 列数, (m_ * n_) 行数
    size_t k_, mu_, m_, n_;
    size_t N_;  // N_ = (m_ * n_) * k_
    row_vector_matrix<FieldT> A_, B_, C_, D_;

    row_vector<FieldT> a0_, b0_, c0_, d0_;      // blinding factor
    row_vector<FieldT> f_plus_s_, f_sub_s_;
    row_vector<FieldT> gr_s_;

    std::vector<FieldT> y_exps_;        // y^0, y^1, y^2, ...
    std::vector<FieldT> y_exps_inv_;    // 
    std::vector<FieldT> xs_; 
    std::vector<FieldT> xs_exps_;   // x0^i0 * x1^i1 * ... * x_(mu-1)^i_(mu-1), for i in [0, 2^(mu)-1], i_t (t in [0, mu-1]) 是 i 的比特拆分

    bool check_related(const relate_t& a_rel,  const relate_t& b_rel, const relate_t& c_rel, const relate_t& d_rel,\
        const relate_t& w_a_rel,  const relate_t& w_b_rel, const relate_t& w_c_rel, const relate_t& w_d_rel) const;

    row_vector<FieldT> w_hat_a(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const;
    row_vector<FieldT> w_hat_b(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const;
    row_vector<FieldT> w_hat_c(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const;
    row_vector<FieldT> w_hat_d(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const;

    void open(row_vector<FieldT>& a_hat, row_vector<FieldT>& b_hat, row_vector<FieldT>& c_hat, \
                row_vector<FieldT>& d_hat, FieldT& e_hat, const FieldT& x) const;

    void get_part_sum(const char choice, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_r) const;
    void get_a_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& a_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_a_r);
    void get_b_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& b_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_b_r);
    void get_c_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& c_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_c_r);
    void get_d_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& d_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_d_r); 
};

template<typename FieldT>
void shift_test();

#endif