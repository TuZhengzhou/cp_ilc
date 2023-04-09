#ifndef BOOTLE_PROD_HPP
#define BOOTLE_PROD_HPP
#include <iostream>
#include "structs.hpp"
#include "structs.tcc"

template<typename FieldT>
class pp_prod {
private:
    size_t mu_, m_, n_, col_num_;
    row_vector_matrix<FieldT> A, B, C;
    
    row_vector<FieldT> a0, b0, c0;
    row_vector_matrix<FieldT> d_plus_s;
    row_vector_matrix<FieldT> d_sub_s;
    row_vector_matrix<FieldT> error_s;

    
    // 返回 y^0, y^1,..., y^(m-1), y^m, y^(m+1), ..., y^(m*n+m-1)
    std::vector<FieldT> get_y_related(const FieldT& y) const;
    // 输入 x, 返回 x^0, x^1, x^2, ..., x^n
    std::vector<FieldT> get_x_related(const FieldT& x) const;

    // 获取 a_related, b_related
    void get_ab_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& a_related, \
                        std::vector<std::vector<std::vector<row_vector<FieldT> > > >& b_related, \
                        const FieldT& x0, const std::vector<FieldT>& y_related) const;
    // 检查 a_related, b_realted 计算是否正确
    bool check_ab_related(const std::vector<std::vector<std::vector<row_vector<FieldT> > > >& a_related, \
                        const std::vector<std::vector<std::vector<row_vector<FieldT> > > >& b_related, \
                        const FieldT& x0, const std::vector<FieldT>& y_related) const;

    // 检查 error_s 是否计算正确
    bool check_error_s(const std::vector<std::vector<std::vector<row_vector<FieldT> > > >& a_related, \
                        const std::vector<std::vector<std::vector<row_vector<FieldT> > > >& b_related) const;

    row_vector<FieldT> open_A(const std::vector<FieldT>& y_related, const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_related) const ;
    row_vector<FieldT> open_B(const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_related) const ;
    row_vector<FieldT> open_C(const std::vector<FieldT>& y_related, const std::vector<FieldT>& compress_xs, const std::vector<FieldT>& x_related) const ;
    row_vector<FieldT> open_d_plus_s(const std::vector<FieldT>& compress_xs) const ;
    row_vector<FieldT> open_d_sub_s(const std::vector<FieldT>& compress_xs) const ;
    row_vector<FieldT> open_error_s(const std::vector<FieldT>& x_related) const ;

public:
    pp_prod() {}
    pp_prod(const size_t mu, const size_t n, const size_t col_num, const row_vector_matrix<FieldT>& A,  const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C);


    // 检查明文是否符合
    bool is_satisfy() const;
    // 证明, 赋值 a0, b0, c0, d_plus_s, d_sub_s, error_s
    bool prove(const FieldT& y, const FieldT& x0);
    // 验证
    bool verify(const FieldT& y, const FieldT& x0, const FieldT& x) const;
};

template<typename FieldT>
void prod_test();

#endif