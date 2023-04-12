#ifndef BOOTLE_PROD_TCC
#define BOOTLE_PROD_TCC
#include <iostream>
#include <cmath>
#include "structs.hpp"
#include "structs.tcc"
#include "prod.hpp"

template<typename FieldT>
pp_prod<FieldT>::pp_prod(const size_t mu, const size_t n, const size_t col_num, const row_vector_matrix<FieldT>& A,  \
                         const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C) :
    mu_(mu), m_(std::pow(2, mu)), n_(n), row_num_(std::pow(2, mu) * n), col_num_(col_num)
{
    assert(A.get_column_num() == col_num && B.get_column_num() == col_num && C.get_column_num() == col_num);
    assert(A.get_row_num() == this->row_num_ && B.get_row_num() == this->row_num_ && C.get_row_num() == this->row_num_);

    printf("mu = \033[32m%ld\033[37m, m = \033[32m%ld\033[37m, n = \033[32m%ld\033[37m, col_num = \033[32m%ld\033[37m\n", this->mu_, this->m_, this->n_, this->col_num_);
    this->A_ = row_vector_matrix<FieldT>(A);
    this->B_ = row_vector_matrix<FieldT>(B);
    this->C_ = row_vector_matrix<FieldT>(C);

    this->y_  = FieldT::random_element();
    this->compress_xs_ = row_vector<FieldT>::random(this->mu_).get_all_items();
}

template<typename FieldT>
std::vector<FieldT> pp_prod<FieldT>::get_y_exps(const FieldT& y) const {
    return get_exps(y, this->m_ * (this->n_ + 1));
}

template<typename FieldT>
std::vector<FieldT> pp_prod<FieldT>::get_x_exps(const FieldT& x) const {
    return get_exps(x, this->n_);
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_A(const std::vector<FieldT>& y_exps, const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_exps) const {
    std::vector<FieldT> linear_combination(this->row_num_);

    for (size_t j = 1; j <= this->n_; j++) {
        for (size_t i = 0; i < this->m_; i++) {
            linear_combination[(j-1) * this->m_ + i] = y_exps[j * this->m_ + i] * compress_xs_related[i] * x_exps[j];
        }
    }
    // std::cout << "openA linear_combination:\n" << linear_combination << std::endl;
    return this->a0 + this->A_.open(linear_combination);
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_B(const std::vector<FieldT>& compress_xs_related, const std::vector<FieldT>& x_exps) const {
    std::vector<FieldT> linear_combination(this->row_num_);

    for (size_t j = 1; j <= this->n_; j++) {
        for (size_t i = 0; i < this->m_; i++) {
            linear_combination[(j-1) * this->m_ + i] = (compress_xs_related[i] * x_exps[j]).inverse();
        }
    }
    // std::cout << "openB linear_combination:\n" << linear_combination << std::endl;
    return this->b0 + this->B_.open(linear_combination);
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_C(const std::vector<FieldT>& y_exps, const std::vector<FieldT>& compress_xs, const std::vector<FieldT>& x_exps) const {
    std::vector<FieldT> linear_combination(this->row_num_);

    for (size_t j = 1; j <= this->n_; j++) {
        for (size_t i = 0; i < this->m_; i++) {
            linear_combination[(j-1) * this->m_ + i] = y_exps[j * this->m_ + i];
        }
    }
    // std::cout << "openC linear_combination:\n" << linear_combination << std::endl;

    row_vector<FieldT> result = this->c0 + this->C_.open(linear_combination);
    result += this->open_d_plus_s(compress_xs);
    result += this->open_d_sub_s(compress_xs);
    result += this->open_error_s(x_exps);
    return result;
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_d_plus_s(const std::vector<FieldT>& compress_xs) const {
    std::vector<FieldT> compress_xs_ = compress_xs;
    // std::cout << "open_d_plus_s compress_xs:\n" << compress_xs_ << std::endl;
    return this->d_plus_s.open(compress_xs_);
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_d_sub_s(const std::vector<FieldT>& compress_xs) const {
    std::vector<FieldT> compress_xs_inverse;
    for (auto compress_x : compress_xs) {
        compress_xs_inverse.emplace_back(compress_x.inverse());
    }
    // std::cout << "open_d_sub_s compress_xs_inverse:\n" << compress_xs_inverse << std::endl;
    // return this->d_sub_s.open(compress_xs);
    return this->d_sub_s.open(compress_xs_inverse);
}

template<typename FieldT>
row_vector<FieldT> pp_prod<FieldT>::open_error_s(const std::vector<FieldT>& x_exps) const {
    std::vector<FieldT> linear_combination(this->n_ * 2 + 1);
    size_t mid = this->n_;
    linear_combination[mid] = FieldT::zero();
    for (size_t i = 0; i < this->n_; i++) {
        linear_combination[mid + i + 1] = x_exps[i+1];
        linear_combination[mid - i - 1] = x_exps[i+1].inverse();
    }
    // std::cout << "open_error_s linear_combination:\n" << linear_combination << std::endl;
    return this->error_s.open(linear_combination);
}

template<typename FieldT>
bool pp_prod<FieldT>::is_satisfy() const {
    if (this->A_ * this->B_ == this->C_) {
        return true;
    }
    return false;
}

template<typename FieldT>
void pp_prod<FieldT>::get_ab_related(relate_t& a_related, relate_t& b_related, const std::vector<FieldT>& y_exps) const {
    size_t dim1, dim2, dim3;
    dim1 = this->n_; 
    dim2 = this->mu_ + 1, 
    dim3 = this->m_;
    
    FieldT pre_challenge, pre_challenge_inverse;

    for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
        for (size_t idx_mu = 0; idx_mu < dim2; idx_mu++) {

            size_t i_step = std::pow(2, idx_mu);
            size_t i_upper_bound = dim3 / i_step;
            a_related[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);
            b_related[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);

            if (idx_mu == 0) {
                for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
                    a_related[idx_n][idx_mu][idx_m] = this->A_.get_row(idx_n * this->m_ + idx_m) * y_exps[(idx_n+1) * this->m_ + idx_m];
                    b_related[idx_n][idx_mu][idx_m] = this->B_.get_row(idx_n * this->m_ + idx_m);
                    // std::cout << "a_related[" << idx_n << "][" << idx_mu << "][" << idx_m << "] = " << a_related[idx_n][idx_mu][idx_m].get_item(0) << \
                    //     ", b_related[" << idx_n << "][" << idx_mu << "][" << idx_m << "] = " << b_related[idx_n][idx_mu][idx_m].get_item(0) << std::endl;
                }
            } else {
                pre_challenge = this->compress_xs_[idx_mu - 1];
                pre_challenge_inverse = pre_challenge.inverse();
                // std::cout << "pre_challenge " << pre_challenge << ", pre_challenge_inverse = " << pre_challenge_inverse << std::endl;
                for (size_t idx_m = 0; idx_m < i_upper_bound; idx_m += 1) {
                    row_vector<FieldT> a_0 = a_related[idx_n][idx_mu-1][2 * idx_m];
                    row_vector<FieldT> a_1 = a_related[idx_n][idx_mu-1][2 * idx_m + 1];
                    row_vector<FieldT> b_0 = b_related[idx_n][idx_mu-1][2 * idx_m];
                    row_vector<FieldT> b_1 = b_related[idx_n][idx_mu-1][2 * idx_m + 1];
                    a_related[idx_n][idx_mu][idx_m] = a_0 + a_1 * pre_challenge;
                    b_related[idx_n][idx_mu][idx_m] = b_0 + b_1 * pre_challenge_inverse;
                    // std::cout << "a_related[" << idx_n << "][" << idx_mu << "][" << idx_m << "] = " << a_related[idx_n][idx_mu][idx_m].get_item(0) << \
                    //     ", b_related[" << idx_n << "][" << idx_mu << "][" << idx_m << "] = " << b_related[idx_n][idx_mu][idx_m].get_item(0) << std::endl;
                }
            }
        }
    }
    return;
}

template<typename FieldT>
bool pp_prod<FieldT>::check_ab_related(const relate_t& a_related, const relate_t& b_related, \
                                       const std::vector<FieldT>& y_exps) const 
{
    FieldT x = FieldT::random_element();
    std::vector<FieldT> x_exps = this->get_x_exps(x);
    std::vector<FieldT> compress_xs_related = get_compress_xs_related(this->compress_xs_);
    row_vector<FieldT> a_dot1 = a0;
    row_vector<FieldT> b_dot1 = b0;
    for (size_t j = 0; j < this->n_; j++) {
        a_dot1 += a_related[j][this->mu_][0] * x_exps[j+1];
        b_dot1 += b_related[j][this->mu_][0] * x_exps[j+1].inverse();
    }
    row_vector<FieldT> a_dot2 = this->open_A(y_exps,  compress_xs_related, x_exps);
    row_vector<FieldT> b_dot2 = this->open_B(compress_xs_related, x_exps);
    if (a_dot1 == a_dot2 && b_dot1 == b_dot2) {
        // printf("ab_related \033[32mpass\033[37m\n");
        return true;
    }
    // printf("ab_related \033[31mfail\033[37m\n");
    return false;
}

template<typename FieldT>
bool pp_prod<FieldT>::check_error_s(const relate_t& a_related, const relate_t& b_related) const 
{
    size_t dim2_idx = this->mu_;
    FieldT x = FieldT::random_element();
    std::vector<FieldT> x_exps = this->get_x_exps(x);

    row_vector<FieldT> error_part_1 = row_vector<FieldT>::all_zero(this->col_num_);
    for (long j_a = 1; j_a <= this->n_; j_a++) {
        for (long j_b = 1; j_b <= this->n_; j_b++) {
            if (j_a == j_b) continue;
            error_part_1 += a_related[j_a-1][dim2_idx][0] * x_exps[j_a] * b_related[j_b-1][dim2_idx][0] * x_exps[j_b].inverse();
        }
    }
    for (long j = 1; j <= this->n_; j++) {
        error_part_1 += this->b0 * a_related[j-1][dim2_idx][0] * x_exps[j];
        error_part_1 += this->a0 * b_related[j-1][dim2_idx][0] * x_exps[j].inverse();
    }
    row_vector<FieldT> error_part_2 = this->open_error_s(x_exps);
    if (error_part_1 == error_part_2) {
        // printf("error_part \033[32mpass\033[37m\n");
        return true;
    }
    // printf("error_part \033[31mfail\033[37m\n");
    return false;
}

size_t findFirstDifferentBit(size_t x, size_t y) {
    size_t len = sizeof(x) * 8;
    for( size_t i = len-1; i >= 0; i-- ) {
        if ((x & (1ul << i)) != (y & (1ul << i))) {
            return i;
        }
    }
    return -1;
}


template<typename FieldT>
bool pp_prod<FieldT>::prove() {
    // 随机生成 a0, b0, c0
    this->a0 = row_vector<FieldT>::random(this->col_num_);
    this->b0 = row_vector<FieldT>::random(this->col_num_);
    this->c0 = a0 * b0;

    // 计算 y_exps
    std::vector<FieldT> y_exps = this->get_y_exps(this->y_);

    // 计算 a_related 和 b_related
    relate_t a_related(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    relate_t b_related(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    this->get_ab_related(a_related, b_related, y_exps);
    assert(this->check_ab_related(a_related, b_related, y_exps) == true);      // 检查 a_related 和 b_related 计算是否正确

    // 计算 d_plus_s 和 d_sub_s
    this->d_plus_s = row_vector_matrix<FieldT>(this->col_num_);
    this->d_sub_s = row_vector_matrix<FieldT>(this->col_num_);

    for (size_t t = 0; t < this->mu_; t++) {
        row_vector<FieldT> d_plus_t = row_vector<FieldT>::all_zero(this->col_num_);
        row_vector<FieldT> d_sub_t  = row_vector<FieldT>::all_zero(this->col_num_);

        size_t step = std::pow(2, t+1);
        size_t upper_bound = this->m_ / step;
        for (size_t j = 0; j < this->n_; j += 1) {
            for (size_t i = 0; i < upper_bound; i += 1) {
                d_plus_t += a_related[j][t][2*i+1] * b_related[j][t][2*i];
                d_sub_t  += b_related[j][t][2*i+1] * a_related[j][t][2*i];
            }
        }
        this->d_plus_s.add_row_vector(d_plus_t);
        this->d_sub_s.add_row_vector(d_sub_t);
        // std::cout << "d_plus_" << t << "= " << d_plus_t.get_item(0) << std::endl;
        // std::cout << "d_sub_"  << t << "= " << d_sub_t.get_item(0)  << std::endl;
    }

    // 计算 errors
    std::vector<row_vector<FieldT> > error_matrix (2 * this->n_ + 1, row_vector<FieldT>::all_zero(this->col_num_));
    size_t dim2_idx = this->mu_;
    for (long j_a = 1; j_a <= this->n_; j_a++) {
        for (long j_b = 1; j_b <= this->n_; j_b++) {
            if (j_a == j_b) continue;
            error_matrix[j_a - j_b + this->n_] += a_related[j_a-1][dim2_idx][0] * b_related[j_b-1][dim2_idx][0];
        }
    }
    for (long j = 1; j <= this->n_; j++) {
        error_matrix[j + this->n_] += this->b0 * a_related[j-1][dim2_idx][0];
        error_matrix[-j + this->n_] += this->a0 * b_related[j-1][dim2_idx][0];
    }
    this->error_s = row_vector_matrix<FieldT>(error_matrix, this->col_num_);
    // for (size_t i = 0; i < this->n_ * 2 + 1; i += 1) {
    //     std::cout << "this->error_s[" << i << "] = " << this->error_s.get_row(i).get_item(0) << std::endl;
    // }
    
    // 检查 errors 是否计算正确
    assert(check_error_s(a_related, b_related) == true);

    return true;
}

template<typename FieldT>
bool pp_prod<FieldT>::verify(const bool output) const {
    
    /* 选取随机挑战 */
    FieldT x = FieldT::random_element();

    std::vector<FieldT> y_exps              = this->get_y_exps(this->y_);
    std::vector<FieldT> compress_xs_related = get_compress_xs_related(this->compress_xs_);
    std::vector<FieldT> x_exps              = this->get_x_exps(x);

    row_vector<FieldT> a_dot = this->open_A(y_exps, compress_xs_related, x_exps);
    row_vector<FieldT> b_dot = this->open_B(compress_xs_related, x_exps);
    row_vector<FieldT> c_dot = this->open_C(y_exps, this->compress_xs_, x_exps);

    if (a_dot * b_dot != c_dot) {
        if (output) printf("pp_prod<FieldT>::verify \033[31mfail\033[37m\n\n");
        return false;
    }
    if (output) printf("pp_prod<FieldT>::verify \033[32mpass\033[37m\n\n");
    return true;
}

template<typename FieldT>
void prod_test() {
    size_t mu, n, col_num, row_num;
    mu = std::rand() % 4 + 4;
    n = std::rand() % 4 + 4;
    col_num = std::pow(2UL, mu) * n;
    row_num = std::pow(2UL, mu) * n;

    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
    row_vector_matrix<FieldT> B = row_vector_matrix<FieldT>::random(row_num, col_num);
    row_vector_matrix<FieldT> C = A * B;

    pp_prod<FieldT> prod(mu, n, col_num, A, B, C);
    if ( prod.is_satisfy() == true ) {
        printf("pp_prod<FieldT>::is_satisfy() \033[032mpass\033[037m\n");
    } else {
        printf("pp_prod<FieldT>::is_satisfy() \033[032mfail\033[037m\n");
    }

    prod.prove();

    FieldT x = FieldT::random_element();   
    prod.verify(true);
}

#endif