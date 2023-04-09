#ifndef BOOTLE_SHIFT_TCC
#define BOOTLE_SHIFT_TCC
#include "shift.hpp"

template<typename FieldT>
pp_shift<FieldT>::pp_shift(size_t k, size_t mu, size_t n, const row_vector_matrix<FieldT>& A, \
    const row_vector_matrix<FieldT>& B, const row_vector_matrix<FieldT>& C, const row_vector_matrix<FieldT>& D) {
    this->k_ = k;
    this->mu_ = mu;
    this->m_ = std::pow(2ul, this->mu_);
    this->n_ = n;
    this->N_ = this->k_ * this->m_ * this->n_;
    printf("k = \033[32m%ld\033[37m, mu = \033[32m%ld\033[37m, m = \033[32m%ld\033[37m, n = \033[32m%ld\033[37m, N = \033[32m%ld\033[37m\n", this->k_, this->mu_, this->m_, this->n_, this->N_);
    this->A_ = row_vector_matrix<FieldT>(A);
    this->B_ = row_vector_matrix<FieldT>(B);
    this->C_ = row_vector_matrix<FieldT>(C);
    this->D_ = row_vector_matrix<FieldT>(D);
}

template<typename FieldT>
bool pp_shift<FieldT>::is_satisfy() {
    if (
    this->A_.get_row(0).get_item(0) != FieldT::one() || this->A_.get_row(0).get_item(0) != FieldT::one()\
    || this->B_.get_row(this->m_ * this->n_ - 1).get_item(this->k_ - 1) 
        != this->D_.get_row(this->m_ * this->n_ - 1).get_item(this->k_ - 1) ) {
        goto is_satisfy_fail;
    }
    for (size_t row = 0; row < this->m_ * this->n_; row++) {
        for (size_t col = 0; col < this->k_ - 1; col++) {
            if (this->A_.get_item(row, col+1) != this->B_.get_item(row, col) \
                || this->C_.get_item(row, col+1) != this->D_.get_item(row, col)) {
                goto is_satisfy_fail;
            }
        }
        if (row + 1 < this->m_ * this->n_
            && (this->A_.get_item(row+1, 0) != this->B_.get_item(row, this->k_-1) \
                || this->C_.get_item(row+1, 0) != this->D_.get_item(row, this->k_-1))) 
        {
            goto is_satisfy_fail;
        }
    }

is_satisfy_pass:
    printf("pp_shift<FieldT>::is_satisfy() \033[32mpass\033[37m\n");
    return true;

is_satisfy_fail:
    printf("pp_shift<FieldT>::is_satisfy() \033[31mfail\033[37m\n");
    return false;
}

template<typename FieldT>
void pp_shift<FieldT>::get_part_sum(const char choice, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_r) const {

    FieldT w_hat_prefix;
    size_t dim1 = this->n_; 
    size_t dim3 = this->m_;
    
    row_vector<FieldT> y_bold;
    // std::cout << "choice = " << choice << std::endl;
    switch ( choice ) {
        case 'A': 
            w_hat_prefix = FieldT::one();
            y_bold = row_vector<FieldT>(std::vector<FieldT>(this->y_exps_.begin(), this->y_exps_.begin()+this->k_));
            for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
                hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3);
                w_hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3); 
                for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
                    hat_r[idx_n][0][idx_m] = this->A_.get_row(idx_n * this->m_ + idx_m);
                    w_hat_r[idx_n][0][idx_m] = w_hat_prefix * y_bold * this->y_exps_[(idx_n * this->m_ + idx_m) * this->k_];
                }
            }
            break;
        case 'B': 
            w_hat_prefix = - this->y_exps_[1];
            y_bold = row_vector<FieldT>(std::vector<FieldT>(this->y_exps_.begin(), this->y_exps_.begin()+this->k_));
            for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
                hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3);
                w_hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3); 
                for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
                    hat_r[idx_n][0][idx_m] = this->B_.get_row(idx_n * this->m_ + idx_m);
                    w_hat_r[idx_n][0][idx_m] = w_hat_prefix * y_bold * this->y_exps_[(idx_n * this->m_ + idx_m) * this->k_];
                }
            }
            break;
        case 'C': 
            w_hat_prefix = - (this->y_exps_[this->N_] * this->y_exps_[this->N_]);
            y_bold = row_vector<FieldT>(std::vector<FieldT>(this->y_exps_inv_.begin(), this->y_exps_inv_.begin()+this->k_));
            for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
                hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3);
                w_hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3); 
                for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
                    hat_r[idx_n][0][idx_m] = this->C_.get_row(idx_n * this->m_ + idx_m);
                    w_hat_r[idx_n][0][idx_m] = w_hat_prefix * y_bold * this->y_exps_inv_[(idx_n * this->m_ + idx_m) * this->k_];
                }
            }
            break;
        case 'D': 
            w_hat_prefix = (this->y_exps_[this->N_] * this->y_exps_[this->N_] * this->y_exps_inv_[1]);
            y_bold = row_vector<FieldT>(std::vector<FieldT>(this->y_exps_inv_.begin(), this->y_exps_inv_.begin()+this->k_));
            for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
                hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3);
                w_hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3); 
                for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
                    hat_r[idx_n][0][idx_m] = this->D_.get_row(idx_n * this->m_ + idx_m);
                    w_hat_r[idx_n][0][idx_m] = w_hat_prefix * y_bold * this->y_exps_inv_[(idx_n * this->m_ + idx_m) * this->k_];
                }
            }
            break;
        default:
            break;
    }

    size_t dim2 = this->mu_ + 1;
    FieldT pre_challenge, pre_challenge_inverse;
    for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
        for (size_t idx_mu = 1; idx_mu < dim2; idx_mu++) {

            size_t i_step = std::pow(2, idx_mu);
            size_t i_upper_bound = dim3 / i_step;
            hat_r[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);
            w_hat_r[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);

            pre_challenge = this->xs_[idx_mu - 1];
            pre_challenge_inverse = pre_challenge.inverse();
            // std::cout << "pre_challenge " << pre_challenge << ", pre_challenge_inverse = " << pre_challenge_inverse << std::endl;
            for (size_t idx_m = 0; idx_m < i_upper_bound; idx_m += 1) {
                row_vector<FieldT>& a_0 = hat_r[idx_n][idx_mu-1][2 * idx_m];
                row_vector<FieldT>& a_1 = hat_r[idx_n][idx_mu-1][2 * idx_m + 1];
                row_vector<FieldT>& b_0 = w_hat_r[idx_n][idx_mu-1][2 * idx_m];
                row_vector<FieldT>& b_1 = w_hat_r[idx_n][idx_mu-1][2 * idx_m + 1];
                hat_r[idx_n][idx_mu][idx_m] = a_0 + a_1 * pre_challenge;
                w_hat_r[idx_n][idx_mu][idx_m] = b_0 + b_1 * pre_challenge_inverse;
            }
        }
    }
    return;
}

template<typename FieldT>
void pp_shift<FieldT>::get_a_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& a_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_a_r) {
    this->get_part_sum('A', a_hat_r, w_hat_a_r);
}

template<typename FieldT>
void pp_shift<FieldT>::get_b_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& b_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_b_r) {
    this->get_part_sum('B', b_hat_r, w_hat_b_r);
}

template<typename FieldT>
void pp_shift<FieldT>::get_c_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& c_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_c_r) {
    this->get_part_sum('C', c_hat_r, w_hat_c_r);
}

template<typename FieldT>
void pp_shift<FieldT>::get_d_hat_related(std::vector<std::vector<std::vector<row_vector<FieldT> > > >& d_hat_r, std::vector<std::vector<std::vector<row_vector<FieldT> > > >& w_hat_d_r) {
    this->get_part_sum('D', d_hat_r, w_hat_d_r);
}



template<typename FieldT>
bool pp_shift<FieldT>::check_related(const relate_t& a_rel,  const relate_t& b_rel, const relate_t& c_rel, const relate_t& d_rel, \
        const relate_t& w_a_rel,  const relate_t& w_b_rel, const relate_t& w_c_rel, const relate_t& w_d_rel) const {
    
    // std::cout << "enter check_related" << std::endl;
    FieldT x = FieldT::random_element();
    std::vector<FieldT> x_exps = get_exps(x, this->n_);
    std::vector<FieldT> x_exps_inv = get_exps(x.inverse(), this->n_);

    row_vector<FieldT> a_hat_1 = this->a0_;
    row_vector<FieldT> b_hat_1 = this->b0_;
    row_vector<FieldT> c_hat_1 = this->c0_;
    row_vector<FieldT> d_hat_1 = this->d0_;
    row_vector<FieldT> w_hat_a_1 = row_vector<FieldT>::all_zero(this->k_);
    row_vector<FieldT> w_hat_b_1 = row_vector<FieldT>::all_zero(this->k_);
    row_vector<FieldT> w_hat_c_1 = row_vector<FieldT>::all_zero(this->k_);
    row_vector<FieldT> w_hat_d_1 = row_vector<FieldT>::all_zero(this->k_);
    for (size_t j = 0; j < this->n_; j++) {
        a_hat_1   += a_rel[j][this->mu_][0]   * x_exps[j+1];
        b_hat_1   += b_rel[j][this->mu_][0]   * x_exps[j+1];
        c_hat_1   += c_rel[j][this->mu_][0]   * x_exps[j+1];
        d_hat_1   += d_rel[j][this->mu_][0]   * x_exps[j+1];
        w_hat_a_1 += w_a_rel[j][this->mu_][0] * x_exps_inv[j+1];
        w_hat_b_1 += w_b_rel[j][this->mu_][0] * x_exps_inv[j+1];
        w_hat_c_1 += w_c_rel[j][this->mu_][0] * x_exps_inv[j+1];
        w_hat_d_1 += w_d_rel[j][this->mu_][0] * x_exps_inv[j+1];
    }

    row_vector<FieldT> a_hat_2(this->k_, FieldT::zero());
    row_vector<FieldT> b_hat_2(this->k_, FieldT::zero());
    row_vector<FieldT> c_hat_2(this->k_, FieldT::zero());
    row_vector<FieldT> d_hat_2(this->k_, FieldT::zero());
    FieldT e_hat_2;

    this->open(a_hat_2, b_hat_2, c_hat_2, d_hat_2, e_hat_2, x);
    row_vector<FieldT> w_hat_a_2 = this->w_hat_a(this->y_exps_[1], x, this->xs_exps_);
    row_vector<FieldT> w_hat_b_2 = this->w_hat_b(this->y_exps_[1], x, this->xs_exps_);
    row_vector<FieldT> w_hat_c_2 = this->w_hat_c(this->y_exps_[1], x, this->xs_exps_);
    row_vector<FieldT> w_hat_d_2 = this->w_hat_d(this->y_exps_[1], x, this->xs_exps_);
    if (
        a_hat_1 == a_hat_2 && w_hat_a_1 == w_hat_a_2 &&
        b_hat_1 == b_hat_2 && w_hat_b_1 == w_hat_b_2 &&
        c_hat_1 == c_hat_2 && w_hat_c_1 == w_hat_c_2 &&
        d_hat_1 == d_hat_2 && w_hat_d_1 == w_hat_d_2
    ) {
        /*
            这里通过，说明 related 的求解和 除 er 外的 open 是没有问题的
            如果这里通过, 且验证不通过, 那么可能错误的地方就是 er 的 open, 或 f_plus_s_/f_sub_s_/gr_s_ 的计算出错
        */ 
        // printf("pp_shift<FieldT>::check_related \033[32mpass\033[37m\n");
        return true;
    }
    // printf("pp_shift<FieldT>::check_related \033[31mfail\033[37m\n");
    return false;
}

template<typename FieldT>
bool pp_shift<FieldT>::prove(const FieldT& y, const FieldT& x0) {
    this->a0_ = row_vector<FieldT>::random(this->k_);
    this->b0_ = row_vector<FieldT>::random(this->k_);
    this->c0_ = row_vector<FieldT>::random(this->k_);
    this->d0_ = row_vector<FieldT>::random(this->k_);

    this->y_exps_ = get_exps(y, this->N_);
    this->y_exps_inv_ = get_exps(y.inverse(), this->N_);
    FieldT dot_product_t = row_vector<FieldT>(this->y_exps_).dot_product(row_vector<FieldT>(this->y_exps_inv_));    // y 的计算没有错
    this->xs_ = fake_xs(x0, this->mu_);
    this->xs_exps_ = get_compress_xs_related(this->xs_);

    std::vector<std::vector<std::vector<row_vector<FieldT> > > > a_hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > b_hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > c_hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > d_hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > w_hat_a_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > w_hat_b_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > w_hat_c_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    std::vector<std::vector<std::vector<row_vector<FieldT> > > > w_hat_d_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
    get_a_hat_related(a_hat_r, w_hat_a_r);
    get_b_hat_related(b_hat_r, w_hat_b_r);
    get_c_hat_related(c_hat_r, w_hat_c_r);
    get_d_hat_related(d_hat_r, w_hat_d_r);

    // 计算 f_plus_s_ 和 f_sub_s_
    this->f_plus_s_ = row_vector<FieldT>(this->mu_, FieldT::zero());
    this->f_sub_s_ = row_vector<FieldT>(this->mu_, FieldT::zero());

    for (size_t t = 0; t < this->mu_; t++) {
        FieldT f_plus_t = FieldT::zero();
        FieldT f_sub_t  = FieldT::zero();

        size_t step = std::pow(2, t+1);
        size_t upper_bound = this->m_ / step;
        for (size_t j = 0; j < this->n_; j += 1) {
            for (size_t i = 0; i < upper_bound; i += 1) {
                f_plus_t += row_vector<FieldT>::dot_product(a_hat_r[j][t][2*i+1], w_hat_a_r[j][t][2*i]);                
                f_plus_t += row_vector<FieldT>::dot_product(b_hat_r[j][t][2*i+1], w_hat_b_r[j][t][2*i]);
                f_plus_t += row_vector<FieldT>::dot_product(c_hat_r[j][t][2*i+1], w_hat_c_r[j][t][2*i]);
                f_plus_t += row_vector<FieldT>::dot_product(d_hat_r[j][t][2*i+1], w_hat_d_r[j][t][2*i]);

                f_sub_t  += row_vector<FieldT>::dot_product(w_hat_a_r[j][t][2*i+1], a_hat_r[j][t][2*i]);
                f_sub_t  += row_vector<FieldT>::dot_product(w_hat_b_r[j][t][2*i+1], b_hat_r[j][t][2*i]);
                f_sub_t  += row_vector<FieldT>::dot_product(w_hat_c_r[j][t][2*i+1], c_hat_r[j][t][2*i]);
                f_sub_t  += row_vector<FieldT>::dot_product(w_hat_d_r[j][t][2*i+1], d_hat_r[j][t][2*i]);
            }
        }
        this->f_plus_s_.set_item(t, f_plus_t);
        this->f_sub_s_.set_item(t, f_sub_t);
    }

    // 计算 errors
    std::vector<FieldT> error_vec (2 * this->n_, FieldT::zero());
    const size_t dim2_idx = this->mu_;
    for (long j_a = 1; j_a <= this->n_; j_a++) {
        for (long j_b = 1; j_b <= this->n_; j_b++) {
            if (j_a == j_b) continue;
            error_vec[j_a - j_b + this->n_] += row_vector<FieldT>::dot_product(a_hat_r[j_a-1][dim2_idx][0], w_hat_a_r[j_b-1][dim2_idx][0]);
            error_vec[j_a - j_b + this->n_] += row_vector<FieldT>::dot_product(b_hat_r[j_a-1][dim2_idx][0], w_hat_b_r[j_b-1][dim2_idx][0]);
            error_vec[j_a - j_b + this->n_] += row_vector<FieldT>::dot_product(c_hat_r[j_a-1][dim2_idx][0], w_hat_c_r[j_b-1][dim2_idx][0]);
            error_vec[j_a - j_b + this->n_] += row_vector<FieldT>::dot_product(d_hat_r[j_a-1][dim2_idx][0], w_hat_d_r[j_b-1][dim2_idx][0]);
        }
    }
    for (long j = 1; j <= this->n_; j++) {
        error_vec[-j + this->n_] += row_vector<FieldT>::dot_product(this->a0_, w_hat_a_r[j-1][dim2_idx][0]);
        error_vec[-j + this->n_] += row_vector<FieldT>::dot_product(this->b0_, w_hat_b_r[j-1][dim2_idx][0]);
        error_vec[-j + this->n_] += row_vector<FieldT>::dot_product(this->c0_, w_hat_c_r[j-1][dim2_idx][0]);
        error_vec[-j + this->n_] += row_vector<FieldT>::dot_product(this->d0_, w_hat_d_r[j-1][dim2_idx][0]);
    }
    this->gr_s_ = row_vector<FieldT>(error_vec);

    assert(check_related(a_hat_r, b_hat_r, c_hat_r, d_hat_r, w_hat_a_r, w_hat_b_r, w_hat_c_r, w_hat_d_r) == true);

    return true;
}

template<typename FieldT>
row_vector<FieldT> pp_shift<FieldT>::w_hat_a(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const {
    std::vector<FieldT> y_exps = get_exps(y, this->N_);
    std::vector<FieldT> x_exps = get_exps(x, this->n_);

    FieldT sum = FieldT::zero();
    for (size_t j = 1; j <= this->n_; j++) {
        size_t base = (j-1) * this->m_;
        FieldT x_exp_j_inverse = x_exps[j].inverse();
        for (size_t i = 0; i <= this->m_-1; i++) {
            sum += y_exps[this->k_ * (i + base)] * this->xs_exps_[i].inverse() * x_exp_j_inverse;
        }
    }
    row_vector<FieldT> y_bold(std::vector<FieldT>(y_exps.begin(), y_exps.begin()+this->k_));  // 不包含 y_exps[k]
    return y_bold * sum;
}

template<typename FieldT>
row_vector<FieldT> pp_shift<FieldT>::w_hat_b(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const {
    return ( - y) * this->w_hat_a(y, x, compress_xs);
}

template<typename FieldT>
row_vector<FieldT> pp_shift<FieldT>::w_hat_c(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const {
    FieldT y_2N = this->y_exps_[this->N_] * this->y_exps_[this->N_];
    return (- y_2N) * this->w_hat_a(y.inverse(), x, compress_xs);
}

template<typename FieldT>
row_vector<FieldT> pp_shift<FieldT>::w_hat_d(const FieldT& y, const FieldT& x, const std::vector<FieldT>& compress_xs) const {
    FieldT y_2N_sub_1 = this->y_exps_[this->N_] * this->y_exps_[this->N_] * this->y_exps_inv_[1];
    return y_2N_sub_1 * this->w_hat_a(y.inverse(), x, compress_xs);
}

template<typename FieldT>
void pp_shift<FieldT>::open(row_vector<FieldT>& a_hat, row_vector<FieldT>& b_hat, row_vector<FieldT>& c_hat, \
row_vector<FieldT>& d_hat, FieldT& e_hat, const FieldT& x) const {
    std::vector<FieldT> x_exps = get_exps(x, this->n_);
    std::vector<FieldT> x_exps_inv = get_exps(x.inverse(), this->n_);
    a_hat.set_zero();
    b_hat.set_zero();
    c_hat.set_zero();
    d_hat.set_zero();
    a_hat += this->a0_;
    b_hat += this->b0_;
    c_hat += this->c0_;
    d_hat += this->d0_;
    for (size_t j = 1; j <= this->n_; j++) {
        FieldT x_exp_j = x_exps[j];
        for (size_t i = 0; i <= this->m_ - 1; i++) {
            a_hat += this->A_.get_row((j-1) * this->m_ + i) * this->xs_exps_[i] * x_exp_j;
            b_hat += this->B_.get_row((j-1) * this->m_ + i) * this->xs_exps_[i] * x_exp_j;
            c_hat += this->C_.get_row((j-1) * this->m_ + i) * this->xs_exps_[i] * x_exp_j;
            d_hat += this->D_.get_row((j-1) * this->m_ + i) * this->xs_exps_[i] * x_exp_j;
        }
    }

    e_hat *= FieldT::zero();
    e_hat += this->f_plus_s_.open(this->xs_);
    e_hat += this->f_sub_s_.open(row_vector<FieldT>(this->xs_).inverse().get_all_items());
    
    std::vector<FieldT> lc_gr(this->n_ * 2);
    lc_gr[this->n_] = FieldT::zero();
    for (size_t i = 0; i < this->n_; i++) {
        lc_gr[i] = x_exps_inv[this->n_-i]; // [1, n]
    }
    for (size_t i = this->n_+1; i < this->n_ * 2; i++) {   
        lc_gr[i] = x_exps[i - this->n_];   // [1, n-1]
    }
    e_hat += this->gr_s_.open(lc_gr);

    return;
}

template<typename FieldT>
bool pp_shift<FieldT>::verify(const FieldT& y, const FieldT& x0, const FieldT& x) const {
    row_vector<FieldT> a_hat(this->k_, FieldT::zero());
    row_vector<FieldT> b_hat(this->k_, FieldT::zero());
    row_vector<FieldT> c_hat(this->k_, FieldT::zero());
    row_vector<FieldT> d_hat(this->k_, FieldT::zero());
    FieldT e_hat;

    this->open(a_hat, b_hat, c_hat, d_hat, e_hat, x);

    row_vector<FieldT> w_hat_a_val = this->w_hat_a(y, x, this->xs_exps_);
    row_vector<FieldT> w_hat_b_val = this->w_hat_b(y, x, this->xs_exps_);
    row_vector<FieldT> w_hat_c_val = this->w_hat_c(y, x, this->xs_exps_);
    row_vector<FieldT> w_hat_d_val = this->w_hat_d(y, x, this->xs_exps_);
    FieldT left = a_hat.dot_product(w_hat_a_val) + b_hat.dot_product(w_hat_b_val) +\
                  c_hat.dot_product(w_hat_c_val) + d_hat.dot_product(w_hat_d_val);

    FieldT y_exp_2N = this->y_exps_[this->N_] * this->y_exps_[this->N_];
    FieldT right = FieldT::one() - y_exp_2N + e_hat;

    if (left != right) {
        printf("pp_shift<FieldT>::verify \033[31mfail\033[37m\n\n");
        return false;
    }
    printf("pp_shift<FieldT>::verify \033[32mpass\033[37m\n\n");
    return true;
}

template<typename FieldT>
void shift_test() {
    size_t mu, n, col_num, row_num;
    mu = std::rand() % 4 + 4;
    n = std::rand() % 4 + 4;
    col_num = std::pow(2UL, mu) * n;
    row_num = std::pow(2UL, mu) * n;

    row_vector_matrix<FieldT> B = row_vector_matrix<FieldT>::linear_grow(row_num, col_num, FieldT::one(), FieldT::one() * 2);
    row_vector_matrix<FieldT> D = row_vector_matrix<FieldT>::linear_grow(row_num, col_num, FieldT::one() * row_num * col_num, FieldT::one());

    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::linear_grow(row_num, col_num, FieldT::one() * (1 - 2), FieldT::one() * 2);
    A.set_item(0, 0, FieldT::one());

    row_vector_matrix<FieldT> C = row_vector_matrix<FieldT>::linear_grow(row_num, col_num, FieldT::one() * (row_num * col_num - 1), FieldT::one());
    C.set_item(0, 0, FieldT::one());

    pp_shift<FieldT> shift = pp_shift<FieldT>(col_num, mu, n, A, B, C, D);
    shift.is_satisfy();

    FieldT y  = FieldT::random_element();
    FieldT x0 = FieldT::random_element();
    FieldT x  = FieldT::random_element();

    shift.prove(y, x0);
    shift.verify(y, x0, x);
}

#endif