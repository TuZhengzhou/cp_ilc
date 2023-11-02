#ifndef BOOTLE_COMEQ_TCC
#define BOOTLE_COMEQ_TCC
#include <iostream>
#include <type_traits>
#include "libff/common/profiling.hpp"
#include "structs.hpp"
#include "comEq.hpp"

template<typename FieldT, typename ppT, typename HType>
pp_comEq<FieldT, ppT, HType>::pp_comEq(const size_t mu, const size_t n, const size_t col_num, const size_t com_num, \
                                const row_vector_matrix<FieldT> &A, const std::vector<FieldT> &rs,\
                                const libff::G1<ppT> &g_base, const HType &h_base, const std::vector<CommitT> &Com) {

  // 矩阵尺寸参数: (m_ * n_) * k_ = (2^mu_ * n_) * k_. k_ 列数, (m_ * n_) 行数
  // size_t mu_, m_, n_, row_num_, col_num_;
  // size_t N_;  // N_ = (m_ * n_) * k_
  this->mu_ = mu;
  this->m_  = std::pow(2, this->mu_);
  this->n_  = n;
  this->row_num_ = this->m_ * this->n_;
  this->col_num_ = col_num;
  this->N_  = this->row_num_ * this->col_num_;
  // A.print();
  // libff::print_indent(); printf("mu = \033[32m%ld\033[37m, m = \033[32m%ld\033[37m, n = \033[32m%ld\033[37m, col_num = \033[32m%ld\033[37m, N = \033[32m%ld\033[37m\n", this->mu_, this->m_, this->n_, this->col_num_, this->N_);

  // /* 实际承诺个数 */
  // size_t com_num_;
  this->com_num_ = com_num;

  /* 检查尺寸是否跟描述一致 */
  assert(A.get_row_num() * A.get_column_num() == this->N_);
  assert(rs.size() == this->com_num_);
  assert(Com.size() == this->com_num_);

  // /* 公开参数：证据矩阵 A_ 中元素的承诺向量 Coms_, 生成承诺使用的两个底数 */
  // std::vector<CommitT> Coms_;
  // libff::G1<ppT>       g_base_;
  // HType       h_base_;
  this->Coms_ = std::move(Com);
  this->g_base_ = g_base;
  this->h_base_ = h_base;

  // /* 证据：输入矩阵 A_, 和生成承诺 Coms_ 时对应使用的随机数向量 rs_ */
  // row_vector_matrix<FieldT> A_;
  // std::vector<FieldT>       rs_;
  this->A_  = std::move(A);
  this->rs_ = std::move(rs);

  /* 随机挑战及其相关值 */
  // FieldT y_;
  // std::vector<FieldT> compress_xs_;
  // std::vector<FieldT> y_exps_;    // y^0, y^1, y^2, ...
  // std::vector<FieldT> xs_exps_;   // x0^i0 * x1^i1 * ... * x_(mu-1)^i_(mu-1), for i in [0, 2^(mu)-1], i_t (t in [0, mu-1]) 是 i 的比特拆分
  this->y_ = FieldT::random_element();
  this->compress_xs_ = row_vector<FieldT>::random(this->mu_).get_all_items();
  this->y_exps_      = get_exps(this->y_, this->N_ + 1);
  this->xs_exps_     = get_compress_xs_related(this->compress_xs_);
}

template<typename FieldT, typename ppT, typename HType>
void pp_comEq<FieldT, ppT, HType>::get_part_sum(relate_t& hat_r, relate_t& w_hat_r) const {
    
    size_t dim1 = this->n_;
    size_t dim3 = this->m_;
    row_vector<FieldT> y_bold = row_vector<FieldT>(std::vector<FieldT>(this->y_exps_.begin() + 1, this->y_exps_.begin()+this->col_num_ + 1));
    
    for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
        hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3);
        w_hat_r[idx_n][0] = std::vector<row_vector<FieldT> >(dim3); 
        for (size_t idx_m = 0; idx_m < dim3; idx_m += 1) {
            hat_r[idx_n][0][idx_m] = this->A_.get_row(idx_n * this->m_ + idx_m);
            w_hat_r[idx_n][0][idx_m] = y_bold * this->y_exps_[(idx_n * this->m_ + idx_m) * this->col_num_];
        }
    }
    size_t dim2 = this->mu_ + 1;
    FieldT pre_challenge, pre_challenge_inverse;
    for (size_t idx_n = 0; idx_n < dim1; idx_n++) {
        for (size_t idx_mu = 1; idx_mu < dim2; idx_mu++) {

            size_t i_step = std::pow(2, idx_mu);
            size_t i_upper_bound = dim3 / i_step;
            hat_r[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);
            w_hat_r[idx_n][idx_mu] = std::vector<row_vector<FieldT> >(i_upper_bound);

            pre_challenge = this->compress_xs_[idx_mu - 1];
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

template<typename FieldT, typename ppT, typename HType>
bool pp_comEq<FieldT, ppT, HType>::prove() {
  libff::enter_block("pp_comEq<FieldT, ppT>::prove()");

  /* 盲化因子 */
  // row_vector<FieldT> a0_;    // blinding factor
  // row_vector<FieldT> e0_;
  this->a0_ = row_vector<FieldT>::random(this->col_num_);
  this->e0_ = row_vector<FieldT>::random(1);
  

  // /* P -> V */
  // libff::G1<ppT> E_;
  this->E_ = this->e0_.get_item(0) * this->g_base_;

  // /* P -> V */
  // HType H_;
  FieldT exp_of_h = FieldT::zero();
  for (size_t i = 0; i < this->com_num_; i++) {
    exp_of_h += this->rs_[i] * this->y_exps_[i+1];
  }
  this->H_ = exp_of_h * this->h_base_;

  relate_t hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));
  relate_t w_hat_r(this->n_, std::vector<std::vector<row_vector<FieldT> > >(this->mu_ + 1));

  this->get_part_sum(hat_r, w_hat_r);

  /* P -> ILC */
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
              f_plus_t += row_vector<FieldT>::dot_product(hat_r[j][t][2*i+1], w_hat_r[j][t][2*i]);                
              f_sub_t  += row_vector<FieldT>::dot_product(w_hat_r[j][t][2*i+1], hat_r[j][t][2*i]);
          }
      }
      this->f_plus_s_.set_item(t, f_plus_t);
      this->f_sub_s_.set_item(t, f_sub_t);
  }

  /* P -> ILC */
  // 计算 errors
  std::vector<FieldT> error_vec (2 * this->n_, FieldT::zero());
  const size_t dim2_idx = this->mu_;
  for (long j_a = 1; j_a <= this->n_; j_a++) {
      for (long j_b = 1; j_b <= this->n_; j_b++) {
          if (j_a == j_b) continue;
          error_vec[j_a - j_b + this->n_] += row_vector<FieldT>::dot_product(hat_r[j_a-1][dim2_idx][0], w_hat_r[j_b-1][dim2_idx][0]);
      }
  }
  for (long j = 1; j <= this->n_; j++) {
      error_vec[-j + this->n_] += row_vector<FieldT>::dot_product(this->a0_, w_hat_r[j-1][dim2_idx][0]);
  }
  this->gr_s_ = row_vector<FieldT>(error_vec);

  // assert(check_related(a_hat_r, b_hat_r, c_hat_r, d_hat_r, w_hat_a_r, w_hat_b_r, w_hat_c_r, w_hat_d_r) == true);
  
  libff::leave_block("pp_comEq<FieldT, ppT>::prove()");

  return true;
}

template<typename FieldT, typename ppT, typename HType>
void pp_comEq<FieldT, ppT, HType>::w_hat(const FieldT& x, row_vector<FieldT>& w_hat) const {
  std::vector<FieldT> y_exps = get_exps(this->y_, this->N_ + 1);
  std::vector<FieldT> x_exps = get_exps(x, this->n_);

  FieldT sum = FieldT::zero();
  for (size_t j = 1; j <= this->n_; j++) {
      size_t base = (j-1) * this->m_;
      FieldT x_exp_j_inverse = x_exps[j].inverse();
      for (size_t i = 0; i <= this->m_-1; i++) {
          sum += y_exps[this->col_num_ * (i + base)] * this->xs_exps_[i].inverse() * x_exp_j_inverse;
      }
  }
  row_vector<FieldT> y_bold(std::vector<FieldT>(y_exps.begin() + 1, y_exps.begin()+this->col_num_ + 1));  // 不包含 y_exps[k]
  w_hat = y_bold * sum;
}


template<typename FieldT, typename ppT, typename HType>
void pp_comEq<FieldT, ppT, HType>::open(row_vector<FieldT>& a_hat, FieldT& e_hat, const FieldT& x) const {

  std::vector<FieldT> x_exps = get_exps(x, this->n_);
  std::vector<FieldT> x_exps_inv = get_exps(x.inverse(), this->n_);
  
  a_hat = row_vector<FieldT>(this->a0_);

  for (size_t j = 1; j <= this->n_; j++) {
      FieldT x_exp_j = x_exps[j];
      for (size_t i = 0; i <= this->m_ - 1; i++) {
          a_hat += this->A_.get_row((j-1) * this->m_ + i) * this->xs_exps_[i] * x_exp_j;
      }
  }

  /*
  e ̂=∑_(t=0)^(μ-1)▒(f_t^+ x_t+f_t^- x_t^(-1) ) +∑_(r=-n,r≠0)^(n-1)▒〖g_r X^r 〗+e_0 y^(N+1)
  */
  e_hat *= FieldT::zero();
  e_hat += this->f_plus_s_.open(this->compress_xs_);
  e_hat += this->f_sub_s_.open(row_vector<FieldT>(this->compress_xs_).inverse().get_all_items());
  
  std::vector<FieldT> lc_gr(this->n_ * 2);
  lc_gr[this->n_] = FieldT::zero();
  for (size_t i = 0; i < this->n_; i++) {
      lc_gr[i] = x_exps_inv[this->n_-i]; // [1, n]
  }
  for (size_t i = this->n_+1; i < this->n_ * 2; i++) {   
      lc_gr[i] = x_exps[i - this->n_];   // [1, n-1]
  }
  e_hat += this->gr_s_.open(lc_gr);

  e_hat += this->e0_.open({this->y_exps_[this->N_ + 1]});

  return; 
}

template<typename FieldT, typename ppT, typename HType>
bool pp_comEq<FieldT, ppT, HType>::verify(const bool output) const {
  libff::enter_block("pp_comEq<FieldT, ppT>::verify()");

  /* 选取随机挑战 */
  FieldT x = FieldT::random_element();

  FieldT e_hat;
  row_vector<FieldT> a_hat(this->col_num_, FieldT::zero());
  this->open(a_hat, e_hat, x);

  row_vector<FieldT> w_hat_val;
  this->w_hat(x, w_hat_val);

  libff::enter_block("calculate left");
  FieldT exp_of_g = a_hat.dot_product(w_hat_val) - e_hat;
  libff::G1<ppT> g_of_left = exp_of_g * this->g_base_ + this->y_exps_[this->N_+1] * this->E_;
  HType h_of_left = this->H_;
  CommitT left = CommitT(g_of_left, h_of_left);
  libff::leave_block("calculate left");

  const size_t comNum = this->Coms_.size();
  libff::enter_block("calculate right");
  CommitT right = CommitT::zero();
  for (size_t i = 1; i <= comNum; i++) {
    right = right + this->y_exps_[i] * this->Coms_[i-1];
  }
  libff::leave_block("calculate right");

  libff::leave_block("pp_comEq<FieldT, ppT>::verify()");

  if (left != right) {
      if (output) {
          libff::print_indent(); printf("pp_comEq<FieldT, ppT>::verify() \033[31mfail\033[37m\n");
      }
      return false;
  }
  if (output) {
      libff::print_indent(); printf("pp_comEq<FieldT, ppT>::verify() \033[32mpass\033[37m\n");
  }
  return true;
}

template<typename FieldT, typename ppT, typename HType>
bool pp_comEq<FieldT, ppT, HType>::is_satisfy() const {
  libff::enter_block("pp_comEq<FieldT, ppT>::is_satisfy()");

  std::vector<FieldT> items = this->A_.flatten();
  for (size_t i = 0; i < this->N_; i++) {
    CommitT com = CommitT(items[i] * this->g_base_, this->rs_[i] * this->h_base_);
    if (com != this->Coms_[i]) {
      libff::leave_block("pp_comEq<FieldT, ppT>::is_satisfy()");
      libff::print_indent(); printf("pp_comEq<FieldT, ppT>::is_satisfy() \033[31mfail\033[37m\n");
      return false;
    }
  }

  // FieldT y;
  // FieldT y_exps = get_exps(y, this->N_+1);
  FieldT exp_of_g = FieldT::zero();
  FieldT exp_of_h = FieldT::zero();
  CommitT right = CommitT::zero();
  for (size_t i = 0; i < this->N_; i++) {
    // CommitT com = CommitT(items[i] * this->g_base_, this->rs_[i] * this->h_base_);
    exp_of_g += items[i] * this->y_exps_[i+1];
    exp_of_h += this->rs_[i] * this->y_exps_[i+1];
    right = right + this->y_exps_[i+1] * this->Coms_[i];
  }
  CommitT left = CommitT(exp_of_g * this->g_base_, exp_of_h * this->h_base_);

  if (left != right) {
    libff::leave_block("pp_comEq<FieldT, ppT>::is_satisfy()");
    libff::print_indent(); printf("pp_comEq<FieldT, ppT>::is_satisfy() \033[31mfail\033[37m\n");
    return false;
  }
  libff::leave_block("pp_comEq<FieldT, ppT>::is_satisfy()");
  libff::print_indent(); printf("pp_comEq<FieldT, ppT>::is_satisfy() \033[32mpass\033[37m\n");
  return true;
}

template<typename FieldT, typename ppT>
void comEq_test() {
  typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
  typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G1<ppT>> CommitT;

  size_t mu, n, col_num, row_num, N;
  mu = std::rand() % 2 + 2;
  n = std::rand() % 2 + 2;
  col_num = std::pow(2UL, mu) * n;
  row_num = std::pow(2UL, mu) * n;
  N = row_num * col_num;

  row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
  row_vector_matrix<FieldT> R = row_vector_matrix<FieldT>::random(row_num, col_num);

  std::vector<FieldT> as = A.flatten();
  std::vector<FieldT> rs = R.flatten();
  std::vector<CommitT> coms = std::vector<CommitT>(N);

  libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  libff::G1<ppT> h_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  
  for (size_t i = 0; i < N; i++) {
    coms[i] = CommitT(as[i] * g_base, rs[i] * h_base);
  }
  
  /* will pass */
  libff::enter_block("comEq test1: designed to \033[32mpass\033[37m");
  pp_comEq<FieldT, ppT, libff::G1<ppT> > comEq1 = pp_comEq<FieldT, ppT, libff::G1<ppT> >(mu, n, col_num, N, A, rs, g_base, h_base, coms);
  comEq1.is_satisfy();
  comEq1.prove();
  bool verify_1 = comEq1.verify(true);
  assert(verify_1 == true);
  libff::leave_block("comEq test1: designed to \033[32mpass\033[37m");

  /* will fail */
  libff::enter_block("comEq test2: designed to \033[31mfail\033[37m");
  A.set_item(0, A.get_item(0)-FieldT::one());
  pp_comEq<FieldT, ppT, libff::G1<ppT> > comEq2 = pp_comEq<FieldT, ppT, libff::G1<ppT> >(mu, n, col_num, N, A, rs, g_base, h_base, coms);
  comEq2.is_satisfy();
  comEq2.prove();
  bool verify_2 = comEq2.verify(true);
  assert(verify_2 == false);
  libff::leave_block("comEq test2: designed to \033[31mfail\033[37m");
}

template<typename FieldT, typename ppT>
void test_single_h_in_G1(const size_t mu, const size_t n, const size_t col_num, double& prove_time, double& verifiy_time) {
  typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
  typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G1<ppT>> CommitT;

  size_t row_num, N;
  row_num = std::pow(2UL, mu) * n;
  N = row_num * col_num;

  row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
  row_vector_matrix<FieldT> R = row_vector_matrix<FieldT>::random(row_num, col_num);

  std::vector<FieldT> as = A.flatten();
  std::vector<FieldT> rs = R.flatten();
  std::vector<CommitT> coms = std::vector<CommitT>(N);

  libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  libff::G1<ppT> h_base = FieldT::random_element() * libff::G1<ppT>::G1_one;

  for (size_t i = 0; i < N; i++) {
    coms[i] = CommitT(as[i] * g_base, rs[i] * h_base);
  }
  
  std::chrono::_V2::steady_clock::time_point start, end;
  std::chrono::duration<int64_t, std::nano>  diff;
  pp_comEq<FieldT, ppT, libff::G1<ppT> > comEq1 = pp_comEq<FieldT, ppT, libff::G1<ppT> >(mu, n, col_num, N, A, rs, g_base, h_base, coms);

  start = std::chrono::steady_clock::now(); // 记录开始时间
  comEq1.prove();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  prove_time = std::chrono::duration<double, std::milli>(diff).count();

  start = std::chrono::steady_clock::now(); // 记录开始时间
  bool verify_1 = comEq1.verify();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

  assert(verify_1 == true);
}

template<typename FieldT, typename ppT>
void test_single_h_in_G2(const size_t mu, const size_t n, const size_t col_num, double& prove_time, double& verifiy_time) {
  typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
  typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT>> CommitT;

  size_t row_num, N;
  row_num = std::pow(2UL, mu) * n;
  N = row_num * col_num;

  row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
  row_vector_matrix<FieldT> R = row_vector_matrix<FieldT>::random(row_num, col_num);

  std::vector<FieldT> as = A.flatten();
  std::vector<FieldT> rs = R.flatten();
  std::vector<CommitT> coms = std::vector<CommitT>(N);

  libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  libff::G2<ppT> h_base = FieldT::random_element() * libff::G2<ppT>::G2_one;

  for (size_t i = 0; i < N; i++) {
    coms[i] = CommitT(as[i] * g_base, rs[i] * h_base);
  }
  
  std::chrono::_V2::steady_clock::time_point start, end;
  std::chrono::duration<int64_t, std::nano>  diff;
  pp_comEq<FieldT, ppT, libff::G2<ppT> > comEq1 = pp_comEq<FieldT, ppT, libff::G2<ppT> >(mu, n, col_num, N, A, rs, g_base, h_base, coms);

  start = std::chrono::steady_clock::now(); // 记录开始时间
  comEq1.prove();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  prove_time = std::chrono::duration<double, std::milli>(diff).count();

  start = std::chrono::steady_clock::now(); // 记录开始时间
  bool verify_1 = comEq1.verify();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

  assert(verify_1 == true);
}

// 承诺输入数目固定为col_num，输入矩阵的其他位置使用0填充
template<typename FieldT, typename ppT>
void test_single_filling_h_in_G1(const size_t mu, const size_t n, const size_t col_num, double& prove_time, double& verifiy_time) {
  typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
  typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G1<ppT> > CommitT;

  size_t row_num, N;
  row_num = std::pow(2UL, mu) * n;
  N = row_num * col_num;

  size_t com_num = col_num;

  row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>(col_num);
  A.add_row_vector(row_vector<FieldT>::random(col_num));
  for ( size_t i = 1; i < row_num; i++ ) {
    A.add_row_vector(row_vector<FieldT>::all_zero(col_num));
  }
  row_vector_matrix<FieldT> R = row_vector_matrix<FieldT>::random(1, col_num);

  std::vector<FieldT> as = A.flatten();
  std::vector<FieldT> rs = R.flatten();
  std::vector<CommitT> coms = std::vector<CommitT>(com_num);

  libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  libff::G1<ppT> h_base = FieldT::random_element() * libff::G1<ppT>::G1_one;

  for (size_t i = 0; i < com_num; i++) {
    coms[i] = CommitT(as[i] * g_base, rs[i] * h_base);
  }
  
  std::chrono::_V2::steady_clock::time_point start, end;
  std::chrono::duration<int64_t, std::nano>  diff;
  pp_comEq<FieldT, ppT, libff::G1<ppT> > comEq1 = pp_comEq<FieldT, ppT, libff::G1<ppT> >(mu, n, col_num, com_num, A, rs, g_base, h_base, coms);

  start = std::chrono::steady_clock::now(); // 记录开始时间
  comEq1.prove();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  prove_time = std::chrono::duration<double, std::milli>(diff).count();

  start = std::chrono::steady_clock::now(); // 记录开始时间
  bool verify_1 = comEq1.verify();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

  assert(verify_1 == true);
}

// 承诺输入数目固定为col_num，输入矩阵的其他位置使用0填充
template<typename FieldT, typename ppT>
void test_single_filling_h_in_G2(const size_t mu, const size_t n, const size_t col_num, double& prove_time, double& verifiy_time) {
  typedef std::vector<std::vector<std::vector<row_vector<FieldT> > > > relate_t;
  typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT> > CommitT;

  size_t row_num, N;
  row_num = std::pow(2UL, mu) * n;
  N = row_num * col_num;

  size_t com_num = col_num;

  row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>(col_num);
  A.add_row_vector(row_vector<FieldT>::random(col_num));
  for ( size_t i = 1; i < row_num; i++ ) {
    A.add_row_vector(row_vector<FieldT>::all_zero(col_num));
  }
  row_vector_matrix<FieldT> R = row_vector_matrix<FieldT>::random(1, col_num);

  std::vector<FieldT> as = A.flatten();
  std::vector<FieldT> rs = R.flatten();
  std::vector<CommitT> coms = std::vector<CommitT>(com_num);

  libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
  libff::G2<ppT> h_base = FieldT::random_element() * libff::G2<ppT>::G2_one;

  for (size_t i = 0; i < com_num; i++) {
    coms[i] = CommitT(as[i] * g_base, rs[i] * h_base);
  }
  
  std::chrono::_V2::steady_clock::time_point start, end;
  std::chrono::duration<int64_t, std::nano>  diff;
  pp_comEq<FieldT, ppT, libff::G2<ppT> > comEq1 = pp_comEq<FieldT, ppT, libff::G2<ppT> >(mu, n, col_num, com_num, A, rs, g_base, h_base, coms);

  start = std::chrono::steady_clock::now(); // 记录开始时间
  comEq1.prove();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  prove_time = std::chrono::duration<double, std::milli>(diff).count();

  start = std::chrono::steady_clock::now(); // 记录开始时间
  bool verify_1 = comEq1.verify();
  end = std::chrono::steady_clock::now(); // 记录结束时间
  diff = end - start; // 计算时间差
  verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

  assert(verify_1 == true);
}

template<typename FieldT, typename ppT, typename HType>
void comEq_tests() {
  size_t mu, n, col_num, repeat;
  double t_prove, t_verify;
  double t_prove_total, t_verify_total;
  double t_prove_ave, t_verify_ave;
  n = 1;
  col_num = 5;
  repeat  = 50;

  printf("\n-----------------------------------------------------------------------------------------------\n");
  if (std::is_same<libff::G1<ppT>, HType>::value == true) {
    printf("ℎ取自𝔾1时𝑝𝑝𝑐𝑜𝑚𝐸𝑞效率测试，重复次数%ld\n", repeat);
  } else {
    printf("ℎ取自𝔾2时𝑝𝑝𝑐𝑜𝑚𝐸𝑞效率测试，重复次数%ld\n", repeat);
  }
  printf("    行数    元素总数    证明生成时间(ms)    证明验证时间(ms)\n");
  for (mu = 0; mu < 8; mu++) {
    t_prove_total  = 0.0;
    t_verify_total = 0.0;
    for (size_t i = 0; i < repeat; i++) {
      if (std::is_same<libff::G1<ppT>, HType>::value == true) {
        test_single_h_in_G1<FieldT, ppT>(mu, n, col_num, t_prove, t_verify);
      } else {
        test_single_h_in_G2<FieldT, ppT>(mu, n, col_num, t_prove, t_verify);
      }

      t_prove_total  += t_prove;
      t_verify_total += t_verify;
    }
    t_prove_ave = t_prove_total / repeat;
    t_verify_ave = t_verify_total / repeat;

    size_t row_num = (size_t)std::pow(2, mu);
    size_t item_num = row_num * col_num;
    printf("%8ld %8ld %16.2lf %20.2lf\n", row_num, item_num, t_prove_ave, t_verify_ave);
    // printf("pp_comEq<FieldT, ppT>::prove : %.2f ms\n", t_prove_ave );
    // printf("pp_comEq<FieldT, ppT>::verify: %.2f ms\n", t_verify_ave);
  }
}

template<typename FieldT, typename ppT, typename HType>
void comEq_filling_tests() {
  size_t mu, n, col_num, repeat;
  double t_prove, t_verify;
  double t_prove_total, t_verify_total;
  double t_prove_ave, t_verify_ave;
  n = 1;
  col_num = 5;
  repeat  = 50;

  printf("\n-----------------------------------------------------------------------------------------------\n");
  if (std::is_same<libff::G1<ppT>, HType>::value == true) {
    printf("承诺输入数目固定效率测试：ℎ取自𝔾1，承诺数目固定为 5，重复次数%ld\n", repeat);
  } else {
    printf("承诺输入数目固定效率测试：ℎ取自𝔾2，承诺数目固定为 5，重复次数%ld\n", repeat);
  }
  printf("    行数    元素总数    证明生成时间(ms)    证明验证时间(ms)\n");
  for (mu = 0; mu < 8; mu++) {
    t_prove_total  = 0.0;
    t_verify_total = 0.0;
    for (size_t i = 0; i < repeat; i++) {
      if (std::is_same<libff::G1<ppT>, HType>::value == true) {
        test_single_filling_h_in_G1<FieldT, ppT>(mu, n, col_num, t_prove, t_verify);
      } else {
        test_single_filling_h_in_G2<FieldT, ppT>(mu, n, col_num, t_prove, t_verify);
      }

      t_prove_total  += t_prove;
      t_verify_total += t_verify;
    }
    t_prove_ave = t_prove_total / repeat;
    t_verify_ave = t_verify_total / repeat;

    size_t row_num = (size_t)std::pow(2, mu);
    size_t item_num = row_num * col_num;
    printf("%8ld %8ld %16.2lf %20.2lf\n", row_num, item_num, t_prove_ave, t_verify_ave);
    // printf("pp_comEq<FieldT, ppT>::prove : %.2f ms\n", t_prove_ave );
    // printf("pp_comEq<FieldT, ppT>::verify: %.2f ms\n", t_verify_ave);
  }
}

#endif