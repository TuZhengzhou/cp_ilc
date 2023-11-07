#ifndef BOOTLE_ILC_TCC
#define BOOTLE_ILC_TCC
#include <random>
#include <chrono>
#include <iostream>
#include <string>
#include <map>
#include "ILC.hpp"
#include "utils.hpp"
#include "libff/common/profiling.hpp"

using namespace std;

size_t COM_NUM = 1;
const size_t G1_TO_G1_TIMES_MIN = 10;
const size_t G1_TO_G1_TIMES_MAX = 100;
const size_t EXP_PER_COM = 2;
const size_t GATES_PER_EXP = (size_t)254 * 65;
const size_t MUL_GATES_PER_EXP = (size_t)254 * 30;
const size_t ADD_GATES_PER_EXP = (size_t)254 * 35;
const size_t WIRE_PER_GATE = 3;

const double ADD_PARTIAL = 6.0;
const double MUL_PARTIAL = 7.0;

//   ret["mu"] = mu;
//   ret[ARG::M]  = m;
//   ret[ARG::N_A] = n_A;
//   ret[ARG::N_M] = n_M;
//   ret[ARG::ROW_NUM] = row_num;
//   ret[ARG::COL_NUM] = col_num;
enum ARG {
    MU,
    M,
    N_A,
    N_M,
    ROW_NUM,
    COL_NUM
};


template <typename FieldT, typename ppT, typename HType>
bool pp_ILC<FieldT, ppT, HType>::get_random_source(const size_t add_num, const size_t multiply_num, FieldT& source_val, size_t& source_idx_in_V) {
    // if (add_num == 0 && multiply_num == 0)
    //     return false;
    assert(add_num > 0 || multiply_num > 0);

    size_t reletive_idx;
    size_t row, col;
    const size_t choose_from_3 = std::rand() % 3;

    // 如果能从加法矩阵里选, 先在加法矩阵里选保证有返回值
    if (add_num > 0) {
        reletive_idx = std::rand() % add_num;

        switch (choose_from_3)
        {
        case 0:
            source_val = this->A_.get_item(reletive_idx);
            break;
        case 1:
            source_val = this->B_.get_item(reletive_idx);
            break;
        default:
            source_val = this->C_.get_item(reletive_idx);
            break;
        } 
        source_idx_in_V = choose_from_3 * this->row_num_A_ * this->col_num_ + reletive_idx;
    } 

    /* 如果乘法有源, 并且没有加法源或者有加法源但选中了乘法源, 则从乘法矩阵里选 */
    if (multiply_num > 0 && (add_num == 0 || std::rand() % 2 == 1)) {
        reletive_idx = std::rand() % multiply_num;

        switch (choose_from_3)
        {
        case 0:
            source_val = this->D_.get_item(reletive_idx);
            break;
        case 1:
            source_val = this->E_.get_item(reletive_idx);
            break;
        default:
            source_val = this->F_.get_item(reletive_idx);
            break;
        } 
        source_idx_in_V = (3 * this->row_num_A_ * this->col_num_) + (choose_from_3 * this->row_num_M_ * this->col_num_) + reletive_idx;
    }
    return true;
}

template <typename FieldT, typename ppT, typename HType>
bool pp_ILC<FieldT, ppT, HType>::mark_copy(std::vector<std::vector<pos_marker> >& markers, const size_t target_idx, const size_t source_idx) {
    size_t row_s, col_s;
    std::tie(row_s, col_s) = to_tuple_pos(source_idx, this->col_num_);

    /* copy 自其他源, 则找到最初源的位置 */
    if (markers[row_s][col_s].is_copy_ == true) {
        std::tie(row_s, col_s) = markers[row_s][col_s].copy_from_;
    }

    /* 将目标位置 添加到源位置 的引用列表 */
    markers[row_s][col_s].add_copy(to_tuple_pos(target_idx, this->col_num_));
    
    size_t row_t, col_t;
    std::tie(row_t, col_t) = to_tuple_pos(target_idx, this->col_num_);
    /* 标记目标位置为引用 */
    markers[row_t][col_t].pos_       = to_tuple_pos(target_idx, this->col_num_);
    markers[row_t][col_t].is_copy_   = true;
    markers[row_t][col_t].copy_from_ = std::make_tuple(row_s, col_s);

    return true;
}

template <typename FieldT, typename ppT, typename HType>
pp_ILC<FieldT, ppT, HType>::pp_ILC(const size_t mu_A, const size_t mu_M, const size_t mu_V, const size_t n_A, const size_t n_M, \
                       const size_t n_V, const size_t col_num, const row_vector_matrix<FieldT>& V, \
                       const permutation<tuple_dim2_t, 2>& PI, const std::vector<size_t>& S, const bool com_input) {
    this->mu_A_ = mu_A;
    this->mu_M_ = mu_M;
    this->mu_V_ = mu_V;
    this->n_A_ = n_A;
    this->n_M_ = n_M;
    this->n_V_ = n_V;
    this->col_num_ = col_num;

    this->m_A_ = std::pow(2, this->mu_A_);
    this->m_M_ = std::pow(2, this->mu_M_);
    this->row_num_A_ = this->m_A_ * this->n_A_;
    this->row_num_M_ = this->m_M_ * this->n_M_;
    this->row_num_V_ = (this->row_num_A_ + this->row_num_M_) * 3;
    // printf("pp_ILC: row_num = %ld, col_num = %ld\n", this->row_num_V_, this->col_num_);

    this->PI_ = permutation<tuple_dim2_t, 2>(PI);
    this->S_ = std::vector<size_t>(S.begin(), S.end());

    this->V_ = row_vector_matrix<FieldT>(V);
    this->A_ = row_vector_matrix<FieldT>(this->col_num_);
    this->B_ = row_vector_matrix<FieldT>(this->col_num_);
    this->C_ = row_vector_matrix<FieldT>(this->col_num_);
    this->D_ = row_vector_matrix<FieldT>(this->col_num_);
    this->E_ = row_vector_matrix<FieldT>(this->col_num_);
    this->F_ = row_vector_matrix<FieldT>(this->col_num_);

    this->U_ = row_vector_matrix<FieldT>(this->col_num_);
    for (size_t idx: S) {
        assert(idx < this->row_num_V_);
        this->U_.add_row_vector(this->V_.get_row(idx));
    }

    size_t i;
    for (i = 0 * this->row_num_A_; i < 1 * this->row_num_A_; i++)     { this->A_.add_row_vector(this->V_.get_row(i)); }
    for (i = 1 * this->row_num_A_; i < 2 * this->row_num_A_; i++)     { this->B_.add_row_vector(this->V_.get_row(i)); }
    for (i = 2 * this->row_num_A_; i < 3 * this->row_num_A_; i++)     { this->C_.add_row_vector(this->V_.get_row(i)); }
    for (i = 3 * this->row_num_A_ + 0 * this->row_num_M_; i < 3 * this->row_num_A_ + 1 * this->row_num_M_; i++) { this->D_.add_row_vector(this->V_.get_row(i)); }
    for (i = 3 * this->row_num_A_ + 1 * this->row_num_M_; i < 3 * this->row_num_A_ + 2 * this->row_num_M_; i++) { this->E_.add_row_vector(this->V_.get_row(i)); }
    for (i = 3 * this->row_num_A_ + 2 * this->row_num_M_; i < 3 * this->row_num_A_ + 3 * this->row_num_M_; i++) { this->F_.add_row_vector(this->V_.get_row(i)); }

    this->eq_ = pp_eq<FieldT>(S.size(), this->col_num_, this->U_);
    this->sum_ = pp_sum<FieldT>(this->row_num_A_, this->col_num_, this->A_, this->B_, this->C_);
    this->prod_ = pp_prod<FieldT>(this->mu_M_, this->n_M_, this->col_num_, this->D_, this->E_, this->F_);
    this->perm_ = pp_perm<FieldT>(this->mu_V_, this->n_V_, this->col_num_, this->V_, this->V_, this->PI_);

    this->com_input_ = com_input;
    if (com_input == true) {
        libff::enter_block("prepare comEq");

        libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
        libff::G1<ppT> h_base = FieldT::random_element() * libff::G1<ppT>::G1_one;

        std::vector<FieldT> as1 = this->A_.get_row(0).get_all_items();
        std::vector<FieldT> as2 = this->D_.get_row(0).get_all_items();
        std::vector<FieldT> rs1 = row_vector<FieldT>::random(this->col_num_).get_all_items();
        std::vector<FieldT> rs2 = row_vector<FieldT>::random(this->col_num_).get_all_items();
        std::vector<CommitT> coms1 = std::vector<CommitT>(this->col_num_);
        std::vector<CommitT> coms2 = std::vector<CommitT>(this->col_num_);
        for (size_t i = 0; i < this->col_num_; i++) {
            coms1[i] = CommitT(as1[i] * g_base, rs1[i] * h_base);
            coms2[i] = CommitT(as2[i] * g_base, rs2[i] * h_base);
        }
        row_vector_matrix<FieldT> com_add (this->col_num_);
        row_vector_matrix<FieldT> com_mult(this->col_num_);

        com_add.add_row_vector (row_vector<FieldT>(std::vector<FieldT>(this->A_.get_row(0).get_all_items())));
        com_mult.add_row_vector(row_vector<FieldT>(std::vector<FieldT>(this->D_.get_row(0).get_all_items())));

        this->comEq_add_  = pp_comEq<FieldT, ppT, HType>(0, 1, this->col_num_, this->col_num_, com_add , rs1, g_base, h_base, coms1);
        this->comEq_mult_ = pp_comEq<FieldT, ppT, HType>(0, 1, this->col_num_, this->col_num_, com_mult, rs2, g_base, h_base, coms2);
        
        libff::leave_block("prepare comEq");
    }

}

template <typename FieldT, typename ppT, typename HType>
bool pp_ILC<FieldT, ppT, HType>::prove() {
    libff::enter_block("pp_ILC<FieldT, ppT>::prove()");
    this->prod_.prove();
    this->perm_.prove();
    
    // libff::inhibit_profiling_counters = false;
    libff::enter_block("pp_ILC<FieldT, ppT>::prove comEq");
    if (this->com_input_) {
        this->comEq_add_.prove();
        this->comEq_mult_.prove();
    }
    libff::leave_block("pp_ILC<FieldT, ppT>::prove comEq");
    // libff::inhibit_profiling_counters = true;
    libff::leave_block("pp_ILC<FieldT, ppT>::prove()");
    return true;
}

template <typename FieldT, typename ppT, typename HType>
bool pp_ILC<FieldT, ppT, HType>::verify(const bool output) {

    libff::enter_block("pp_ILC<FieldT, ppT>::verify()");
    row_vector_matrix<FieldT> pub_U = row_vector_matrix<FieldT>(this->col_num_);
    for (size_t idx: this->S_) {
        assert(idx < this->row_num_V_);
        pub_U.add_row_vector(this->V_.get_row(idx));
    }

    bool output_mid_result = false;
    bool check_result = true;
    check_result &= this->eq_.verify(pub_U, output_mid_result);
    check_result &= this->sum_.verify(row_vector<FieldT>::all_zero(this->col_num_), output_mid_result);
    check_result &= this->prod_.verify(output_mid_result);
    check_result &= this->perm_.verify(output_mid_result);

    // libff::inhibit_profiling_counters = false;
    libff::enter_block("pp_ILC<FieldT, ppT>::verify comEq");
    if (this->com_input_) {
        check_result &= this->comEq_add_.verify(output_mid_result);
        check_result &= this->comEq_mult_.verify(output_mid_result);
    }
    libff::leave_block("pp_ILC<FieldT, ppT>::verify comEq");
    // libff::inhibit_profiling_counters = true;

    libff::leave_block("pp_ILC<FieldT, ppT>::verify()");
    if (check_result == true) {
        if (output) {
            libff::print_indent(); printf("pp_ILC<FieldT, ppT>::verify() \033[32mpass\033[37m\n");
        }
        return true;
    }
    if (output) {
        libff::print_indent(); printf("pp_ILC<FieldT, ppT>::verify() \033[31mfail\033[37m\n");
    }
    return false;
}

inline size_t get_num_gate_before(size_t exp_iter) {
    return (pow(2, exp_iter) * (ADD_PARTIAL + MUL_PARTIAL)) * 90;
}

inline size_t get_num_gate_with_commit(size_t exp_iter, size_t num_commit) {
    return get_num_gate_before(exp_iter) + num_commit * EXP_PER_COM * (MUL_GATES_PER_EXP + ADD_GATES_PER_EXP);
}

map<size_t, size_t> generate_circuit_size_arguments(size_t exp_iter, size_t commit_num) {
  size_t m, mu, n_A, n_M, row_num, col_num;
  size_t num_gate_commit;
  double num_gate_before, num_gate, num_gate_apos, num_wire, sq, row_A, row_M, n;

  /** first cal the size of origin circuit */
  num_gate_before = get_num_gate_before(exp_iter);
  num_gate = get_num_gate_with_commit(exp_iter, commit_num);

  num_wire = num_gate * WIRE_PER_GATE;
  sq       = sqrt(num_wire);
  row_A    = sq * (ADD_PARTIAL) / (ADD_PARTIAL + MUL_PARTIAL);
  row_M    = sq * (MUL_PARTIAL) / (ADD_PARTIAL + MUL_PARTIAL);
  n     = sq / std::pow(2, exp_iter);
  n_A   = row_A / std::pow(2, exp_iter); 
  n_M   = row_M / std::pow(2, exp_iter);

  /** minimize nlogn + mn */
  solveOptimizationProblem(row_A, row_M, n_A, n_M, mu);
  m = (size_t)std::pow(2, mu);

  row_num = (n_A + n_M) * m;
  col_num = (size_t)std::ceil(num_gate / row_num);
  num_gate_apos = row_num * col_num;

  map<size_t, size_t> ret;
  ret[ARG::MU] = mu;
  ret[ARG::M]  = m;
  ret[ARG::N_A] = n_A;
  ret[ARG::N_M] = n_M;
  ret[ARG::ROW_NUM] = row_num;
  ret[ARG::COL_NUM] = col_num;

  return ret;
}

template <typename FieldT, typename ppT, typename HType>
pp_ILC<FieldT, ppT, HType> pp_ILC<FieldT, ppT, HType>::random(const size_t mu_A, const size_t mu_M, const size_t mu_V, \
                                                const size_t n_A, const size_t n_M, const size_t n_V, \
                                                const size_t col_num, const bool com_input) {
    size_t m_A, m_M, m_V;
    m_A = std::pow(2, mu_A);
    m_M = std::pow(2, mu_M);
    m_V = std::pow(2, mu_V);
    assert((m_A * n_A + m_M * n_M) * 3 == m_V * n_V);

    pp_ILC<FieldT, ppT, HType> result;
    result.mu_A_ = mu_A;    
    result.mu_M_ = mu_M;
    result.mu_V_ = mu_V;
    result.n_A_  = n_A;
    result.n_M_  = n_M;
    result.n_V_  = n_V;
    result.col_num_ = col_num;

    // libff::print_indent(); printf("row_num = \033[32m%ld\033[37m, col_num = \033[32m%ld\033[37m\n", (size_t)(std::pow(2, result.mu_V_) * result.n_V_), result.col_num_);


    result.m_A_ = std::pow(2, result.mu_A_);
    result.m_M_ = std::pow(2, result.mu_M_);
    result.row_num_A_ = result.m_A_ * result.n_A_;
    result.row_num_M_ = result.m_M_ * result.n_M_;
    result.row_num_V_ = (result.row_num_A_ + result.row_num_M_) * 3;

    const size_t input_row_A = 1;
    const size_t input_row_M = 1;
    result.A_ = row_vector_matrix<FieldT>::random(result.row_num_A_, result.col_num_);
    result.B_ = row_vector_matrix<FieldT>::random(result.row_num_A_, result.col_num_);
    result.D_ = row_vector_matrix<FieldT>::random(result.row_num_M_, result.col_num_);
    result.E_ = row_vector_matrix<FieldT>::random(result.row_num_M_, result.col_num_);

    // const size_t num_commit = result.col_num_;
    const size_t num_commit = std::min(COM_NUM, result.col_num_);
    /* 修改承诺输入为加法乘法各 num_commit 个, 对应行的其他位置填充 0 */
    for (size_t i = num_commit; i < result.col_num_; i++) {
        result.A_.set_item(i, FieldT::zero());
        result.D_.set_item(i, FieldT::zero());
    }
    result.C_ = result.A_ + result.B_;
    result.F_ = result.D_ * result.E_;

    /* 初始化 markers */
    std::vector<std::vector<pos_marker> > markers = std::vector<std::vector<pos_marker> >(result.row_num_V_, std::vector<pos_marker>(result.col_num_));
    for (size_t i = 0; i < result.row_num_V_; i++) {
        for (size_t j = 0; j < result.col_num_; j++) {
            markers[i][j] = pos_marker(std::make_tuple(i, j));
        }
    }

    bool coin;    // 决定加法还是乘法
    // size_t idx_A, idx_M;
    std::vector<size_t> idxs         = {input_row_A * result.col_num_, input_row_M * result.col_num_};
    std::vector<size_t> upper_bounds = {result.row_num_A_ * result.col_num_, result.row_num_M_ * result.col_num_};
    // const std::vector<size_t> bases  = {0, result.row_num_A_ * 3};
    // const std::vector<size_t> gaps   = {result.row_num_A_ * result.col_num_, result.row_num_M_ * result.col_num_};

    /* 生成随机电路 */
    while( idxs[0] < upper_bounds[0] || idxs[1] < upper_bounds[1]) {
        // std::cout << "??" << std::endl;
        if (idxs[0] < upper_bounds[0]) {
            coin = false;
        }
        if (idxs[1] < upper_bounds[1] && (idxs[0] >= upper_bounds[0] || std::rand() % 2)) {
            coin = true;
        }
        
        FieldT left_val, right_val;
        size_t left_source_idx_in_V, right_source_idx_in_V;
        result.get_random_source(idxs[0], idxs[1], left_val , left_source_idx_in_V );
        result.get_random_source(idxs[0], idxs[1], right_val, right_source_idx_in_V);
        
        size_t left_target_idx_in_V, right_target_idx_in_V;
        // size_t output_idx_in_V, row_o, col_o;
        if (coin == false) {    /* 加法 */
            result.A_.set_item(idxs[0], left_val);
            result.B_.set_item(idxs[0], right_val);
            result.C_.set_item(idxs[0], left_val + right_val);

            left_target_idx_in_V  = 0 + idxs[0];
            right_target_idx_in_V = result.row_num_A_ * result.col_num_ + idxs[0];
            result.mark_copy(markers, left_target_idx_in_V , left_source_idx_in_V );
            result.mark_copy(markers, right_target_idx_in_V, right_source_idx_in_V);

            // output_idx_in_V = result.row_num_A_ * result.col_num_ * 2 + idxs[0];
            // std::tie(row_o, col_o) = to_tuple_pos(output_idx_in_V);
            // markers[row_o][col_o].is_copy_ = false;
        } else {                /* 乘法 */
            result.D_.set_item(idxs[1], left_val);
            result.E_.set_item(idxs[1], right_val);
            result.F_.set_item(idxs[1], left_val * right_val);

            left_target_idx_in_V  = result.row_num_A_ * result.col_num_ * 3 + 0 + idxs[1];
            right_target_idx_in_V = result.row_num_A_ * result.col_num_ * 3 + result.row_num_M_ * result.col_num_ + idxs[1];
            result.mark_copy(markers, left_target_idx_in_V , left_source_idx_in_V );
            result.mark_copy(markers, right_target_idx_in_V, right_source_idx_in_V);
        }
        idxs[coin] += 1;
    }
    
    /* 把矩阵转移到 V_ */
    result.V_ = row_vector_matrix<FieldT>(result.col_num_);
    for (size_t i = 0; i < result.row_num_A_; i++) result.V_.add_row_vector(result.A_.get_row(i));
    for (size_t i = 0; i < result.row_num_A_; i++) result.V_.add_row_vector(result.B_.get_row(i));
    for (size_t i = 0; i < result.row_num_A_; i++) result.V_.add_row_vector(result.C_.get_row(i));
    for (size_t i = 0; i < result.row_num_M_; i++) result.V_.add_row_vector(result.D_.get_row(i));
    for (size_t i = 0; i < result.row_num_M_; i++) result.V_.add_row_vector(result.E_.get_row(i));
    for (size_t i = 0; i < result.row_num_M_; i++) result.V_.add_row_vector(result.F_.get_row(i));
    // result.V_.print();

    /* 求解排列 PI */
    std::vector<std::vector<tuple_dim2_t> > PI;
    for (size_t i = 0; i < result.row_num_V_; i++) {
        for (size_t j = 0; j < result.col_num_; j++) {
            if (markers[i][j].is_copy_ == false && markers[i][j].copies_.size() >= 1) {

                // std::cout << "cycle: " << std::endl;
                std::vector<tuple_dim2_t> cycle_tmp = {markers[i][j].pos_};
                // std::cout << "[" << std::get<0>(markers[i][j].pos_) << ", " << std::get<1>(markers[i][j].pos_) << "], ";
                for (size_t k = 0; k < markers[i][j].copies_.size(); k++) {
                    cycle_tmp.push_back(markers[i][j].copies_[k]);
                    // std::cout << "[" << std::get<0>(markers[i][j].copies_[k]) << ", " << std::get<1>(markers[i][j].copies_[k]) << "], ";
                }
                // std::cout << std::endl;

                PI.push_back(cycle_tmp);

            }
        }
    }
    result.PI_ = permutation<tuple_dim2_t, 2>(PI);

    /* 随机设置指派值 */
    result.U_ = row_vector_matrix<FieldT>(result.col_num_);
    size_t S_len = std::max(1ul, result.row_num_V_ / 6ul);
    for (size_t i = 0; i < S_len; i++) {
        result.S_.push_back(std::rand() % result.row_num_V_);
        result.U_.add_row_vector(result.V_.get_row(result.S_.back()));
    }

    /* 初始化子协议 */
    result.eq_ = pp_eq<FieldT>(result.S_.size(), result.col_num_, result.U_);
    result.sum_ = pp_sum<FieldT>(result.row_num_A_, result.col_num_, result.A_, result.B_, result.C_);
    result.prod_ = pp_prod<FieldT>(result.mu_M_, result.n_M_, result.col_num_, result.D_, result.E_, result.F_);
    result.perm_ = pp_perm<FieldT>(result.mu_V_, result.n_V_, result.col_num_, result.V_, result.V_, result.PI_);

    result.com_input_ = com_input;
    if (result.com_input_) {
        libff::enter_block("prepare comEq");

        libff::G1<ppT> g_base = FieldT::random_element() * libff::G1<ppT>::G1_one;
        libff::G1<ppT> h_base = FieldT::random_element() * libff::G1<ppT>::G1_one;

        std::vector<FieldT> as1 = result.A_.get_row(0).get_all_items();
        std::vector<FieldT> as2 = result.D_.get_row(0).get_all_items();
        std::vector<FieldT> rs1 = row_vector<FieldT>::random(num_commit).get_all_items();
        std::vector<FieldT> rs2 = row_vector<FieldT>::random(num_commit).get_all_items();
        std::vector<CommitT> coms1 = std::vector<CommitT>(num_commit);
        std::vector<CommitT> coms2 = std::vector<CommitT>(num_commit);
        for (size_t i = 0; i < num_commit; i++) {
            coms1[i] = CommitT(as1[i] * g_base, rs1[i] * h_base);
            coms2[i] = CommitT(as2[i] * g_base, rs2[i] * h_base);
        }
        row_vector_matrix<FieldT> com_add (result.col_num_);
        row_vector_matrix<FieldT> com_mult(result.col_num_);

        com_add.add_row_vector (row_vector<FieldT>(std::vector<FieldT>(result.A_.get_row(0).get_all_items())));
        com_mult.add_row_vector(row_vector<FieldT>(std::vector<FieldT>(result.D_.get_row(0).get_all_items())));

        result.comEq_add_  = pp_comEq<FieldT, ppT, HType>(0, 1, result.col_num_, num_commit, com_add , rs1, g_base, h_base, coms1);
        result.comEq_mult_ = pp_comEq<FieldT, ppT, HType>(0, 1, result.col_num_, num_commit, com_mult, rs2, g_base, h_base, coms2);
        
        libff::leave_block("prepare comEq");
    }

    return result;
}

template <typename FieldT, typename ppT, typename HType>
bool ILC_test_origin(const size_t exp_iter, std::map<size_t, size_t> circuit_arguments) {
    size_t mu, m, n_A, n_M, n_V, row_num, col_num;
    size_t row_A, row_M;  /* 加法行向量个数, 乘法行向量个数, 行向量总数 */
  
    permutation<tuple_dim2_t, 2> PI;
    row_vector_matrix<FieldT> U;
  
    row_vector_matrix<FieldT> V;
    row_vector_matrix<FieldT> A, B, C, D, E, F;
    
    mu = circuit_arguments[ARG::MU];
    m  = circuit_arguments[ARG::M];
    n_A = circuit_arguments[ARG::N_A];
    n_M = circuit_arguments[ARG::N_M];
    row_num = circuit_arguments[ARG::ROW_NUM];
    col_num = circuit_arguments[ARG::COL_NUM];

    row_A = m * n_A;
    row_M = m * n_M;
    n_V   = 3 * (n_A + n_M);

    A = row_vector_matrix<FieldT>::random(row_A, col_num);
    B = row_vector_matrix<FieldT>::random(row_A, col_num);
    C = A + B;

    D = row_vector_matrix<FieldT>::random(row_M, col_num);
    E = row_vector_matrix<FieldT>::random(row_M, col_num);
    F = D * E;

    V = row_vector_matrix<FieldT>(col_num);
    for (size_t i = 0; i < row_A; i++) V.add_row_vector(A.get_row(i));
    for (size_t i = 0; i < row_A; i++) V.add_row_vector(B.get_row(i));
    for (size_t i = 0; i < row_A; i++) V.add_row_vector(C.get_row(i));
    for (size_t i = 0; i < row_M; i++) V.add_row_vector(D.get_row(i));
    for (size_t i = 0; i < row_M; i++) V.add_row_vector(E.get_row(i));
    for (size_t i = 0; i < row_M; i++) V.add_row_vector(F.get_row(i));

    size_t cycle_len = 1;   /* 设置为 1, 即不进行排列变换 */
    size_t cycle_num = row_num / 2;
    std::vector<size_t> dim_limits = {row_num, col_num};
    PI = permutation<tuple_dim2_t, 2>::random_permutation(cycle_len, cycle_num, 2, dim_limits);

    // // 将 A 矩阵中在同一 cycle 位置的元素设置为相同值
    // size_t row, col;
    // for (size_t i = 0; i < cycle_num; i++) {
    //     FieldT val = FieldT::random_element();
    //     const std::vector<tuple_dim2_t> cyc_vec = PI.get_cycle(i).get_content_vector();
    //     for (size_t j = 0; j < cycle_len; j++) {
    //         std::tie(row, col) = cyc_vec[j];
    //         A.set_item(row, col, val);
    //     }
    // }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, row_num);
    std::vector<size_t> S;
    for (size_t i = 0; i < row_num / 3; i++) {
        S.emplace_back(dis(gen) % row_num);
    } 

    // pp_ILC<FieldT, ppT, HType> ILC = pp_ILC<FieldT, ppT, HType>(mu_A, mu_M, mu_V, n_A, n_M, n_V, col_num, V, PI, S, false);
    pp_ILC<FieldT, ppT, HType> ILC = pp_ILC<FieldT, ppT, HType>(mu, mu, mu, n_A, n_M, n_V, col_num, V, PI, S, false);

    libff::enter_block("0 - NO-CP TOTAL");

    libff::enter_block("0 - NO-CP PROVE");
    ILC.prove();
    libff::leave_block("0 - NO-CP PROVE");

    libff::enter_block("0 - NO-CP VERIFY");
    ILC.verify();
    libff::leave_block("0 - NO-CP VERIFY");

    libff::leave_block("0 - NO-CP TOTAL");


    return true;
    
}

template <typename FieldT, typename ppT, typename HType>
bool ILC_test_with_com_input(const size_t exp_iter, std::map<size_t, size_t> circuit_arguments) {
    size_t mu, m, col_num;
    size_t row_A, n_A;
    size_t row_M, n_M;
    size_t row_num, n_V;

    mu = circuit_arguments[ARG::MU];
    m  = circuit_arguments[ARG::M];
    col_num = circuit_arguments[ARG::COL_NUM];

    n_A = circuit_arguments[ARG::N_A];
    n_M = circuit_arguments[ARG::N_M];
    n_V   = 3 * (n_A + n_M);
    row_A = m * n_A;
    row_M = m * n_M;
    row_num = circuit_arguments[ARG::ROW_NUM];

    std::chrono::_V2::steady_clock::time_point start, end;
    std::chrono::duration<int64_t, std::nano>  diff;

    pp_ILC<FieldT, ppT, HType> ILC_random = pp_ILC<FieldT, ppT, HType>::random(mu, mu, mu, n_A, n_M, n_V, col_num, true);

    libff::enter_block("1 - CP TOTAL");

    libff::enter_block("1 - CP PROVE");
    ILC_random.prove();
    libff::leave_block("1 - CP PROVE");

    libff::enter_block("1 - CP VERIFY");
    ILC_random.verify();
    libff::leave_block("1 - CP VERIFY");

    libff::leave_block("1 - CP TOTAL");

    return true;
}

string effeciency_print_no_cp() {
  libff::inhibit_profiling_info = false;
  libff::print_cumulative_time_entry("0 - NO-CP TOTAL");
  libff::print_cumulative_time_entry("0 - NO-CP PROVE");
  libff::print_cumulative_time_entry("0 - NO-CP VERIFY");
  // libff::print_cumulative_op_counts();
  // libff::print_time("0 - NO-CP TOTAL");
  // libff::print_time("0 - NO-CP PROVE");
  // libff::print_time("0 - NO-CP VERIFY");
  libff::inhibit_profiling_info = true;

  return string("place holder");
}

string effeciency_print_cp() {
  libff::inhibit_profiling_info = false;
  libff::print_cumulative_time_entry("1 - CP TOTAL");
  libff::print_cumulative_time_entry("1 - CP PROVE");
  libff::print_cumulative_time_entry("1 - CP VERIFY");
  // libff::print_cumulative_op_counts();
  // libff::print_time("1 - CP TOTAL");
  // libff::print_time("1 - CP PROVE");
  // libff::print_time("1 - CP VERIFY");
  libff::inhibit_profiling_info = true;

  return string("place holder");
}
/* 
  承诺数量固定为 1，电路变大
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_1() {

  size_t exp_iter, repeat, num_commit;
  size_t num_gate_before, num_gate_after, num_gate_apos;
  std::map<size_t, size_t> circuit_arguments;

  repeat = 10;
  num_commit = 1;
  printf("\n-----------------------------------------------------------------------------------------------\n");
  printf("情形 1 下对比实验结果， 重复次数为 %ld\n", repeat);

  for (exp_iter = 1; exp_iter < 6; exp_iter++) {

    printf("iter = %lu\n", exp_iter);

    num_gate_before = get_num_gate_before(exp_iter);
    num_gate_after  = get_num_gate_with_commit(exp_iter, num_commit);
    printf("num_commit = %8lu num_gate_before = %8lu num_gate_after = %8lu\n", num_commit, num_gate_before, num_gate_after);

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, num_commit);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_origin<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_no_cp();

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, 0);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_with_com_input<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_cp();
  }
}

/* 
  电路不变，承诺数量倍数增加
  原电路门数量固定为
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_2() {

  size_t exp_iter, repeat, num_commit;
  size_t num_gate_before, num_gate_after, num_gate_apos;
  std::map<size_t, size_t> circuit_arguments;

  repeat = 10;
  exp_iter = 4;
  printf("\n-----------------------------------------------------------------------------------------------\n");
  printf("情形 2 下对比实验结果， 重复次数为 %ld\n", repeat);

  for (num_commit = 1; num_commit <= 32; num_commit *= 2) {

    printf("iter = %lu\n", exp_iter);

    num_gate_before = get_num_gate_before(exp_iter);
    num_gate_after  = get_num_gate_with_commit(exp_iter, num_commit);
    printf("num_commit = %8lu num_gate_before = %8lu num_gate_after = %8lu\n", num_commit, num_gate_before, num_gate_after);

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, num_commit);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_origin<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_no_cp();

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, 0);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_with_com_input<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_cp();
  }
}

/*
  电路和承诺数量比值不变（行数不变，承诺数量与列数相同）
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_3() {

  size_t exp_iter, repeat, num_commit;
  size_t num_gate_before, num_gate_after, num_gate_apos;
  std::map<size_t, size_t> circuit_arguments;

  repeat = 10;
  exp_iter = 2;
  printf("\n-----------------------------------------------------------------------------------------------\n");
  printf("情形 3 下对比实验结果， 重复次数为 %ld\n", repeat);

  for (exp_iter = 1; exp_iter <= 5; exp_iter++) {

    num_commit = pow(2, exp_iter);

    printf("iter = %lu\n", exp_iter);

    

    num_gate_before = get_num_gate_before(exp_iter);
    num_gate_after  = get_num_gate_with_commit(exp_iter, num_commit);
    printf("num_commit = %8lu num_gate_before = %8lu num_gate_after = %8lu\n", num_commit, num_gate_before, num_gate_after);

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, num_commit);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_origin<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_no_cp();

    libff::clear_profiling_counters();
    circuit_arguments = generate_circuit_size_arguments(exp_iter, 0);
    for (size_t i = 0; i < repeat; i++) {
      ILC_test_with_com_input<FieldT, ppT, HType>(exp_iter, circuit_arguments);
    }
    effeciency_print_cp();
  }
}
#endif