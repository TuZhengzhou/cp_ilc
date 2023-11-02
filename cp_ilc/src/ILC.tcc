#ifndef BOOTLE_ILC_TCC
#define BOOTLE_ILC_TCC
#include <random>
#include <chrono>
#include "ILC.hpp"
#include "libff/common/profiling.hpp"

size_t COM_NUM = 1;
const size_t G1_TO_G1_TIMES_MIN = 10;
const size_t G1_TO_G1_TIMES_MAX = 100;

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

    // const size_t com_num = result.col_num_;
    const size_t com_num = std::min(COM_NUM, result.col_num_);
    /* 修改承诺输入为加法乘法各 com_num 个, 对应行的其他位置填充 0 */
    for (size_t i = com_num; i < result.col_num_; i++) {
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
        std::vector<FieldT> rs1 = row_vector<FieldT>::random(com_num).get_all_items();
        std::vector<FieldT> rs2 = row_vector<FieldT>::random(com_num).get_all_items();
        std::vector<CommitT> coms1 = std::vector<CommitT>(com_num);
        std::vector<CommitT> coms2 = std::vector<CommitT>(com_num);
        for (size_t i = 0; i < com_num; i++) {
            coms1[i] = CommitT(as1[i] * g_base, rs1[i] * h_base);
            coms2[i] = CommitT(as2[i] * g_base, rs2[i] * h_base);
        }
        row_vector_matrix<FieldT> com_add (result.col_num_);
        row_vector_matrix<FieldT> com_mult(result.col_num_);

        com_add.add_row_vector (row_vector<FieldT>(std::vector<FieldT>(result.A_.get_row(0).get_all_items())));
        com_mult.add_row_vector(row_vector<FieldT>(std::vector<FieldT>(result.D_.get_row(0).get_all_items())));

        result.comEq_add_  = pp_comEq<FieldT, ppT, HType>(0, 1, result.col_num_, com_num, com_add , rs1, g_base, h_base, coms1);
        result.comEq_mult_ = pp_comEq<FieldT, ppT, HType>(0, 1, result.col_num_, com_num, com_mult, rs2, g_base, h_base, coms2);
        
        libff::leave_block("prepare comEq");
    }

    return result;
}

template <typename FieldT, typename ppT, typename HType>
bool ILC_test_origin(const size_t mu, const size_t times_of_col, double& prove_time, double& verifiy_time) {
    size_t mu_A, mu_M, mu_V, n_A, n_M, n_V, row_num, col_num;
    size_t m_A, m_M, m;  /* 加法行向量个数, 乘法行向量个数, 行向量总数 */
  
    permutation<tuple_dim2_t, 2> PI;
    row_vector_matrix<FieldT> U;
  
    row_vector_matrix<FieldT> V;
    row_vector_matrix<FieldT> A, B, C, D, E, F;

    /* 加法-乘法比例与 椭圆曲线上点乘的加法-乘法比例相同，为 16: 2 */
    mu_A = mu;
    n_A  = 16;  /* 2^3 * 16 = 128 */
    mu_M = mu;
    n_M  = 2;   /* 2^3 * 2 = 16 */
    mu_V = mu;
    n_V  = 18 * 3;  /* 2^3 * 54 = 432 */
    row_num = 54 * std::pow(2, mu);
    col_num = 50 * std::pow(2, mu) * times_of_col + (size_t)std::ceil(FieldT::size_in_bits() * 18 * (1 + 1) * COM_NUM * 3 / row_num);
    m_A = std::pow(2, mu_A) * n_A;
    m_M = std::pow(2, mu_M) * n_M;
    m = m_A +m_M;

    A = row_vector_matrix<FieldT>::random(m_A, col_num);
    B = row_vector_matrix<FieldT>::random(m_A, col_num);
    C = A + B;

    D = row_vector_matrix<FieldT>::random(m_M, col_num);
    E = row_vector_matrix<FieldT>::random(m_M, col_num);
    F = D * E;

    V = row_vector_matrix<FieldT>(col_num);
    for (size_t i = 0; i < m_A; i++) V.add_row_vector(A.get_row(i));
    for (size_t i = 0; i < m_A; i++) V.add_row_vector(B.get_row(i));
    for (size_t i = 0; i < m_A; i++) V.add_row_vector(C.get_row(i));
    for (size_t i = 0; i < m_M; i++) V.add_row_vector(D.get_row(i));
    for (size_t i = 0; i < m_M; i++) V.add_row_vector(E.get_row(i));
    for (size_t i = 0; i < m_M; i++) V.add_row_vector(F.get_row(i));

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
    // std::cout << S << std::endl;
    // std::cout << "V_row: " << V.get_row_num() << ", V_col: " << V.get_column_num() << std::endl;

    pp_ILC<FieldT, ppT, HType> ILC = pp_ILC<FieldT, ppT, HType>(mu_A, mu_M, mu_V, n_A, n_M, n_V, col_num, V, PI, S, false);

    std::chrono::_V2::steady_clock::time_point start, end;
    std::chrono::duration<int64_t, std::nano>  diff;

    start = std::chrono::steady_clock::now(); // 记录开始时间
    ILC.prove();
    end = std::chrono::steady_clock::now(); // 记录结束时间
    diff = end - start; // 计算时间差
    prove_time = std::chrono::duration<double, std::milli>(diff).count();
    
    start = std::chrono::steady_clock::now(); // 记录开始时间
    ILC.verify();
    end = std::chrono::steady_clock::now(); // 记录结束时间
    diff = end - start; // 计算时间差
    verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

    return true;
    
}

template <typename FieldT, typename ppT, typename HType>
bool ILC_test_with_com_input(const size_t mu, const size_t times_of_col, double& prove_time, double& verifiy_time) {
    size_t mu_A, mu_M, mu_V, n_A, n_M, n_V, row_num, col_num;

    // mu_A = 4;
    // n_A  = 3;  /* 2^4 * 3 = 48 */
    // mu_M = 4;
    // n_M  = 9;  /* 2^4 * 9 = 144 */
    // mu_V = 6;
    // n_V  = 9;  /* 2^6 * 3 = 192 */
    // // row_num = 576;
    // col_num = 192;


    /* 加法-乘法比例与 椭圆曲线上点乘的加法-乘法比例相同，为 16: 2 */
    // mu_A = 3;
    // n_A  = 16;  /* 2^3 * 16 = 128 */
    // mu_M = 3;
    // n_M  = 2;   /* 2^3 * 2 = 16 */
    // mu_V = 3;
    // n_V  = 18 * 3;  /* 2^3 * 54 = 432 */
    // // row_num = 432;
    // col_num = 192;

    mu_A = mu;
    n_A  = 16;  /* 2^3 * 16 = 128 */
    mu_M = mu;
    n_M  = 2;   /* 2^3 * 2 = 16 */
    mu_V = mu;
    n_V  = 18 * 3;  /* 2^0 * 54 = 54 */
    // row_num = 54;
    col_num = 50 * std::pow(2, mu)  * times_of_col;

    std::chrono::_V2::steady_clock::time_point start, end;
    std::chrono::duration<int64_t, std::nano>  diff;

    start = std::chrono::steady_clock::now(); // 记录开始时间
    pp_ILC<FieldT, ppT, HType> ILC_random = pp_ILC<FieldT, ppT, HType>::random(mu_A, mu_M, mu_V, n_A, n_M, n_V, col_num, true);
    end = std::chrono::steady_clock::now(); // 记录结束时间
    diff = end - start; // 计算时间差
    // std::cout << "pp_ILC<FieldT, ppT>::random Time taken: " 
    //           << std::chrono::duration<double, std::micro>(diff).count() 
    //           << " us" << std::endl; // 输出时间差（单位为毫秒）

    start = std::chrono::steady_clock::now(); // 记录开始时间
    ILC_random.prove();
    end = std::chrono::steady_clock::now(); // 记录结束时间
    diff = end - start; // 计算时间差
    prove_time = std::chrono::duration<double, std::milli>(diff).count();
    
    start = std::chrono::steady_clock::now(); // 记录开始时间
    ILC_random.verify();
    end = std::chrono::steady_clock::now(); // 记录结束时间
    diff = end - start; // 计算时间差
    verifiy_time = std::chrono::duration<double, std::milli>(diff).count();

    return true;
}

/* 
  承诺数量固定为 1，电路变大（行数固定，增加列数）
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_1() {

    size_t mu, times_of_col, repeat;
    double t1_prove, t1_verify, t1_prove_total, t1_verify_total;
    double t2_prove, t2_verify, t2_prove_total, t2_verify_total;

    size_t row_num, col_before, col_after;
    double gates_before, gates_after;

    repeat = 1;
    mu = 1;

    printf("\n-----------------------------------------------------------------------------------------------\n");
    printf("情形 1 下对比实验结果， 重复次数为 %ld\n", repeat);
    printf("承诺输入数量   原电路门数量(× 1000)  电路门总数(× 1000)  新协议证明时间(ms) 朴素方案证明时间(ms)  新协议验证时间(ms)  朴素方案验证时间(ms)\n");


    for (times_of_col = 1; times_of_col < 5; times_of_col++) {
        COM_NUM = 1;
        t1_prove_total  = 0.0;
        t1_verify_total = 0.0;
        t2_prove_total  = 0.0;
        t2_verify_total = 0.0;
        for (size_t i = 0; i < repeat; i++) {
            ILC_test_origin<FieldT, ppT, HType>(mu, times_of_col, t1_prove, t1_verify);
            ILC_test_with_com_input<FieldT, ppT, HType>(mu, times_of_col, t2_prove, t2_verify);
            t1_prove_total  += t1_prove;
            t1_verify_total += t1_verify;
            t2_prove_total  += t2_prove;
            t2_verify_total += t2_verify;
        }
        t1_prove_total  /= repeat;
        t1_verify_total /= repeat;
        t2_prove_total  /= repeat;
        t2_verify_total /= repeat;

        row_num = 54 * std::pow(2, mu);
        col_before = 50 * std::pow(2, mu) * times_of_col;
        col_after  = col_before + (size_t)std::ceil(FieldT::size_in_bits() * 18 * (1 + 1) * COM_NUM * 3 / row_num);

        gates_before = row_num * col_before / 3.0 / 1000.0;
        gates_after  = row_num * col_after  / 3.0 / 1000.0;
        printf("%8ld %16.2lf %20.2lf %20.2lf %20.2lf %20.2lf %20.2lf\n", COM_NUM, gates_before, gates_after, t2_prove_total, t1_prove_total, t2_verify_total, t1_verify_total);
        // printf("ILC_test_origin::             prove : %.2f ms, verify: %.2f ms\n", t1_prove_total, t1_verify_total);
        // printf("ILC_test_with_com_input::     prove : %.2f ms, verify: %.2f ms\n", t2_prove_total, t2_verify_total);
    }
}

/* 
  电路不变，承诺数量倍数增加
  原电路门数量固定为
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_2() {

    size_t mu, times_of_col, com_num, repeat;
    double t1_prove, t1_verify, t1_prove_total, t1_verify_total;
    double t2_prove, t2_verify, t2_prove_total, t2_verify_total;

    size_t row_num, col_before, col_after;
    double gates_before, gates_after;

    repeat = 1;
    mu = 1;
    times_of_col = 1;

    printf("\n-----------------------------------------------------------------------------------------------\n");
    printf("情形 2 下对比实验结果， 重复次数为 %ld\n", repeat);
    printf("承诺输入数量   原电路门数量(× 1000)  电路门总数(× 1000)  新协议证明时间(ms) 朴素方案证明时间(ms)  新协议验证时间(ms)  朴素方案验证时间(ms)\n");
    
    for (com_num = 1; com_num < 64; com_num *= 2) {
        COM_NUM = com_num;
        t1_prove_total  = 0.0;
        t1_verify_total = 0.0;
        t2_prove_total  = 0.0;
        t2_verify_total = 0.0;
        for (size_t i = 0; i < repeat; i++) {
            ILC_test_origin<FieldT, ppT, HType>(1, 1, t1_prove, t1_verify);
            ILC_test_with_com_input<FieldT, ppT, HType>(1, 1, t2_prove, t2_verify);
            t1_prove_total  += t1_prove;
            t1_verify_total += t1_verify;
            t2_prove_total  += t2_prove;
            t2_verify_total += t2_verify;
        }
        t1_prove_total  /= repeat;
        t1_verify_total /= repeat;
        t2_prove_total  /= repeat;
        t2_verify_total /= repeat;

        row_num = 54 * std::pow(2, mu);
        col_before = 50 * std::pow(2, mu) * times_of_col;
        col_after  = col_before + (size_t)std::ceil(FieldT::size_in_bits() * 18 * (1 + 1) * COM_NUM * 3 / row_num);

        gates_before = row_num * col_before / 3.0 / 1000.0;
        gates_after  = row_num * col_after  / 3.0 / 1000.0;
        printf("%8ld %16.2lf %20.2lf %20.2lf %20.2lf %20.2lf %20.2lf\n", COM_NUM, gates_before, gates_after, t2_prove_total, t1_prove_total, t2_verify_total, t1_verify_total);
        // printf("%8ld %16ld %20ld %20.2lf %20.2lf %20.2lf %20.2lf\n", COM_NUM, 0, 0, t2_prove_total, t1_prove_total, t2_verify_total, t1_verify_total);
        // printf("ILC_test_origin::             prove : %.2f ms, verify: %.2f ms\n", t1_prove_total, t1_verify_total);
        // printf("ILC_test_with_com_input::     prove : %.2f ms, verify: %.2f ms\n", t2_prove_total, t2_verify_total);
    }
}

/*
  电路和承诺数量比值不变（行数不变，承诺数量与列数相同）
*/
template <typename FieldT, typename ppT, typename HType>
void ILC_test_compare_3() {

    size_t mu, times_of_col, com_num, repeat;
    double t1_prove, t1_verify, t1_prove_total, t1_verify_total;
    double t2_prove, t2_verify, t2_prove_total, t2_verify_total;

    size_t row_num, col_before, col_after;
    double gates_before, gates_after;

    repeat = 1;
    mu = 1;

    printf("\n-----------------------------------------------------------------------------------------------\n");
    printf("情形 3 下对比实验结果， 重复次数为 %ld\n", repeat);
    printf("承诺输入数量   原电路门数量(× 1000)  电路门总数(× 1000)  新协议证明时间(ms) 朴素方案证明时间(ms)  新协议验证时间(ms)  朴素方案验证时间(ms)\n");

    for (com_num = 1; com_num < 6; com_num++) {
        COM_NUM = com_num * 5;
        t1_prove_total  = 0.0;
        t1_verify_total = 0.0;
        t2_prove_total  = 0.0;
        t2_verify_total = 0.0;
        for (size_t i = 0; i < repeat; i++) {
            ILC_test_origin<FieldT, ppT, HType>(mu, com_num, t1_prove, t1_verify);
            ILC_test_with_com_input<FieldT, ppT, HType>(mu, com_num, t2_prove, t2_verify);
            t1_prove_total  += t1_prove;
            t1_verify_total += t1_verify;
            t2_prove_total  += t2_prove;
            t2_verify_total += t2_verify;
        }
        t1_prove_total  /= repeat;
        t1_verify_total /= repeat;
        t2_prove_total  /= repeat;
        t2_verify_total /= repeat;

        times_of_col = com_num;
        row_num = 54 * std::pow(2, mu);
        col_before = 50 * std::pow(2, mu) * times_of_col;
        col_after  = col_before + (size_t)std::ceil(FieldT::size_in_bits() * 18 * (1 + 1) * COM_NUM * 3 / row_num);

        gates_before = row_num * col_before / 3.0 / 1000.0;
        gates_after  = row_num * col_after  / 3.0 / 1000.0;
        printf("%8ld %16.2lf %20.2lf %20.2lf %20.2lf %20.2lf %20.2lf\n", COM_NUM, gates_before, gates_after, t2_prove_total, t1_prove_total, t2_verify_total, t1_verify_total);
        // printf("%8ld %16ld %20ld %20.2lf %20.2lf %20.2lf %20.2lf\n", COM_NUM, 0, 0, t2_prove_total, t1_prove_total, t2_verify_total, t1_verify_total);
        // printf("ILC_test_origin::             prove : %.2f ms, verify: %.2f ms\n", t1_prove_total, t1_verify_total);
        // printf("ILC_test_with_com_input::     prove : %.2f ms, verify: %.2f ms\n", t2_prove_total, t2_verify_total);
    }
}
#endif