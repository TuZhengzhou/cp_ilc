#ifndef BOOTLE_SAME_PROD_TCC
#define BOOTLE_SAME_PROD_TCC
#include "prod.hpp"
#include "prod.tcc"
#include "shift.hpp"
#include "same_prod.hpp"

template<typename FieldT>
bool pp_same_prod<FieldT>::is_satisfy() const {
    if (this->A_.partial_products(true).get_item(this->m_ * this->n_ - 1, this->col_num_ - 1)
        == this->B_.partial_products(true).get_item(this->m_ * this->n_ - 1, this->col_num_ - 1)
    ) {
        printf("pp_same_prod<FieldT>::is_satisfy() \033[032mpass\033[037m\n");
        return true;
    }
    printf("pp_same_prod<FieldT>::is_satisfy() \033[031mfail\033[037m\n");
    return false;
}

template<typename FieldT>
pp_same_prod<FieldT>::pp_same_prod(const size_t mu, const size_t n, const size_t k, const row_vector_matrix<FieldT>& A, const row_vector_matrix<FieldT>& B) {
    this->mu_ = mu;
    this->n_  = n;
    this->col_num_ = k;
    this->m_  = std::pow(2, this->mu_);

    this->A_  = row_vector_matrix<FieldT>(A);
    this->A1_ = this->A_.partial_products(false);
    this->A2_ = this->A1_ * this->A_;

    this->B_  = row_vector_matrix<FieldT>(B);
    this->B1_ = this->B_.partial_products(false);
    this->B2_  = this->B1_ * this->B_;

    this->pp_prod_A_ = pp_prod<FieldT>(this->mu_, this->n_, this->col_num_, this->A_, this->A1_, this->A2_);
    this->pp_prod_B_ = pp_prod<FieldT>(this->mu_, this->n_, this->col_num_, this->B_, this->B1_, this->B2_);
    this->pp_shift_AB_ = pp_shift<FieldT>(this->col_num_, this->mu_, this->n_, this->A1_, this->A2_, this->B1_, this->B2_);
}

template<typename FieldT>
bool pp_same_prod<FieldT>::prove(const std::map<std::string, FieldT> &mid_challenges) {
    this->pp_prod_A_.prove(mid_challenges.at("pp_prod_A_y"), mid_challenges.at("pp_prod_A_x0"));
    this->pp_prod_B_.prove(mid_challenges.at("pp_prod_B_y"), mid_challenges.at("pp_prod_B_x0"));
    this->pp_shift_AB_.prove(mid_challenges.at("pp_shift_AB_y"), mid_challenges.at("pp_shift_AB_x0"));
    return true;
};

template<typename FieldT>
bool pp_same_prod<FieldT>::verify(const std::map<std::string, FieldT> &mid_challenges, const std::map<std::string, FieldT> &challenges, const bool output) const {
    bool result = true;
    result &= this->pp_prod_A_.verify(mid_challenges.at("pp_prod_A_y"), mid_challenges.at("pp_prod_A_x0"), challenges.at("pp_prod_A_x"));
    result &= this->pp_prod_B_.verify(mid_challenges.at("pp_prod_B_y"), mid_challenges.at("pp_prod_B_x0"), challenges.at("pp_prod_B_x"));
    result &= this->pp_shift_AB_.verify(mid_challenges.at("pp_shift_AB_y"), mid_challenges.at("pp_shift_AB_x0"), challenges.at("pp_shift_AB_x"));
    
    if ( result ) {
        if (output) printf("pp_same_prod<FieldT>::verify() \033[032mpass\033[037m\n\n");
        return true;
    }
    if (output) printf("pp_same_prod<FieldT>::verify() \033[031mfail\033[037m\n\n");
    return false;
};

template<typename FieldT>
bool same_prod_test() {
    size_t mu, n, col_num, row_num;
    mu = std::rand() % 4 + 3;
    n = std::rand() % 4 + 3;
    // mu = 2;
    // n = 2;
    col_num = std::pow(2UL, mu) * n;
    row_num = std::pow(2UL, mu) * n;

    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::random(row_num, col_num);
    row_vector_matrix<FieldT> B = A.shuffle();

    pp_same_prod<FieldT> same_prod(mu, n, col_num, A, B);
    same_prod.is_satisfy();

    /*
        prove part
    */
    std::map<std::string, FieldT> mid_challenges;
    mid_challenges["pp_prod_A_y"]   = FieldT::random_element();
    mid_challenges["pp_prod_A_x0"]  = FieldT::random_element();
    mid_challenges["pp_prod_B_y"]   = FieldT::random_element();
    mid_challenges["pp_prod_B_x0"]  = FieldT::random_element();
    mid_challenges["pp_shift_AB_y"] = FieldT::random_element();
    mid_challenges["pp_shift_AB_x0"] = FieldT::random_element();

    same_prod.prove(mid_challenges);

    /*
        verify part
    */
    std::map<std::string, FieldT> challenges;
    challenges["pp_prod_A_x"]   = FieldT::random_element();
    challenges["pp_prod_B_x"]   = FieldT::random_element();
    challenges["pp_shift_AB_x"] = FieldT::random_element();

    same_prod.verify(mid_challenges, challenges, true);
    return true;
}
#endif