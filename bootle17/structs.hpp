#ifndef BOOTLE_STRUCTS_HPP
#define BOOTLE_STRUCTS_HPP
#include <iostream>
#include <vector>
#include <libff/common/utils.hpp>

template<typename FieldT>
class row_vector_matrix; // 前置声明

typedef long integer_coeff_t;

template<typename FieldT>
class row_vector {
private:
    std::vector<FieldT> contents;
public:
    row_vector() {}
    row_vector(const size_t size);
    row_vector(const size_t size, const FieldT& value);
    row_vector(const std::vector<FieldT> &contents);
    row_vector(const row_vector<FieldT>& other);
    void add_item(const FieldT& item);
    void set_item(const size_t idx, const FieldT& item);

    size_t size() const;
    const FieldT& get_item(const size_t idx) const;
    const std::vector<FieldT>& get_all_items() const;

    row_vector<FieldT> operator+(const row_vector<FieldT> &other) const;
    row_vector<FieldT> operator-(const row_vector<FieldT> &other) const;
    row_vector<FieldT> operator*(const row_vector<FieldT> &other) const;
    row_vector<FieldT> operator*(const FieldT &field_coeff) const;
    row_vector<FieldT> operator*(const integer_coeff_t integer_coeff) const;
    void operator+=(const row_vector<FieldT> &other);
    void operator*=(const FieldT &field_coeff);
    void operator*=(const row_vector<FieldT> &other);
    bool operator==(const row_vector<FieldT> &other) const;
    bool operator!=(const row_vector<FieldT> &other) const;

    row_vector<FieldT> inverse() const;
    FieldT dot_product(const row_vector<FieldT> &other) const;
    // std::ostream& operator<<(std::ostream &out) const;

    static row_vector<FieldT> random(const size_t len);     // 生成随机行向量
    static row_vector<FieldT> all_one(const size_t len);    // 全 1 的行向量
    static row_vector<FieldT> all_zero(const size_t len);   // 全 0 的行向量
    static row_vector<FieldT> linear_grow(const size_t len, const FieldT& start, const FieldT& step);   // 线性增长
    static FieldT dot_product(const row_vector<FieldT> &row_vec1, const row_vector<FieldT> &row_vec2);      // 内积运算
    
    void set_zero();    // 将向量设置为全 0
    FieldT open(const std::vector<FieldT>& linear_combination) const;

    friend void row_vector_matrix<FieldT>::set_item(const size_t row_num, const size_t col_num, const FieldT& val);
};

template<typename FieldT>
row_vector<FieldT> operator*(const FieldT &field_coeff, const row_vector<FieldT> &row_vec);

template<typename FieldT>
class row_vector_matrix {
private:
    size_t num_column;
    std::vector<row_vector<FieldT> > matrix;
public:
    row_vector_matrix() {}
    row_vector_matrix(size_t num_column): num_column(num_column) {};
    row_vector_matrix(std::vector<row_vector<FieldT> >& matrix, size_t num_column);
    row_vector_matrix(const row_vector_matrix& other);

    void set_num_column(const size_t num_column);
    void add_row_vector(const row_vector<FieldT>& row_vec);

    // row_idx 和 col_idx 从 0 开始, 可选范围为分别为 [0, row_num-1], [0, col_num-1]
    void set_item(const size_t row_idx, const size_t col_idx, const FieldT& val);
    

    size_t get_row_num() const;
    size_t get_column_num() const;

    const FieldT& get_item(const size_t row_idx, const size_t col_idx) const;
    const row_vector<FieldT>& get_row(const size_t idx) const;
    row_vector<FieldT> open(const std::vector<FieldT>& linear_combination) const;

    row_vector_matrix<FieldT> operator+(const row_vector_matrix<FieldT> &other) const;
    // row_vector_matrix<FieldT> operator-(const row_vector_matrix<FieldT> &other) const;
    row_vector_matrix<FieldT> operator*(const row_vector_matrix<FieldT> &other) const;
    // row_vector_matrix<FieldT> operator*(const row_vector_matrix &field_coeff) const;
    // row_vector_matrix<FieldT> operator*(const row_vector_matrix integer_coeff) const;
    // void operator+=(const row_vector_matrix<FieldT> &other);
    bool operator==(const row_vector_matrix<FieldT> &other) const;
    bool operator!=(const row_vector_matrix<FieldT> &other) const;

    /* Consider the matrices A and B, which have mn rows, given respectively by
    vectors ai,j , bi,j ∈ F^k, with 0 ≤ i ≤ m−1, 1 ≤ j ≤ n. The top-left element of A
    is a 1. Columns 2 up to k of A are equal to columns 1 up to k − 1 of B. Further,
    we can obtain the first column of A from the last column of B by deleting the
    last entry c and appending a 1. In this case, A is said to be the shift of B.*/
    row_vector_matrix<FieldT> shift() const;

    // 随机生成一个 row_vector_matrix
    static row_vector_matrix<FieldT> random(const size_t row_num, const size_t col_num);
    // 全 1 的 row_vector_matrix
    static row_vector_matrix<FieldT> all_one(const size_t row_num, const size_t col_num);
    // 全 0 的 row_vector_matrix
    static row_vector_matrix<FieldT> all_zero(const size_t row_num, const size_t col_num);
    /*
        从 1 开始顺序增长的 row_vector_matrix
        1   2   3   ... k-1     k
        k+1 k+2 k+3 ... 2k-1    2k
        ...
    */
    static row_vector_matrix<FieldT> linear_grow(const size_t row_num, const size_t col_num, const FieldT& start, const FieldT& step);
};

template<typename FieldT>
FieldT next_challenge(FieldT &challenge) {
    return challenge * 17 - FieldT::one() * 16;
}

// 输入 x0, 返回 x0, 2x0, 3x0, ..., mu_x0
template<typename FieldT>
std::vector<FieldT> fake_xs(const FieldT& x0, const size_t len);

// 输入 x, n, 返回 x^0, x^1, x^2, ..., x^n
template<typename FieldT>
std::vector<FieldT> get_exps(const FieldT& x, const size_t n);

// 输入 x0, x1, ..., x_(mu-1), 返回 i 取值范围为 [0, 2^mu - 1] 时所有的 (x0^i0)(x1^i1)...(x_(mu-1)^i_(mu-1))
template<typename FieldT>
std::vector<FieldT> get_compress_xs_related(const std::vector<FieldT>& compress_xs);

template<typename FieldT>
std::string fr_to_string(const FieldT& fr_item);

template<typename FieldT>
std::string frs_to_string(const std::vector<FieldT>& fr_vector);

template<typename FieldT>
std::string sha256_to_hex_string(const std::string& pre_image);

libff::bit_vector hexToBin(const std::string& str);

template<typename FieldT>
FieldT sha256_hex_digest_to_Fr(const std::string hex_digest);



void size_t_test();

template <typename FieldT>
void group_test();

template <typename FieldT>
void field_test();

template <typename FieldT>
void matrix_test();
#endif