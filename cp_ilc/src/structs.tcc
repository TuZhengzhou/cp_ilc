#ifndef BOOTLE_STRUCTS_TCC
#define BOOTLE_STRUCTS_TCC
#include <vector>
#include <map>
#include <algorithm>
#include <openssl/sha.h>
#include "structs.hpp"

namespace py_pr{
	template<typename T>
	inline std::ostream& out_put(std::ostream& o,const T & x){
		return o<<x;
	}
	inline std::ostream& out_put(std::ostream& o,const std::string& x){
		return o<<"\""<<x<<"\"";
	}
	inline std::ostream& out_put(std::ostream& o,const char* & x){
		return o<<"\""<<x<<"\"";
	}
	inline std::ostream& out_put(std::ostream& o,const char & x){
		return o<<"\""<<x<<"\"";
	}
	template<typename T1,typename T2>
	inline std::ostream& out_put(std::ostream& o,const std::pair<T1,T2> & x){
		out_put(o,x.first);
		o<<": ";
		out_put(o,x.second);
		return o;
	}			
}

template<typename T>
std::ostream& operator<<(std::ostream &o,std::vector<T> &x){
    o<<"[";
    for(auto i=x.begin();i<x.end();++i){
		//可以直接for(auto i:x)，但是我不知道怎么特判第一个来控制","
        if(i!=x.begin()) o<<", ";
       	// py_pr::out_put(o,*i);
        o << *i;
    }
    o<<"]"<<std::endl;
    return o;
}

template<typename T1,typename T2>
std::ostream& operator<<(std::ostream &o,std::map<T1,T2> &x){
    o<<"{";
	//类似python的格式
    for(auto i=x.begin();i!=x.end();++i){
        if(i!=x.begin()) o<<", ";
        py_pr::out_put(o,*i);
        // o << *i;
    }
    o<<"}"<<std::endl;
    return o;
}

template<typename FieldT>
row_vector<FieldT>::row_vector(const size_t size, const FieldT& value) {
    this->contents = std::vector<FieldT>(size, value);
}

template<typename FieldT>
row_vector<FieldT>::row_vector(const size_t size) {
    this->contents = std::vector<FieldT>(size, FieldT::zero());
}

template<typename FieldT>
row_vector<FieldT>::row_vector(const std::vector<FieldT> &contents) :
    contents(contents)
{
}

template<typename FieldT>
row_vector<FieldT>::row_vector(const row_vector<FieldT>& other)
{
    this->contents = other.get_all_items();
}

template<typename FieldT>
size_t row_vector<FieldT>::size() const {
    return this->contents.size();
}

template<typename FieldT>
void row_vector<FieldT>::add_item(const FieldT& item){
    this->contents.emplace_back(item);
}

template<typename FieldT>
void row_vector<FieldT>::set_item(const size_t idx, const FieldT& item) {
    assert( idx < this->contents.size() );
    this->contents[idx] = item;
}

template<typename FieldT>
const FieldT& row_vector<FieldT>::get_item(const size_t idx) const {
    assert(idx <= this->contents.size());
    return this->contents[idx];
}

template<typename FieldT>
const std::vector<FieldT>& row_vector<FieldT>::get_all_items() const {
    // return std::vector<FieldT>(this->contents);
    return this->contents;
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator+(const row_vector<FieldT> &other) const {
    assert(this->size() == other.size());
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(this->contents[i] + other.contents[i]);
    }
    return row_vector(contents);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator-(const row_vector<FieldT> &other) const {
    assert(this->size() == other.size());
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(this->contents[i] - other.contents[i]);
    }
    return row_vector(contents);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator-() const {
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(- this->contents[i]);
    }
    return row_vector(contents);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator*(const row_vector<FieldT> &other) const {
    assert(this->size() == other.size());
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(this->contents[i] * other.contents[i]);
    }
    return row_vector(contents);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator*(const integer_coeff_t integer_coeff) const {
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(this->contents[i] * integer_coeff);
    }
    return row_vector(contents);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::operator*(const FieldT &field_coeff) const {
    std::vector<FieldT> contents;
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back(this->contents[i] * field_coeff);
    }
    return row_vector(contents);
}

template<typename FieldT>
void row_vector<FieldT>::operator+=(const row_vector<FieldT> &other) {
    assert(this->size() == other.size());
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        this->contents[i] += other.contents[i];
    }
}

template<typename FieldT>
void row_vector<FieldT>::operator*=(const FieldT &field_coeff) {
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        this->contents[i] *= field_coeff;
    }
}

template<typename FieldT>
void row_vector<FieldT>::operator*=(const row_vector<FieldT> &other) {
    assert(this->size() == other.size());
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        this->contents[i] *= other.contents[i];
    }
}

template<typename FieldT>
bool row_vector<FieldT>::operator==(const row_vector<FieldT> &other) const {
    assert(this->size() == other.size());
    return (this->contents == other.contents);
}

template<typename FieldT>
bool row_vector<FieldT>::operator!=(const row_vector<FieldT> &other) const {
    return !((*this) == other);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::inverse() const {
    std::vector<FieldT> result_t;
    for( auto item: this->contents ) {
        assert(item != FieldT::zero());
        result_t.emplace_back(item.inverse());
    }
    return row_vector<FieldT>(result_t);
}

// template<typename FieldT>
// std::ostream& operator<<(std::ostream &out, const std::vector<FieldT> &vec) {
//     for (auto item : vec) {
//         out << item << std::endl;
//     }
//     return out;
// }

// template<typename FieldT>
// std::ostream& row_vector<FieldT>::operator<<(std::ostream &out) const {
//     out << "size = %ld" << this->size() << std::endl;
//     out << this->contents << std::endl;
//     return out;
// }
template<typename FieldT>
FieldT row_vector<FieldT>::dot_product(const row_vector<FieldT> &other) const {
    assert(this->size() == other.size());
    FieldT result = FieldT::zero();
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        result += this->get_item(i) * other.get_item(i);
    }
    return result;
}

template<typename FieldT>
FieldT row_vector<FieldT>::dot_product(const row_vector<FieldT> &row_vec1, const row_vector<FieldT> &row_vec2) {
    return row_vec1.dot_product(row_vec2);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::random(const size_t len) {
    std::vector<FieldT> tmp;
    for (size_t i = 0; i < len; i++) {
        tmp.emplace_back(FieldT::random_element());
    }
    return row_vector<FieldT>(tmp);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::all_one(const size_t len) {
    std::vector<FieldT> tmp;
    for (size_t i = 0; i < len; i++) {
        tmp.emplace_back(FieldT::one());
    }
    return row_vector<FieldT>(tmp);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::all_zero(const size_t len) {
    std::vector<FieldT> tmp;
    for (size_t i = 0; i < len; i++) {
        tmp.emplace_back(FieldT::zero());
    }
    return row_vector<FieldT>(tmp);
}

template<typename FieldT>
row_vector<FieldT> row_vector<FieldT>::linear_grow(const size_t len, const FieldT& start, const FieldT& step) {
    std::vector<FieldT> tmp;
    for (size_t i = 0; i < len; i++) {
        tmp.emplace_back(start + step * i);
        // std::cout << tmp.back() << " ";
    }
    // std::cout << std::endl;
    return row_vector<FieldT>(tmp);
}

template<typename FieldT>
void row_vector<FieldT>::set_zero() {
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        this->contents[i] = FieldT::zero();
    }
}

template<typename FieldT>
FieldT row_vector<FieldT>::open(const std::vector<FieldT>& linear_combination) const {
    assert(this->size() == linear_combination.size());
    FieldT result = FieldT::zero();
    size_t upper_bound = this->size();
    for (size_t i = 0; i < upper_bound; i++) {
        result += this->contents[i] * linear_combination[i];
    }
    return result;
}

template<typename FieldT>
row_vector<FieldT> operator*(const FieldT &field_coeff, const row_vector<FieldT> &row_vec) {
    std::vector<FieldT> contents;
    size_t upper_bound = row_vec.size();
    for (size_t i = 0; i < upper_bound; i++) {
        contents.emplace_back((row_vec.get_item(i)) * field_coeff);
    }
    return row_vector<FieldT>(contents);
}

template<typename FieldT>
row_vector_matrix<FieldT>::row_vector_matrix(const std::vector<FieldT>& vec, const size_t num_row, const size_t num_column) {
    assert(vec.size() == num_row * num_column);
    this->num_column = num_column;
    this->matrix.clear();

    for (size_t i = 0; i < num_row; i++) {
        this->add_row_vector(row_vector<FieldT>(std::vector<FieldT>(vec.begin() + i * num_column, vec.begin() + (i+1) * num_column)));
    }
}


template<typename FieldT>
row_vector_matrix<FieldT>::row_vector_matrix(std::vector<row_vector<FieldT> >& matrix, size_t num_column) {
    size_t upper_bound = matrix.size();
    for (size_t i = 0; i < upper_bound; i++) {
        assert(matrix[i].size() == num_column);
    }
    this->num_column = num_column;
    this->matrix = matrix;
};

template<typename FieldT>
row_vector_matrix<FieldT>::row_vector_matrix(const row_vector_matrix& other) {
    this->num_column = other.get_column_num();

    this->matrix.clear();
    size_t upper_bound = other.get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        this->matrix.emplace_back(row_vector<FieldT>(other.get_row(i)));
    }
};

template<typename FieldT>
void row_vector_matrix<FieldT>::set_num_column(const size_t num_column) {
    assert(this->num_column == 0);
    this->num_column = num_column;
}

template<typename FieldT>
void row_vector_matrix<FieldT>::add_row_vector(const row_vector<FieldT>& row_vec) {
    assert(this->num_column == row_vec.size());
    this->matrix.emplace_back(row_vector<FieldT>(row_vec));
}

template<typename FieldT>
void row_vector_matrix<FieldT>::set_item(const size_t row_idx, const size_t col_idx, const FieldT& val) {
    assert(row_idx < this->get_row_num() && col_idx < this->get_column_num());
    this->matrix[row_idx].contents[col_idx] = val;
}

template<typename FieldT>
void row_vector_matrix<FieldT>::set_item(const size_t idx, const FieldT& val) {
    size_t row_idx = idx / this->get_column_num();
    size_t col_idx = idx % this->get_column_num();
    assert(row_idx < this->get_row_num() && col_idx < this->get_column_num());
    this->matrix[row_idx].contents[col_idx] = val;
}

template<typename FieldT>
size_t row_vector_matrix<FieldT>::get_row_num() const {
    return this->matrix.size();
}

template<typename FieldT>
size_t row_vector_matrix<FieldT>::get_column_num() const {
    return this->num_column;
}

template<typename FieldT>
const FieldT& row_vector_matrix<FieldT>::get_item(const size_t row_idx, const size_t col_idx) const {
    assert(row_idx < this->get_row_num());
    return this->matrix[row_idx].get_item(col_idx);
}

template<typename FieldT>
const FieldT& row_vector_matrix<FieldT>::get_item(const size_t idx) const {
    size_t row_idx = idx / this->get_column_num();
    size_t col_idx = idx % this->get_column_num();
    assert(row_idx < this->get_row_num());
    return this->matrix[row_idx].get_item(col_idx);
}

template<typename FieldT>
const row_vector<FieldT>& row_vector_matrix<FieldT>::get_row(const size_t idx) const {
    assert(idx < this->get_row_num());
    return this->matrix[idx];
}

template<typename FieldT>
row_vector<FieldT> row_vector_matrix<FieldT>::open(const std::vector<FieldT>& linear_combination) const {
    assert(this->get_row_num() == linear_combination.size());

    row_vector<FieldT> result = row_vector<FieldT>::all_zero(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        result += this->matrix[i] * linear_combination[i];
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::operator+(const row_vector_matrix<FieldT> &other) const {
    assert(this->num_column == other.get_column_num());
    assert(this->get_row_num() == other.get_row_num());

    row_vector_matrix<FieldT> result(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        result.add_row_vector(this->get_row(i) + other.get_row(i));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::operator-(const row_vector_matrix<FieldT> &other) const {
    assert(this->num_column == other.get_column_num());
    assert(this->get_row_num() == other.get_row_num());

    row_vector_matrix<FieldT> result(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        result.add_row_vector(this->get_row(i) - other.get_row(i));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::operator*(const row_vector_matrix<FieldT> &other) const {
    assert(this->num_column == other.get_column_num());
    assert(this->get_row_num() == other.get_row_num());

    row_vector_matrix<FieldT> result(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        result.add_row_vector(this->get_row(i) * other.get_row(i));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::operator*(const FieldT &field_coeff) const {

    row_vector_matrix<FieldT> result(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        result.add_row_vector(this->get_row(i) * field_coeff);
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::operator*(const integer_coeff_t integer_coeff) const {
    return (*this) * (FieldT::one() * integer_coeff);
}

template<typename FieldT>
bool row_vector_matrix<FieldT>::operator==(const row_vector_matrix<FieldT> &other) const {
    assert(this->num_column == other.get_column_num());
    assert(this->get_row_num() == other.get_row_num());

    row_vector_matrix<FieldT> result(this->num_column);
    size_t upper_bound = this->get_row_num();
    for (size_t i = 0; i < upper_bound; i++) {
        if (this->get_row(i) != other.get_row(i)) {
            return false;
        }
    }
    return true;
}

template<typename FieldT>
bool row_vector_matrix<FieldT>::operator!=(const row_vector_matrix<FieldT> &other) const {
    return !((*this) == other);
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::random(const size_t row_num, const size_t col_num) {
    row_vector_matrix<FieldT> result(col_num);
    for (size_t i = 0; i < row_num; i++) {
        result.add_row_vector(row_vector<FieldT>::random(col_num));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::shift() const {
    
    size_t row_num = this->get_row_num();
    size_t col_num = this->get_column_num();
    row_vector_matrix<FieldT> result = row_vector_matrix<FieldT>::all_one(row_num, col_num);
    
    for (size_t row = 0; row < row_num; row++) {
        for (size_t col = 0; col < col_num - 1; col++) {
            result.set_item(row, col+1, this->get_item(row, col));
        }
        if (row + 1 < row_num) {
            result.set_item(row+1, 0, this->get_item(row, col_num-1));
        }
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::partial_products(const bool choice) const {
    size_t row_num = this->get_row_num();
    size_t col_num = this->get_column_num();
    row_vector_matrix<FieldT> result = row_vector_matrix<FieldT>::all_one(row_num, col_num);

    if (choice == false) {
        for (size_t row = 0; row < row_num; row++) {
            for (size_t col = 0; col < col_num - 1; col++) {
                result.set_item(row, col+1, result.get_item(row, col) * this->get_item(row, col));
            }
            if (row + 1 < row_num) {
                result.set_item(row+1, 0, result.get_item(row, col_num-1)* this->get_item(row, col_num-1));
            }
        }
    } else {
        FieldT pre_partial_sum = FieldT::one();
        for (size_t row = 0; row < row_num; row++) {
            for (size_t col = 0; col < col_num; col++) {
                result.set_item(row, col, pre_partial_sum * this->get_item(row, col));
                pre_partial_sum = result.get_item(row, col);
            }
        }
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> row_vector_matrix<FieldT>::flatten() const {
    std::vector<FieldT> result;
    size_t row_num = this->get_row_num();
    for (size_t i = 0; i < row_num; i++) {
        result.insert(result.end(), this->get_row(i).get_all_items().begin(), this->get_row(i).get_all_items().end());
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::shuffle() const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<FieldT> flattened = this->flatten();
    std::shuffle(flattened.begin(), flattened.end(), gen);

    return row_vector_matrix<FieldT>(flattened, this->get_row_num(), this->get_column_num());
}

template <typename FieldT, size_t N>
row_vector_matrix<FieldT> apply_permutation(const row_vector_matrix<FieldT>& matrix, const permutation<std::tuple<size_t, size_t>, N>& perm) {
    
    row_vector_matrix<FieldT> result = row_vector_matrix<FieldT>(matrix);
    
    size_t cycle_num = perm.num_cycles();
    for (size_t i = 0; i < cycle_num; i++) {
        const std::vector<std::tuple<size_t, size_t>>& cycle_tmp = perm.get_cycle(i).get_content_vector();
        size_t cycle_len = cycle_tmp.size();

        size_t row_1, col_1, row_2, col_2;
        std::tie(row_1, col_1) = cycle_tmp[0];
        FieldT head_val = result.get_item(row_1, col_1);

        for (size_t j = 1; j < cycle_len; j++) {
            std::tie(row_2, col_2) = cycle_tmp[j];
            result.set_item(row_1, col_1, result.get_item(row_2, col_2));
            row_1 = row_2;
            col_1 = col_2;
        }
        result.set_item(row_1, col_1, head_val);
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::all_one(const size_t row_num, const size_t col_num) {
    row_vector_matrix<FieldT> result(col_num);
    for (size_t i = 0; i < row_num; i++) {
        result.add_row_vector(row_vector<FieldT>::all_one(col_num));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::all_zero(const size_t row_num, const size_t col_num) {
    row_vector_matrix<FieldT> result(col_num);
    for (size_t i = 0; i < row_num; i++) {
        result.add_row_vector(row_vector<FieldT>::all_zero(col_num));
    }
    return result;
}

template<typename FieldT>
row_vector_matrix<FieldT> row_vector_matrix<FieldT>::linear_grow(const size_t row_num, const size_t col_num, const FieldT& start, const FieldT& step) {
    row_vector_matrix<FieldT> result(col_num);
    for (size_t i = 0; i < row_num; i++) {
        result.add_row_vector(row_vector<FieldT>::linear_grow(col_num, start + step * col_num * i, step));
    }
    return result;
}

template<typename FieldT>
void row_vector_matrix<FieldT>::print() const {
    size_t row_num = this->get_row_num();
    size_t col_num = this->get_column_num();
    for (size_t i = 0; i < row_num; i++) {
        for(size_t j = 0; j < col_num; j++) {
            std::cout << this->get_row(i).get_item(j) << " ";
        }
        std::cout << std::endl;
    }
}

template<typename FieldT>
std::vector<FieldT> fake_xs(const FieldT& x0, const size_t len) {
    assert(x0 != FieldT::zero());
    std::vector<FieldT> result = {x0};
    for (size_t i = 2; i <= len; i++) {
        result.emplace_back(x0 * i);
    }
    return result;
}

template<typename FieldT>
std::vector<FieldT> get_exps(const FieldT& x, const size_t n) {
    std::vector<FieldT> result(n + 1);
    result[0] = FieldT::one();

    size_t upper_bound = n + 1;
    for (size_t i = 1; i < upper_bound; i++) {
        result[i] = result[i-1] * x;
    }
    return result;
}

template<typename FieldT>
std::vector<FieldT> get_compress_xs_related(const std::vector<FieldT>& compress_xs) {
    size_t mu_ = compress_xs.size();
    size_t m_ = std::pow(2ul, mu_);
    std::vector<FieldT> result(m_, FieldT::one());
    for (size_t i = 0; i < mu_; i++) {
        size_t step = std::pow(2, i+1);
        size_t duration = step / 2;
        for (size_t j = duration; j < m_; j += step) {
            for (size_t k = 0; k < duration; k++) {
                result[j+k] *= compress_xs[i];
            }
        }
    }
    return result;
}

template<typename FieldT>
std::string fr_to_string(const FieldT& fr_item) {
    std::stringstream ss;
    ss << fr_item.as_bigint();
    return ss.str();
}

template<typename FieldT>
std::string frs_to_string(const std::vector<FieldT>& fr_vector) {
    std::string result;
    for (auto fr: fr_vector) {
        result += fr_to_string<FieldT>(fr);
    }
    return result;
}

std::string hexToChar(const char c) {
    switch(tolower(c))
    {
        case '0': return "0000";
        case '1': return "0001";
        case '2': return "0010";
        case '3': return "0011";
        case '4': return "0100";
        case '5': return "0101";
        case '6': return "0110";
        case '7': return "0111";
        case '8': return "1000";
        case '9': return "1001";
        case 'a': return "1010";
        case 'b': return "1011";
        case 'c': return "1100";
        case 'd': return "1101";
        case 'e': return "1110";
        case 'f': return "1111";
    }
    return "";
}

libff::bit_vector hexToBin(const std::string& str) {
    libff::bit_vector res;
    for (auto item : str) {
        std::string hexItem = hexToChar(item);
        res.push_back(hexItem[0] == '1' ? true : false);
        res.push_back(hexItem[1] == '1' ? true : false);
        res.push_back(hexItem[2] == '1' ? true : false);
        res.push_back(hexItem[3] == '1' ? true : false);
    }
    return res;
}

std::string sha256_to_hex_string(const std::string& pre_image) {
    // 使用 SHA256 哈希函数计算消息的哈希值
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256((const unsigned char*)pre_image.c_str(), pre_image.length(), hash);

    // 将哈希值打印为十六进制字符串
    std::string hash_hex;
    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        char buf[3];
        sprintf(buf, "%02x", hash[i]);
        hash_hex += buf;
    }
    return hash_hex;
}

template<typename FieldT>
FieldT sha256_hex_digest_to_Fr(const std::string hex_digest) {
    assert(hex_digest.size() == 64UL);
    assert(hex_digest.size() <= FieldT::size_in_bits());    // size_in_bits 在 <libff/algebra/fields/fp.hpp> 
    libff::bit_vector hash_bits = hexToBin(hex_digest);
    return libff::convert_bit_vector_to_field_element<FieldT>(hash_bits);
}


void size_t_test() {
    size_t len = rand() % 10101;
    integer_coeff_t coeff = rand() % 10101;
    // printf("len = %ld, coeff = %ld\n", len, coeff);

    std::vector<size_t> example;
    std::vector<size_t> sum;
    std::vector<size_t> sub;
    std::vector<size_t> mult_coeff;
    std::vector<size_t> hadmard_product;
    for (size_t i = 0; i < len; i++) {
        example.emplace_back(1 * i);        // {0, 1, ..., len-1}
        sum.emplace_back(1 * i + 1 * i);    // {0, 2, ..., 2(len-1)}
        sub.emplace_back(0);
        mult_coeff.emplace_back(1 * i * coeff);     // {0, coeff, ..., coeff(len-1)}
        hadmard_product.emplace_back(1 * i * i);    // {0, 1*1, ..., (len-1)*(len-1)}
    }
    row_vector<size_t> vec1 = row_vector<size_t>(example);
    row_vector<size_t> vec2 = row_vector<size_t>(example);

    assert(vec1 == vec2);
    assert((vec1+vec2) == row_vector<size_t>(sum));
    assert((vec1-vec2) == row_vector<size_t>(sub));
    assert((vec1*coeff) == row_vector<size_t>(mult_coeff));
    assert((vec1*vec2) == row_vector<size_t>(hadmard_product));
    std::cout << "size_t_test \033[32mpass\033[37m" << std::endl;
}

template <typename FieldT>
void group_test() {
    size_t len = rand() % 10101;
    FieldT coeff = rand() % 10101;
    // printf("len = %ld, coeff = %ld\n", len, coeff);
    std::cout << "len = " << len << ", coeff = " << coeff << std::endl;

    std::vector<FieldT> example;
    std::vector<FieldT> sum;
    std::vector<FieldT> sub;
    std::vector<FieldT> mult_coeff;
    for (size_t i = 0; i < len; i++) {
        example.emplace_back(FieldT::one() * i);
        sum.emplace_back(FieldT::one() * i + FieldT::one() * i);
        sub.emplace_back(FieldT::zero());
        mult_coeff.emplace_back(FieldT::one() * i * coeff);
    }

    row_vector<FieldT> vec1 = row_vector<FieldT>(example);
    row_vector<FieldT> vec2 = row_vector<FieldT>(example);
    
    assert(vec1 == vec2);
    assert((vec1+vec2) == row_vector<FieldT>(sum));
    assert((vec1-vec2) == row_vector<FieldT>(sub));
    assert((vec1*coeff) == row_vector<FieldT>(mult_coeff));
    assert((coeff*vec1) == row_vector<FieldT>(mult_coeff));
    std::cout << "group_test \033[32mpass\033[37m" << std::endl;
}

template <typename FieldT>
void field_test() {
    size_t len = rand() % 1011;
    FieldT coeff = rand() % 1011;
    // printf("len = %ld, coeff = %ld\n", len, coeff);
    std::cout << "len = " << len << ", coeff = " << coeff << std::endl;

    std::vector<FieldT> example;
    std::vector<FieldT> sum;
    std::vector<FieldT> sub;
    std::vector<FieldT> mult_coeff;
    std::vector<FieldT> hadmard_product;
    for (size_t i = 0; i < len; i++) {
        example.emplace_back(FieldT::one() * i);
        sum.emplace_back(FieldT::one() * i + FieldT::one() * i);
        sub.emplace_back(FieldT::zero());
        mult_coeff.emplace_back(FieldT::one() * i * coeff);
        hadmard_product.emplace_back(FieldT::one() * i * FieldT::one()* i);
    }

    assert(FieldT::one() == size_t(1));

    row_vector<FieldT> vec1 = row_vector<FieldT>(example);
    row_vector<FieldT> vec2 = row_vector<FieldT>(example);
    
    assert(vec1 == vec2);
    assert((vec1+vec2) == row_vector<FieldT>(sum));
    assert((vec1-vec2) == row_vector<FieldT>(sub));
    assert((vec1*coeff) == row_vector<FieldT>(mult_coeff));
    assert((coeff*vec1) == row_vector<FieldT>(mult_coeff));
    assert((vec1*vec2) == row_vector<FieldT>(hadmard_product));
    std::cout << "field_test \033[32mpass\033[37m" << std::endl;
}

template <typename FieldT>
void matrix_test() {
    size_t row_num = std::rand() % 101;
    size_t col_num = std::rand() % 101;

    row_vector_matrix<FieldT> example(col_num);
    row_vector<FieldT> result(std::vector<FieldT>(col_num, FieldT::zero()));
    std::vector<FieldT> linear_combination;

    for (size_t i = 0; i < row_num; i++) {
        row_vector<FieldT> rand_row_vec = row_vector<FieldT>::random(col_num);
        example.add_row_vector(rand_row_vec);
        size_t coeff = std::rand();;
        linear_combination.emplace_back(FieldT::one() * coeff);
        result += rand_row_vec * coeff;
    }

    size_t coeff = std::rand();
    assert(result == example.open(linear_combination));
    std::cout << "matrix_test \033[32mpass\033[37m" << std::endl;
}

template <typename T, size_t N>
cycle<T, N>::cycle(const std::vector<T>& contents) {
    this->cycle_ = std::vector<T>(contents.begin(), contents.end());
}

template <typename T, size_t N>
cycle<T, N>::cycle(const cycle& other) {
    const std::vector<T>& other_contents = other.get_content_vector();
    this->cycle_ = std::vector<T>(other_contents.begin(), other_contents.end());
}

template <typename T, size_t N>
const std::vector<T>& cycle<T, N>::get_content_vector() const {
    return this->cycle_;
}

template <typename T, size_t N>
cycle<T, N> cycle<T, N>::random_cycle(const size_t cycle_len, const size_t dim, const std::vector<size_t>& dim_limits) {

}

template <typename T, size_t N>
permutation<T, N>::permutation(const std::vector<std::vector<T>> content) {
    size_t cycle_num = content.size();
    for (size_t i = 0; i < cycle_num; i++) {
        this->permutation_.push_back(cycle<T, N>(content[i]));
    }
}

template <typename T, size_t N>
permutation<T, N>::permutation(const permutation<T, N>& other) {
    size_t cycle_num = other.num_cycles();
    for (size_t i = 0; i < cycle_num; i++) {
        this->permutation_.push_back(cycle<T, N>(other.get_cycle(i)));
    }
}

template <typename T, size_t N>
const size_t permutation<T, N>::num_cycles() const {
    return this->permutation_.size();
}

template <typename T, size_t N>
const cycle<T, N>& permutation<T, N>::get_cycle(const size_t idx) const {
    assert(idx < this->num_cycles());
    return this->permutation_.at(idx);
}

template<typename T, size_t... Is>
auto vector_to_tuple_helper(const std::vector<T>& v, std::index_sequence<Is...>) {
    return std::make_tuple(v[Is]...);
}

template<typename T, size_t N>
auto vector_to_tuple(const std::vector<T>& v) {
    return vector_to_tuple_helper(v, std::make_index_sequence<N>{});
}

template <typename T, std::size_t N>
permutation<T, N> permutation<T, N>::random_permutation(const size_t cycle_len, const size_t cycle_num, \
                                                    const size_t dim_num, const std::vector<size_t>& dim_limits) {
    assert(dim_num >= 1);
    assert(dim_num == dim_limits.size());

    size_t limit_max = 0;
    size_t total_item_num = 1ul;
    for (size_t limit: dim_limits) { 
        total_item_num *= limit; 
        limit_max = limit > limit_max ? limit : limit_max;
    };
    assert(cycle_len * cycle_num <= total_item_num);
    assert(total_item_num < SIZE_MAX);

    permutation<T, N> result;

    size_t i, j, k;

    /*
        部分乘积, 表示每一维度对应元素数目
    */
    std::vector<size_t> part_prod_inv;
    for (i = 0; i < dim_num; i++) {
        if (i == 0) {
            part_prod_inv.emplace_back(total_item_num / dim_limits[0]);
        } else {
            part_prod_inv.emplace_back(part_prod_inv[i-1] / dim_limits[i]);
        }
    }
    // std::cout << "part_prod_inv: \n" << part_prod_inv << std::endl;

    /* 设置随机数生成器 */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, limit_max);

    size_t mark_idx;
    std::vector<bool> mark_board = std::vector<bool>(total_item_num, false);
    for ( i = 0; i < cycle_num; i++) {

        std::vector<T> cycle_vector;
        for ( j = 0; j < cycle_len; j++) {

            mark_idx = 0;
            std::vector<size_t> dim_idxs(dim_num);
            for ( k = 0; k < dim_num; k++) {
                dim_idxs[k] = dis(gen) % dim_limits[k];
                mark_idx += dim_idxs[k] * part_prod_inv[k];
            }

            if (mark_board[mark_idx] == false) {
                mark_board[mark_idx] = true;
                cycle_vector.push_back(vector_to_tuple<size_t, N>(dim_idxs));
            } else {
                j -= 1;
            }
        }

        result.permutation_.emplace_back(cycle<T, N>(cycle_vector));
    }
    return result;
}

#endif