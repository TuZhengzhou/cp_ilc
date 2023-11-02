// #define CURVE_ALT_BN128
#include <bits/stdc++.h>
#include "structs.tcc"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

int main() {
    typedef libff::Fr<libff::mnt6_pp> FieldT;
    libff::mnt6_pp::init_public_params();

    size_t row, col, start, step;
    row = 3; col = 3; start = 4; step = 1;
    row_vector_matrix<FieldT> tmp = row_vector_matrix<FieldT>::linear_grow(row, col, start, step);
    row_vector_matrix<FieldT> shifted = tmp.shift();
    for (size_t i = 0; i < row; i++) {
        for(size_t j = 0; j < col; j++) {
            std::cout << shifted.get_row(i).get_item(j) << " ";
        }
        std::cout << std::endl;
    }

    return 1;
}