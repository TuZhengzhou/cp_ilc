#include <bits/stdc++.h>
#include "structs.tcc"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

int main() {
    typedef libff::Fr<libff::mnt6_pp> FieldT;
    libff::mnt6_pp::init_public_params();

    size_t row = 3;
    size_t col = 3;
    row_vector_matrix<FieldT> tmp = row_vector_matrix<FieldT>::all_one(row, col) * 2;

    tmp.print();

    tmp.partial_products(false).print();

    tmp.partial_products(true).print();

    std::vector<FieldT> flattened = tmp.partial_products(false).flatten();
    std::cout << flattened;

    row_vector_matrix<FieldT> from_vec(flattened, tmp.get_row_num(), tmp.get_column_num());
    from_vec.print();

    from_vec.shuffle().print();
}