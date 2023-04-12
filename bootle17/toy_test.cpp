#include <iostream>
#include <unordered_map>
#include <string> 
#include <vector>
#include "structs.tcc"
#include "permutation.tcc"

#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
using namespace std;

string linearEncode(string message) {
    unordered_map<char, int> charCount;
    for (char c : message) {
        charCount[c]++;
    }

    string code;
    for (char c : message) {
        if (charCount[c] == 0) continue; // already encoded
        code += to_string(charCount[c]); // add count to code
        code += c; // add character to code
        charCount[c] = 0; // mark as encoded
    }

    return code;
}


using namespace std;

class Cycle {
public:
    vector<int> elements;

    Cycle() {}

    Cycle(const vector<int>& v) : elements(v) {}

    bool contains(int x) const {
        for (int i = 0; i < elements.size(); i++) {
            if (elements[i] == x) {
                return true;
            }
        }
        return false;
    }

    Cycle apply(const Cycle& c) const {
        vector<int> result(elements.size());
        for (int i = 0; i < elements.size(); i++) {
            int index = c.indexOf(elements[i]);
            if (index != -1) {
                result[i] = c.elements[(index + 1) % c.elements.size()];
            } else {
                result[i] = elements[i];
            }
        }
        return Cycle(result);
    }

    int indexOf(int x) const {
        for (int i = 0; i < elements.size(); i++) {
            if (elements[i] == x) {
                return i;
            }
        }
        return -1;
    }

    friend ostream& operator<<(ostream& os, const Cycle& c) {
        os << "(";
        for (int i = 0; i < c.elements.size(); i++) {
            os << c.elements[i];
            if (i < c.elements.size() - 1) {
                os << " ";
            }
        }
        os << ")";
        return os;
    }
};


int main() {
    typedef libff::Fr<libff::mnt6_pp> FieldT;
    libff::mnt6_pp::init_public_params();

    string message = "abbcccddddeeeee";
    string code = linearEncode(message);
    cout << "Original message: " << message << endl;
    cout << "Encoded message: " << code << endl;

    Cycle c1({1, 2, 3});
    Cycle c2({2, 3, 1});
    Cycle c3 = c1.apply(c2);
    cout << c1 << " apply " << c2 << " = " << c3 << endl;

    map<tuple<int, int>, string> myMap_tuple;
    myMap_tuple[make_tuple(1, 2)] = "Hello";
    myMap_tuple[make_tuple(3, 4)] = "World";
    cout << myMap_tuple[make_tuple(1, 2)] << endl;
    cout << myMap_tuple[make_tuple(3, 4)] << endl;

    map<pair<int, int>, string> myMap_pair;
    myMap_pair[make_pair(1, 2)] = "Hello";
    myMap_pair[make_pair(3, 4)] = "World";
    cout << myMap_pair[make_pair(1, 2)] << endl;
    cout << myMap_pair[make_pair(3, 4)] << endl;


    std::vector<tuple<int, int> > vec = {make_tuple(1,2), make_tuple(3,4)};
    cycle<tuple<int, int>, 2> myCycle = cycle<tuple<int, int>, 2>(vec);

    std::cout << myCycle << std::endl;

    // std::vector<std::vector<tuple<size_t, size_t> > > permutation_vec = {{make_tuple(0,0), make_tuple(1,1)}, {make_tuple(0,1), make_tuple(1,0)}};
    // permutation<tuple<size_t, size_t>, 2 > myPermtation = permutation<tuple<size_t, size_t>, 2 >(permutation_vec);

    // std::cout << myPermtation.get_cycle(0) << std::endl;
    // std::cout << myPermtation.get_cycle(1) << std::endl;

    // permutation<tuple<size_t, size_t>, 2 > myPermtation1 = permutation<tuple<size_t, size_t>, 2 >(myPermtation);
    // std::cout << myPermtation1.get_cycle(0) << std::endl;
    // std::cout << myPermtation1.get_cycle(1) << std::endl;

    // row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::linear_grow(2, 2, 1, 1);
    // row_vector_matrix<FieldT> B = apply_permutation(A, myPermtation);
    // B.print();

    typedef tuple<size_t, size_t, size_t> dim3;
    size_t cycle_len = 3;
    size_t cycle_num = 3;
    
    // const size_t dim_num = 3;
    // std::vector<size_t> dim_limits = {6,7,8};
    // permutation<dim3, dim_num> rand_perm = permutation<dim3, dim_num>::random_permutation(cycle_len, cycle_num, dim_num, dim_limits);
    // for (size_t i = 0; i < cycle_num; i++) {
    //     std::cout << rand_perm.get_cycle(i) << std::endl;
    // }

    typedef tuple<size_t, size_t> dim2;
    cycle_len = 2;
    cycle_num = 1;
    const size_t dim_num = 2;
    std::vector<size_t> dim_limits = {6,7};
    permutation<dim2, dim_num> rand_perm = permutation<dim2, dim_num>::random_permutation(cycle_len, cycle_num, dim_num, dim_limits);
    for (size_t i = 0; i < cycle_num; i++) {
        std::cout << rand_perm.get_cycle(i) << std::endl;
    }

    row_vector_matrix<FieldT> A = row_vector_matrix<FieldT>::linear_grow(dim_limits[0], dim_limits[1], 1, 1);
    row_vector_matrix<FieldT> B = apply_permutation(A, rand_perm);
    B.print();

    pp_perm<FieldT> perm(7, 1, 3, A, A, rand_perm);

    return 0;
}
