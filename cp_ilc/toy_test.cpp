#include <iostream>
#include <unordered_map>
#include <string> 
#include <vector>
#include <cassert>
#include <assert.h>
// #include "../libsnark/depends/libff/libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
// #include <libsnark/common/default_types/r1cs_ppzkadsnark_pp.hpp>
#include "../libsnark/libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp"
#include "../libsnark/libsnark/knowledge_commitment/knowledge_commitment.hpp"
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include "structs.tcc"
#include "permutation.tcc"

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


// int main() {
//     typedef libff::Fr<libff::mnt6_pp> FieldT;
//     libff::mnt6_pp::init_public_params();

//     string message = "abbcccddddeeeee";
//     string code = linearEncode(message);
//     cout << "Original message: " << message << endl;
//     cout << "Encoded message: " << code << endl;

//     Cycle c1({1, 2, 3});
//     Cycle c2({2, 3, 1});
//     Cycle c3 = c1.apply(c2);
//     cout << c1 << " apply " << c2 << " = " << c3 << endl;

//     map<tuple<int, int>, string> myMap_tuple;
//     myMap_tuple[make_tuple(1, 2)] = "Hello";
//     myMap_tuple[make_tuple(3, 4)] = "World";
//     cout << myMap_tuple[make_tuple(1, 2)] << endl;
//     cout << myMap_tuple[make_tuple(3, 4)] << endl;

//     map<pair<int, int>, string> myMap_pair;
//     myMap_pair[make_pair(1, 2)] = "Hello";
//     myMap_pair[make_pair(3, 4)] = "World";
//     cout << myMap_pair[make_pair(1, 2)] << endl;
//     cout << myMap_pair[make_pair(3, 4)] << endl;


//     std::vector<tuple<int, int> > vec = {make_tuple(1,2), make_tuple(3,4)};
//     cycle<tuple<int, int>, 2> myCycle = cycle<tuple<int, int>, 2>(vec);

//     std::cout << myCycle << std::endl;

    
// }


// #include <libsnark/common/utils.hpp>
// #include "../libsnark/libsnark/common/"

// using namespace libsnark;

int main()
{
    
    libff::mnt6_pp::init_public_params();
    typedef libff::mnt6_pp ppT;
    typedef libff::Fr<ppT> FieldT;
    typedef libsnark::sha256_two_to_one_hash_gadget<FieldT> HashT;
    typedef libsnark::knowledge_commitment<libff::G1<ppT>, libff::G2<ppT> > CommitT;

    const libff::G1<ppT> g1_1 = libff::G1<ppT>::G1_one;


    FieldT r_fr = FieldT::random_element();
    const libff::G1<ppT> g1_r1 = r_fr * libff::G1<ppT>::G1_one;      // 支持 域元素 和 bigint
    const libff::G1<ppT> g1_r2 = r_fr * r_fr.inverse() * libff::G1<ppT>::G1_one;      // 支持 域元素 和 bigint

    std::cout << "g1_1" << std::endl;   g1_1.print();
    std::cout << "g1_r1" << std::endl;  g1_r1.print();
    std::cout << "g1_r2" << std::endl;  g1_r2.print();
    

    libff::bit_vector elt_bits = libff::convert_field_element_to_bit_vector<FieldT>(r_fr);
    for (auto bit: elt_bits) {
        std::cout << bit;
    }
    std::cout << std::endl;

    const libff::G2<ppT> g2_1 = libff::G2<ppT>::G2_one;
    const libff::G2<ppT> g2_r1 = r_fr * g2_1;
    const libff::G2<ppT> g2_r2 = r_fr * r_fr.inverse() * g2_1;
    std::cout << "g2_1" << std::endl;   g2_1.print();
    std::cout << "g2_r1" << std::endl;  g2_r1.print();
    std::cout << "g2_r2" << std::endl;  g2_r2.print();

    // // 计算哈希值 h
    // const libff::bit_vector h = HashT::get_hash(libff::convert_field_element_to_bit_vector<FieldT>(r_fr));
    const libff::bit_vector h = HashT::get_hash(libff::bit_vector(512, false)); // 需要比特位为 512 位

    // // 生成承诺
    const CommitT kc(g1_r1, g2_r1);
    kc.print();

    FieldT coef = FieldT::random_element();
    std:: cout << (coef * kc == CommitT(coef * g1_r1, coef * g2_r1)) << std::endl;
    // // 输出承诺信息
    // std::cout << "Commitment: " << kc << std::endl;

    return 0;
}

