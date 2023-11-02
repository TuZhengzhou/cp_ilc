project(test)

find_package(OpenSSL REQUIRED)

include_directories(..)
include_directories(../../depends/libsnark/depends/libff)
link_directories(/usr/local/lib)

add_executable(example1 main.cpp)
target_link_libraries(
    example1
    
    ff      # libff 库生成的 lib 的 lib 名
    gmpxx
    OpenSSL::Crypto
)

add_executable(example2 main.cpp)
target_link_libraries(
    example2
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(structs_test structs_test.cpp)
target_link_libraries(
    structs_test
    
    ff
    gmpxx
    OpenSSL::Crypto
)


add_executable(test_eq test_eq.cpp)
target_link_libraries(
    test_eq
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_sum test_sum.cpp)
target_link_libraries(
    test_sum
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_prod test_prod.cpp)
target_link_libraries(
    test_prod
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_shift test_shift.cpp)
target_link_libraries(
    test_shift
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_same_prod test_same_prod.cpp)
target_link_libraries(
    test_same_prod
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_permutation test_permutation.cpp)
target_link_libraries(
    test_permutation
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(test_comEq test_comEq.cpp)
target_link_libraries(
    test_comEq
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_test(NAME example1
        COMMAND example1)
set_tests_properties(example1 PROPERTIES PASS_REGULAR_EXPRESSION "1 4 5 \n6 7 8 \n9 10 11 \n")

add_test(NAME example2
        COMMAND example2)
set_tests_properties(example2 PROPERTIES PASS_REGULAR_EXPRESSION "1 4 5 \n6 7 8 \n9 10 11 \n")

add_test(NAME structs_test
        COMMAND structs_test)
set_tests_properties(structs_test PROPERTIES PASS_REGULAR_EXPRESSION 
"2 2 2 
2 2 2 
2 2 2 
1 2 4 
8 16 32 
64 128 256 
2 4 8 
16 32 64 
128 256 512 
"
)



add_test(NAME test_eq
        COMMAND test_eq)

add_test(NAME test_sum
        COMMAND test_sum)

add_test(NAME test_prod
        COMMAND test_prod)

add_test(NAME test_shift
        COMMAND test_shift)

add_test(NAME test_same_prod
        COMMAND test_same_prod)

add_test(NAME test_permutation
        COMMAND test_permutation)

add_test(NAME test_comEq
        COMMAND test_comEq)

add_dependencies(bootle_check example1)
add_dependencies(bootle_check example2)
add_dependencies(bootle_check structs_test)
add_dependencies(bootle_check test_comEq)
add_dependencies(bootle_check test_eq)
add_dependencies(bootle_check test_sum)
# add_dependencies(bootle_check test_sum)
# add_dependencies(bootle_check test_)
# add_dependencies(bootle_check test_sum)
# add_dependencies(bootle_check test_sum)