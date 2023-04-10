project(test)

find_package(OpenSSL REQUIRED)

include_directories(..)
include_directories(../../libsnark/depends/libff)
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

add_dependencies(bootle_check example1)
add_dependencies(bootle_check example2)
add_dependencies(bootle_check structs_test)