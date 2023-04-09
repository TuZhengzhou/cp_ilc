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


add_test(NAME example1
        COMMAND example1)
set_tests_properties(example1 PROPERTIES PASS_REGULAR_EXPRESSION "1 4 5 \n6 7 8 \n9 10 11 \n")

add_test(NAME example2
        COMMAND example2)
set_tests_properties(example2 PROPERTIES PASS_REGULAR_EXPRESSION "1 4 5 \n6 7 8 \n9 10 11 \n")
