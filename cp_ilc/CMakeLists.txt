set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_BUILD_TYPE "Release")


# add_definitions(-DCURVE=BN128)
find_package(OpenSSL REQUIRED)

# include_directories(/usr/local/include)
include_directories(../depends/libsnark)
# include_directories(../libsnark/depends/libfqfft/)
include_directories(../depends/libsnark/depends/libff/)
include_directories(./)
include_directories(./src/)

link_directories(/usr/local/lib)

add_executable(
    bootle17

    main.cpp
)
target_link_libraries(
    bootle17 
    
    ff
    gmpxx
    # gmp
    # procps
    OpenSSL::Crypto
)

add_executable(
    cred

    cred.cpp
)
target_link_libraries(
    cred 
    
    ff
    gmpxx
    # gmp
    # procps
    OpenSSL::Crypto
)

add_executable(
    toy_test
    
    toy_test.cpp
    ./src/eq.tcc
    ./src/sum.tcc
    ./src/same_prod.tcc
)
target_link_libraries(
    toy_test 
    
    ff
    gmpxx
    OpenSSL::Crypto
)

add_custom_target(bootle_check COMMAND ${CMAKE_CTEST_COMMAND})

option(BUILD_TEST "Whether or not to build the tests" ON)
if (${BUILD_TEST})
    message(STATUS "Building tests")
    enable_testing()  # Enable testing only works in root scope
    add_subdirectory(tests)
endif ()
