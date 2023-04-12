set(CMAKE_CXX_STANDARD 14)

# add_definitions(-DCURVE=BN128)
find_package(OpenSSL REQUIRED)

# include_directories(/usr/local/include)
# include_directories(../libsnark)
# include_directories(../libsnark/depends/libfqfft/)
include_directories(../libsnark/depends/libff/)
include_directories(./)

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
    toy_test
    
    toy_test.cpp
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
