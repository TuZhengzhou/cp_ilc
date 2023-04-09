set(CMAKE_CXX_STANDARD 14)

# add_definitions(-DCURVE=BN128)
find_package(OpenSSL REQUIRED)
include_directories(/usr/local/include)
include_directories(../libsnark)
include_directories(../libsnark/depends/libfqfft/)
include_directories(../libsnark/depends/libff/)
link_directories(/usr/local/lib)

add_executable(bootle17 
    main.cpp)

target_link_libraries(bootle17 ff)
target_link_libraries(bootle17 gmpxx)
target_link_libraries(bootle17 gmp)
target_link_libraries(bootle17 procps)
target_link_libraries(bootle17 OpenSSL::Crypto)

include_directories(/usr/local/include)
include_directories(../libsnark)
include_directories(../libsnark/depends/libfqfft/)
link_directories(/usr/local/lib)

add_executable(bootle_test test.cpp)

target_link_libraries(bootle17 ff)
target_link_libraries(bootle17 gmpxx)
target_link_libraries(bootle17 gmp)
target_link_libraries(bootle17 procps)
