set(CMAKE_CXX_STANDARD 14)

add_definitions(-DCURVE=BN128)

include_directories(/usr/local/include)
include_directories(../libsnark)
include_directories(../libsnark/depends/libfqfft/)
link_directories(/usr/local/lib)

add_executable(qap main.cpp)

target_link_libraries(qap ff)
target_link_libraries(qap gmpxx)
target_link_libraries(qap gmp)
target_link_libraries(qap procps)
