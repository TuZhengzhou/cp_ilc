language: cpp
sudo: required
dist: trusty
compiler:
  - gcc
os:
  - linux

addons:
  apt:
    sources:
      - ubuntu-toolchain-r-test
    packages:
      - gcc-5
      - g++-5

before_install:
  - sudo apt-get install build-essential cmake git libboost-all-dev cmake libgmp3-dev libssl-dev libprocps3-dev pkg-config gnuplot-x11 python-markdown

before_script:
  - git submodule init && git submodule update
  - mkdir build
  - cd build
  - cmake ..

script:
  - make
  - make check

notifications:
  email: false

报错：用 make 编译可执行文件时报错 error: #error C++ versions less than C++14 are not supported.
解决：前一步骤 cmake 时，不执行 cmake ..，而是执行 cmake -DCMAKE_CXX_STANDARD=14 ..
参考：https://blog.csdn.net/hit0803107/article/details/118441830

参数：EXCLUDE_FROM_ALL
作用：cmake 的 add_library, add_executable, add_subdirectory 等命令都有一个 EXCLUDE_FROM_ALL 参数.
根据cmake官网的解释，如果某个 target 或 subdirectory 被设置为 EXCLUDE_FROM_ALL 属性,那么这个 target (或这个 subdirectory 中的所有 target )就会被排除在 all target 列表之外，这样，当执行默认的 make(或nmake) 时，这个 target(或这个subdirectory中的所有target) 就不会被编译
需要时：手动 make [target] 编译 [target]
例如：make gadgetlib1_simple_test