set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DCRED_DEBUG)
endif()

find_package(OpenSSL REQUIRED)

include_directories(/usr/include)
include_directories(/usr/include/x86_64-linux-gnu)
include_directories(../../libsnark)
include_directories(../../libsnark/depends/libff)
include_directories(../includes)
include_directories(../src)

link_directories(/usr/local/lib)
link_directories(/usr/lib)

add_executable(
  cred 

  main.cpp
  structs.cpp
  prover.cpp
  verifier.cpp
  ipa.cpp
  kzg.cpp
  circuit_generator.cpp
)

target_link_libraries(
  cred 

  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_naive 

  test_naive.cpp
  structs.cpp
  prover.cpp
  verifier.cpp
  ipa.cpp
  kzg.cpp
  circuit_generator.cpp
)

target_link_libraries(
  test_naive 

  ff
  gmpxx
  OpenSSL::Crypto
)