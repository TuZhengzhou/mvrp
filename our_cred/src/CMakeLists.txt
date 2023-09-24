set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

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
)

target_link_libraries(
  cred 

  ff
  gmpxx
)

add_executable(
  test_naive 

  test_naive.cpp
  structs.cpp
  prover.cpp
  verifier.cpp
  ipa.cpp
)

target_link_libraries(
  test_naive 

  ff
  gmpxx
)